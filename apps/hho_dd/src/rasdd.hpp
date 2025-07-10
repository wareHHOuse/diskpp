/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include "diskpp/common/timecounter.hpp"
#include "diskpp/output/silo.hpp"
#include "mumps.hpp"

enum class bicgstab_status {
    undefined,
    converged,
    not_converged,
    diverged,
};

template<typename T>
struct bicgstab_io {
    bicgstab_status     status;
    size_t              iterations;
    T                   rr_min;
    T                   rr_max;
    T                   rr;
    bool                verbose;
    std::vector<T>      history;
    disk::dynamic_vector<T> residual;

    bicgstab_io() : status(bicgstab_status::undefined),
        iterations(1000), rr_min(1e-6), rr_max(1e8),
        verbose(false)
    {}
};

template<typename T, typename Functor>
disk::dynamic_vector<T>
bicgstab(const Eigen::SparseMatrix<T>& A, const disk::dynamic_vector<T>& b,
    const Functor& precond, bicgstab_io<T>& bio)
{
    using dv = disk::dynamic_vector<T>;
    dv ret = dv::Zero( b.size() );

    dv x = ret;
    dv r = b - A*x;
    dv rhat = b;
    T rho = rhat.dot(r);
    dv p = r;

    T nr, nr0;
    nr = nr0 = r.norm();
    bio.rr = nr/nr0;

    for (size_t i = 0; i < bio.iterations; i++) {
        dv y = precond(p);
        dv v = A*y;
        T alpha = rho/rhat.dot(v);
        dv h = x + alpha*y;
        dv s = r - alpha*v;
        dv z = precond(s);
        dv t = A*z;
        dv Mt = precond(t);
        T omega = Mt.dot(z)/Mt.dot(Mt);
        x = h + omega*z;
        r = s - omega*t;

        nr = r.norm();
        bio.rr = nr/nr0;
        if ( bio.rr < bio.rr_min ) {
            if (bio.verbose) {
                std::cout << "BiCGStab converged with RR = ";
                std::cout << bio.rr << std::endl;
            }
            bio.iterations = i+1;
            bio.residual = r;
            bio.status = bicgstab_status::converged;
            return x;
        }
        if ( bio.rr > bio.rr_max ) {
            if (bio.verbose) {
                std::cout << "BiCGStab DIVERGED with RR = ";
                std::cout << bio.rr << std::endl;
            }
            bio.iterations = i+1;
            bio.residual = r;
            bio.status = bicgstab_status::diverged;
            return x;
        }

        T rho_prev = rho;
        rho = rhat.dot(r);
        T beta = (rho/rho_prev)*(alpha/omega);
        p = r + beta*(p - omega*v);
        bio.history.push_back(bio.rr);
        if (bio.verbose) {
            std::cout << "BiCGStab iteration " << i+1 << "/";
            std::cout << bio.iterations << ", RR = ";
            std::cout << bio.rr << std::endl;
        }
    }

    bio.status = bicgstab_status::not_converged;
    return x;
}

template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class rasdd {
    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using CoordT = typename Mesh::coordinate_type;
    using spmat_t = Eigen::SparseMatrix<ScalT>;
    using flagmap_type = std::map<size_t, std::vector<bool>>;
    using weightmap_type = std::map<size_t, disk::dynamic_vector<CoordT>>;
    using Rmatmap_type = std::map<size_t, std::pair<spmat_t, spmat_t>>;

#ifdef HHODD_USE_EIGEN_SPARSELU
    using lu_type = Eigen::SparseLU<spmat_t>;
#endif
#ifdef HHODD_USE_PARDISO
    using lu_type = Eigen::PardisoLU<spmat_t>;
#endif
#ifdef HHODD_USE_MUMPS
    using lu_type = mumps_solver<ScalT>;
#endif

    using lumap_type = std::map<size_t, lu_type>;

    flagmap_type    subdomain_cells;
    flagmap_type    subdomain_faces;
    weightmap_type  sub_to_global;
    lumap_type      invAs;

    Rmatmap_type    Rmats;

    
    const mesh_type&    msh;
    size_t              sizeF;

    disk::neighbour_connectivity<mesh_type>   cvf;

    void make_flagmaps(void)
    {
        subdomain_cells.clear();
        subdomain_faces.clear();
        /* for each subdomain make a bitmap indicating if
         * a {cell,face} belongs to that subdomain */
        for (size_t clnum = 0; clnum < msh.cells_size(); clnum++) 
        {
            auto cl = msh.cell_at(clnum);
            auto di = msh.domain_info(cl);
            auto& cell_flags = subdomain_cells[di.tag()];
            if (cell_flags.size() == 0)
                cell_flags.resize( msh.cells_size() );
            cell_flags[clnum] = true; 

            auto& face_flags = subdomain_faces[di.tag()];
            if (face_flags.size() == 0)
                face_flags.resize( msh.faces_size() );
            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
                auto bi = msh.boundary_info(fc);
                if (bi.is_boundary() /*and bi.is_internal()*/)
                    continue;
                face_flags[offset(msh, fc)] = true;
            }
        }
    }

    void make_overlapping_subdomains(size_t overlap_layers)
    {
        make_flagmaps();

        for (size_t ol = 0; ol < overlap_layers; ol++) {
            std::map<size_t, std::set<size_t>> layer_cells;
            /* Iterate on the currently identified subdomains */
            for (auto& [tag, present] : subdomain_cells) {
                /* for each cell */
                for (size_t clnum = 0; clnum < present.size(); clnum++) {
                    if (not present[clnum])
                        continue;

                    auto cl = msh.cell_at(clnum);
                    auto fcs = faces(msh, cl);
                    for (auto& fc : fcs) {
                        /* for each neighbour */
                        auto neigh = cvf.neighbour_via(msh, cl, fc);
                        if (not neigh)
                            continue;

                        auto neighnum = offset(msh, neigh.value());
                        if (not present[neighnum])
                            layer_cells[tag].insert(neighnum);
                    }
                }
            }

            /* move the found neighbours to the subdomain */
            for (auto& [tag, cellnums] : layer_cells) {
                for (auto& cellnum : cellnums)
                    subdomain_cells[tag][cellnum] = true;
            }
        }

        for (auto& [tag, cell_present] : subdomain_cells) {
            for(size_t clnum = 0; clnum < cell_present.size(); clnum++) {
                if (not cell_present[clnum])
                    continue;
                auto cl = msh.cell_at(clnum);
                auto fcs = faces(msh, cl);
                for (auto& fc : fcs) {
                    auto neigh = cvf.neighbour_via(msh, cl, fc);
                    if (not neigh)
                        continue;

                    auto neigh_id = offset(msh, neigh.value());
                    if (cell_present[neigh_id] /*or include_dd_bnd*/) {
                        subdomain_faces[tag][offset(msh, fc)] = true;
                    }
                }
            }
        }
    }

    void make_sub_to_global_weights(void)
    {
        sub_to_global.clear();
        for (auto& cl : msh) {
            auto di = msh.domain_info(cl);
            auto tag = di.tag();
            auto& vec = sub_to_global[tag];
            if (vec.size() == 0)
                vec = disk::dynamic_vector<ScalT>::Zero(msh.faces_size());
            auto fcs = faces(msh, cl);
            for (auto& fc : fcs) {
                auto fcnum = offset(msh, fc);
                auto bi = msh.boundary_info(fc);
                if (not bi.is_boundary())
                    vec[fcnum] = 1.0;

                if (bi.is_boundary() and bi.is_internal())
                    vec[fcnum] = 0.5;
            }
        }
    }

public:

    void prepare(const spmat_t& A, const std::vector<bool>& dirichlet_faces)
    {
        using vec_of_triplets_t = std::vector<Eigen::Triplet<ScalT>>;

        auto dff = dirichlet_faces;
        std::vector<bool> ndff(dff.size());
        std::transform(dff.begin(), dff.end(), ndff.begin(),
            [](bool d) {return not d;});

        auto num_non_dirichlet = std::count(ndff.begin(), ndff.end(), true);

        for (auto& [tag, Rj_faces] : subdomain_faces) {    
            auto num_fcs_Rj = std::count(Rj_faces.begin(), Rj_faces.end(), true);
            std::cout << tag << ": " << num_fcs_Rj << std::endl;
            vec_of_triplets_t Rj_trip;
            spmat_t Rj(num_fcs_Rj*sizeF, num_non_dirichlet*sizeF);
            vec_of_triplets_t Rtj_trip;
            spmat_t Rtj(num_fcs_Rj*sizeF, num_non_dirichlet*sizeF);

            size_t row = 0;
            size_t col = 0;
            for (size_t i = 0; i < ndff.size(); i++) {
                if (not ndff[i])
                    continue;

                if (Rj_faces[i]) {
                    assert(row < num_fcs_Rj);
                    assert(col < num_non_dirichlet);
                    for (size_t k = 0; k < sizeF; k++)
                        Rj_trip.push_back({sizeF*row+k, sizeF*col+k, 1.0});
                    if (sub_to_global[tag][i] > 0.0)
                        for (size_t k = 0; k < sizeF; k++)
                            Rtj_trip.push_back({sizeF*row+k, sizeF*col+k, sub_to_global[tag][i]});

                    row++;
                }
                col++;
            }

            Rj.setFromTriplets(Rj_trip.begin(), Rj_trip.end());
            Rtj.setFromTriplets(Rtj_trip.begin(), Rtj_trip.end());
            
            Rmats[tag] = std::pair(Rj, Rtj);
        }

#ifdef HHODD_USE_EIGEN_SPARSELU
        std::cout << "Factorizing local matrices (SparseLU)..." << std::endl;
#endif
#ifdef HHODD_USE_PARDISO
        std::cout << "Factorizing local matrices (PARDISO)..." << std::endl;
#endif
#ifdef HHODD_USE_MUMPS
        std::cout << "Factorizing local matrices (MUMPS)..." << std::endl;
#endif
        for (auto& [tag, mats] : Rmats) {
            timecounter tc;
            tc.tic();
            std::cout << " - Subdomain " << tag << "..." << std::flush;
            auto& [Rj, Rtj] = mats;
            spmat_t Aj = Rj*A*Rj.transpose();
            invAs[tag].compute(Aj);
            std::cout << "done in " << tc.toc() << " seconds" << std::endl;
        }
    }


    rasdd(const Mesh& pmsh, size_t psizeF, size_t overlap)
        : msh(pmsh), sizeF(psizeF), cvf(connectivity_via_faces(pmsh))
    {
        make_overlapping_subdomains(overlap);
        make_sub_to_global_weights();
    }

    disk::dynamic_vector<ScalT>
    operator()(const disk::dynamic_vector<ScalT>& v) const
    {
        using dv = disk::dynamic_vector<ScalT>;
        dv ret = dv::Zero(v.size());
        for (auto& [tag, mats] : Rmats) {
            auto& [Rj, Rtj] = mats;
            dv w0 = Rj*v;
            dv w1 = invAs.at(tag).solve(w0);
            ret = ret + Rtj.transpose() * w1;
        }
        return ret;
    }

    void save_debug_data(void) const {
        disk::silo_database silo_db;
        silo_db.create("ras_debug.silo");
        silo_db.add_mesh(msh, "mesh");

        for (auto& [tag, cell_present] : subdomain_cells)
        {
            std::vector<double> yesno(msh.cells_size());
            std::transform(cell_present.begin(), cell_present.end(),
                yesno.begin(), [](bool x) { return double(x); } );

            std::stringstream ss;
            ss << "ext_domain" << tag;
            silo_db.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
        }

        std::vector<double> initial_dom;
        initial_dom.reserve( msh.cells_size() );
        for (auto& cl : msh) {
            auto di = msh.domain_info(cl);
            initial_dom.push_back( di.tag() );
        }
        silo_db.add_variable("mesh", "orig_doms", initial_dom, disk::zonal_variable_t);

        if constexpr (Mesh::dimension == 2) {
            for (auto& [tag, ifcs] : subdomain_faces ) {
                std::stringstream ss;
                ss << "ras_debug_faces_" << tag << ".m";
                std::ofstream ofs(ss.str());
                for (size_t i = 0; i < ifcs.size(); i++) {
                    if (not ifcs[i])
                        continue;
                    auto fc = msh.face_at(i);
                    auto pts = points(msh, fc);
                    ofs << "line([" << pts[0].x() << "," << pts[1].x() << "],";
                    ofs << "[" << pts[0].y() << "," << pts[1].y() << "]);" << std::endl;
                }
            }
        }
    }
};

