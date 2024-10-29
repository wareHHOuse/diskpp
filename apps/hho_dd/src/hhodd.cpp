#include <iostream>
#include <regex>
#include <set>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"

template<typename T>
disk::dynamic_vector<T>
bicgstab(const Eigen::SparseMatrix<T>& A, const disk::dynamic_vector<T>& b,
    std::map<size_t, std::pair<Eigen::SparseMatrix<T>, Eigen::SparseMatrix<T>>>& Rmats,
    std::map<size_t, Eigen::SparseLU<Eigen::SparseMatrix<T>>>& invAs)
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

    auto apply_precond = [&Rmats, &invAs](const dv& v) {
        dv ret = dv::Zero(v.size());
        for (auto& [tag, mats] : Rmats) {
            auto& [Rj, Rtj] = mats;
            dv w1 = invAs[tag].solve(Rj*v);
            ret = ret + Rtj.transpose() * w1;
        }
        return ret;
    };

    for (size_t i = 0; i < 50; i++) {
        std::cout << "BiCGStab " << i << " " << nr << std::endl;
        dv y = apply_precond(p);
        dv v = A*y;
        T alpha = rho/rhat.dot(v);
        dv h = x + alpha*y;
        dv s = r - alpha*v;
        //if( h suff accurate ) {
        //}
        dv z = apply_precond(s);
        dv t = A*z;
        dv Mt = apply_precond(t);
        T omega = Mt.dot(z)/Mt.dot(Mt);
        x = h + omega*z;
        r = s - omega*t;
        nr = r.norm();
        if ( (nr/nr0) < 1e-7 ) {
            return x;
        }

        T rho_prev = rho;
        rho = rhat.dot(r);
        T beta = (rho/rho_prev)*(alpha/omega);
        p = r + beta*(p - omega*v);
    }


    return x;
}

namespace disk {


template<typename Mesh>
struct source;

template<mesh_2D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return 2.0*M_PI*M_PI*sx*sy;
    }
};


template<mesh_3D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return 3.0*M_PI*M_PI*sx*sy*sz;
    }
};


template<typename Mesh>
struct solution;

template<mesh_2D Mesh>
struct solution<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return sx*sy;
    }
};


template<mesh_3D Mesh>
struct solution<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return sx*sy*sz;
    }
};


}

template<typename Mesh>
using sdmap_t = std::map<size_t, std::vector<typename Mesh::cell_type>>;

using flagmap_t = std::map<size_t, std::vector<bool>>;

template<typename Mesh>
static flagmap_t
make_cell_flagmap(const Mesh& msh) {
    flagmap_t ret;
    for (size_t clnum = 0; clnum < msh.cells_size(); clnum++) {
        auto cl = msh.cell_at(clnum);
        auto di = msh.domain_info(cl);
        auto& flags = ret[di.tag()];
        if (flags.size() == 0)
            flags.resize( msh.cells_size() );
        flags[clnum] = true; 
    }

    return ret;
}

template<typename Mesh>
static flagmap_t
make_face_flagmap(const Mesh& msh) {
    flagmap_t ret;
    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        auto& flags = ret[di.tag()];
        if (flags.size() == 0)
            flags.resize( msh.faces_size() );
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs) {
            auto bi = msh.boundary_info(fc);
            if (bi.is_boundary() /*and bi.is_internal()*/)
                continue;
            flags[offset(msh, fc)] = true;
        } 
    }

    return ret;
}

template<typename Mesh>
auto
make_overlapping_subdomains(const Mesh& msh, size_t overlap_layers)
{
    using cell_type = typename Mesh::cell_type;
    using cvi = typename std::vector<cell_type>::iterator;
    auto cvf = connectivity_via_faces(msh);

    auto subdomain_cells = make_cell_flagmap(msh);
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

    auto subdomain_faces = make_face_flagmap(msh);
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

    return std::pair(subdomain_cells, subdomain_faces);
}

template<typename Mesh>
void
hhodd(const Mesh& msh, size_t levels)
{
    auto [sd_cells, sd_faces] = make_overlapping_subdomains(msh, levels);
    
    disk::silo_database silo;
    silo.create("hhodd.silo");
    silo.add_mesh(msh, "mesh");

    for (auto& [tag, cell_present] : sd_cells)
    {
        std::vector<double> yesno(msh.cells_size());
        std::transform(cell_present.begin(), cell_present.end(),
            yesno.begin(), [](bool x) { return double(x); } );

        std::stringstream ss;
        ss << "domain" << tag;
        silo.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
    }

    silo.close();

    if constexpr (Mesh::dimension == 2) {
        for (auto& [tag, ifcs] : sd_faces ) {
            std::stringstream ss;
            ss << "faces_" << tag << ".m";
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

template<typename Mesh>
auto domfaces(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;
    std::map<size_t, disk::dynamic_vector<T>> ret;

    auto cvf = connectivity_via_faces(msh);

    for (auto& cl : msh) {
        auto di = msh.domain_info(cl);
        auto tag = di.tag();
        auto& vec = ret[tag];
        if (vec.size() == 0)
            vec = disk::dynamic_vector<T>::Zero(msh.faces_size());
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

    return ret;
}

template<typename Mesh>
void
diffusion_solver(const Mesh& msh, size_t degree, size_t levels, size_t iters)
{
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;
    using face_basis_type = typename hho_space<mesh_type>::face_basis_type;

    degree_info di(degree, degree);

    disk::source<Mesh> f;

    auto assm = make_assembler(msh, di);

    auto sizeF = face_basis_type::size_of_degree(di.face);

    std::cout << "Assembling global matrix..." << std::endl;
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
    }
    assm.finalize();

    std::cout << "Unknowns: " << assm.LHS.rows() << " ";
    std::cout << "Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    auto dff = assm.dirichlet_faces_flags();
    std::vector<bool> ndff(dff.size());
    std::transform(dff.begin(), dff.end(), ndff.begin(),
        [](bool d) {return not d;});

    auto num_non_dirichlet = std::count(ndff.begin(), ndff.end(), true);

    auto [Rj_cells_bytag, Rj_faces_bytag] = make_overlapping_subdomains(msh, levels);
    auto Rtj_faces = domfaces(msh);

    std::map<size_t, std::pair<Eigen::SparseMatrix<T>, Eigen::SparseMatrix<T>>> Rmats;

    for (auto& [tag, Rj_faces] : Rj_faces_bytag) {    
        auto num_fcs_Rj = std::count(Rj_faces.begin(), Rj_faces.end(), true);
        std::cout << tag << ": " << num_fcs_Rj << std::endl;
        std::vector<Eigen::Triplet<T>> Rj_trip;
        Eigen::SparseMatrix<T> Rj(num_fcs_Rj*sizeF, num_non_dirichlet*sizeF);
        std::vector<Eigen::Triplet<T>> Rtj_trip;
        Eigen::SparseMatrix<T> Rtj(num_fcs_Rj*sizeF, num_non_dirichlet*sizeF);

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
                if (Rtj_faces[tag][i] > 0.0)
                    for (size_t k = 0; k < sizeF; k++)
                        Rtj_trip.push_back({sizeF*row+k, sizeF*col+k, Rtj_faces[tag][i]});

                row++;
            }
            col++;
        }

        Rj.setFromTriplets(Rj_trip.begin(), Rj_trip.end());
        Rtj.setFromTriplets(Rtj_trip.begin(), Rtj_trip.end());
        
        Rmats[tag] = std::pair(Rj, Rtj);
    }

    std::cout << "Factorizing local matrices..." << std::endl;

    std::map<size_t, Eigen::SparseLU<Eigen::SparseMatrix<T>>> invAs;
    for (auto& [tag, mats] : Rmats) {
        std::cout << " - Subdomain " << tag << std::endl;
        auto& [Rj, Rtj] = mats;
        Eigen::SparseMatrix<T> Aj = Rj*assm.LHS*Rj.transpose();
        invAs[tag].compute(Aj);
    }

    /*
    std::vector<double> errors;
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
    for (size_t niter = 0; niter < iters; niter++) {
        disk::dynamic_vector<T> w = assm.RHS - assm.LHS*sol;
        std::cout << "Iteration " << niter << std::endl;
        for (auto& [tag, mats] : Rmats) {
            auto& [Rj, Rtj] = mats;
            disk::dynamic_vector<T> w1 = invAs[tag].solve(Rj*w);
            sol = sol + Rtj.transpose() * w1;
        }

        std::vector<T> u_data;

        T error = 0.0;
        disk::solution<Mesh> u_sol;

        std::cout << "Postprocessing..." << std::endl;
        for (auto& cl : msh)
        {
            auto [R, A] = local_operator(msh, cl, di);
            auto S = local_stabilization(msh, cl, di, R);       


            disk::dynamic_matrix<T> lhs = A+S;

            auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
            disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
        
            disk::dynamic_vector<T> sol_ana = local_reduction(msh, cl, di, u_sol);

            auto locsolF = assm.take_local_solution(msh, cl, sol);

            disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            u_data.push_back(locsol(0));

            disk::dynamic_vector<T> diff = locsol - sol_ana;
            error += diff.dot(lhs*diff);
        }
        errors.push_back(std::sqrt(error));
        std::cout << "A-norm error: " << std::sqrt(error) << std::endl;

        std::stringstream ss;
        ss << "hhodd_iter_" << niter << ".silo";

        disk::silo_database silo_db;
        silo_db.create(ss.str());
        silo_db.add_mesh(msh, "mesh");

        silo_db.add_variable("mesh", "u", u_data, disk::zonal_variable_t);

        for (auto& [tag, cell_present] : Rj_cells_bytag)
        {
            std::vector<double> yesno(msh.cells_size());
            std::transform(cell_present.begin(), cell_present.end(),
                yesno.begin(), [](bool x) { return double(x); } );

            std::stringstream ss;
            ss << "domain" << tag;
            silo_db.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
        }
    }

    std::ofstream ofs("error.txt");
    for (auto& e : errors)
        ofs << e << std::endl;
    */

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
    sol = bicgstab(assm.LHS, assm.RHS, Rmats, invAs);

    std::vector<T> u_data;

    T error = 0.0;
    disk::solution<Mesh> u_sol;

    std::cout << "Postprocessing..." << std::endl;
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);       


        disk::dynamic_matrix<T> lhs = A+S;

        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
    
        disk::dynamic_vector<T> sol_ana = local_reduction(msh, cl, di, u_sol);

        auto locsolF = assm.take_local_solution(msh, cl, sol);

        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        u_data.push_back(locsol(0));

        disk::dynamic_vector<T> diff = locsol - sol_ana;
        error += diff.dot(lhs*diff);
    }

    std::cout << "A-norm error: " << std::sqrt(error) << std::endl;
}

int main(int argc, char **argv)
{
    int levels = 1;
    int degree = 0;
    int iters = 5;
    int opt;
    while ((opt = getopt(argc, argv, "l:k:i:")) != -1) {
        switch (opt) {
        case 'l': {
            int l = atoi(optarg);
            if (l < 1) {
                std::cerr << "Levels must be positive, resetting to 1." << std::endl;
                l = 1;
            }
            levels = l;
            } break;
        
        case 'k': {
            int k = atoi(optarg);
            if (k < 0) {
                std::cerr << "Degree must be greater than 0, resetting to 0." << std::endl;
                k = 0;
            }
            degree = k;
            } break;
        case 'i': {
            int i = atoi(optarg);
            if (i < 1) {
                std::cerr << "Iterations must be positive resetting to 1." << std::endl;
                i = 1;
            }
            iters = i;
            } break;
        }
    }
    argc -= optind;
    argv += optind;

    if (argc < 1) {
        std::cout << "missing filename" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[0];



    using T = double;

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,2>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        diffusion_solver(msh, degree, levels, iters);

        return 0;
    }

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,3>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        diffusion_solver(msh, degree, levels, iters);

        return 0;
    }

}
