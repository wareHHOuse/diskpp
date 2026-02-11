/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <cstddef>
#include <iostream>
#include <regex>
#include <set>
#include <string>
#include <filesystem>

#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/loaders/loader_gmsh.hpp"
#include "diskpp/loaders/loader_utils.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"
#include "rasdd.hpp"
#include "common.hpp"
#include "solvers.hpp"
#include "diskpp_git_revision.h"


template<typename FineMesh>
auto make_projectors(const FineMesh& fmsh, const coarse_mesh_t<FineMesh>& cmsh,
    cc2ff_t<FineMesh>& cc2ff, const std::vector<int>& occs, size_t coarse_degree, size_t fine_degree)
{
    using T = typename FineMesh::coordinate_type;
    using fine_mesh_type = FineMesh;
    using coarse_mesh_type = coarse_mesh_t<FineMesh>;
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using SM = Eigen::SparseMatrix<T>;
    using triplet = Eigen::Triplet<T>;

    auto cbs = disk::scalar_basis_size(coarse_degree, coarse_mesh_type::dimension);
    auto fbs = disk::scalar_basis_size(fine_degree, fine_mesh_type::dimension-1);

    auto nrows = fbs * fmsh.faces_size();
    auto ncols = cbs * cmsh.cells_size();

    /* Map from a DG space to an HHO space with zero Dirichlet BCs*/
    std::vector<std::optional<size_t>> cmap(fmsh.faces_size());
    size_t ci = 0;
    for (size_t i = 0; i < cmap.size(); i++) {
        auto ffc = fmsh.face_at(i);
        auto bi = fmsh.boundary_info(ffc);
        if ( bi.is_boundary() and not bi.is_internal() )
            continue;
        
        cmap[i] = ci++;
    }

    SM gproj = SM(ci*fbs, ncols);

    std::ofstream ofs("face_coeffs.m");

    std::vector<triplet> triplets;
    for (auto& ccl : cmsh) {
        /* Coarse cell basis */
        auto cphi = disk::basis::scaled_monomial_basis(cmsh, ccl, coarse_degree);

        auto col = cbs * offset(cmsh, ccl);

        const auto& ffcs = cc2ff.at(ccl);
        for (const auto& ffc : ffcs) {

            auto bi = fmsh.boundary_info(ffc);
            if ( bi.is_boundary() and not bi.is_internal() )
                continue;

            /* Fine face basis */
            auto fphi = disk::basis::scaled_monomial_basis(fmsh, ffc, fine_degree);

            auto fcofs_full = offset(fmsh, ffc);
            auto fcofs = cmap[fcofs_full].value();
            auto row = fbs * fcofs;

            auto C2F = integrate(fmsh, ffc, cphi, fphi);
            auto mass = integrate(fmsh, ffc, fphi, fphi);
            decltype(C2F) P = mass.ldlt().solve(C2F)/occs[fcofs_full];


            /*
            auto pts = points(fmsh, ffc);
            
            ofs << "line([" << pts[0].x() << "," << pts[1].x() << "], [" << pts[0].y() << "," << pts[1].y() << "], ";
            if (occs[fcofs_full] == 2)
                ofs << "'r' );" << std::endl;
            else
                ofs << "'k' );" << std::endl;
            */

            
            for (size_t i = 0; i < fbs; i++) {
                for (size_t j = 0; j < cbs; j++) {
                    triplets.push_back( {int(row+i), int(col+j), P(i,j)} );
                }
            }
        }
    }
    gproj.setFromTriplets( triplets.begin(), triplets.end() );
    //std::ofstream ofs("proj.txt");
    //for (auto& t : triplets)
    //    ofs << t.row()+1 << " " << t.col()+1 << " " << t.value() << std::endl;

    triplets.clear();

    return gproj;
}











enum class hho_variant {
    equal_order,
    mixed_order
};

enum class ras_mode {
    bicgstab,
    iterate
};

enum class mesh_type {
    triangles,
    cartesian,
    hexas,
    tetras,
    undefined,
};
struct solver_config {
    size_t          overlap;
    size_t          ras_maxiter;
    size_t          degree;
    bool            ras_debug;
    hho_variant     variant;
    ras_mode        mode;
    std::string     fn_bicg_hist;
    std::string     fn_silo;
    std::string     fn_err_hist;
    std::string     outdir;
    mesh_type       meshtype;
    size_t          imesh_reflevels;
    size_t          imesh_partitions;
    bool            use_twolevel;
    size_t          dgdegree;
};

struct iterdata {
    double error;
    double residual;
};

std::ostream&
operator<<(std::ostream& os, const iterdata& id)
{
    os << id.error << " " << id.residual;
    return os;
};

template<typename T>
void make_submesh_from_subdomain(const disk::simplicial_mesh<T,2>& msh,
    disk::simplicial_mesh<T,2>& submsh, size_t subdom_id,
    std::vector<size_t>& cell_sub_to_orig,
    std::vector<size_t>& face_sub_to_orig
    )
{
    auto stor = msh.backend_storage();
    auto substor = submsh.backend_storage();
    
    using cell_type = disk::simplicial_mesh<T,2>::cell_type;
    using face_type = disk::simplicial_mesh<T,2>::face_type;
    using node_type = disk::simplicial_mesh<T,2>::node_type;

    struct face_with_bndid {
        face_type                   fc;
        disk::boundary_descriptor   bndid;
        size_t                      orig_num;
    
        bool operator<(const face_with_bndid& other) {
            return fc < other.fc;
        }

        bool operator==(const face_with_bndid& other) {
            return fc == other.fc;
        }
    };

    std::vector<face_with_bndid> fwb;
    fwb.reserve( 3 * msh.cells_size() );
    size_t max_index = 0;
    for (size_t cli = 0; cli < msh.cells_size(); cli++) {
        const auto& cl = msh.cell_at(cli);
        /* cycle on all the elements and copy those of
         * the specified subdomain. */
        auto di = msh.domain_info(cl);
        if (di.tag() != subdom_id)
            continue;

        /* remember the original numbering of this element */    
        cell_sub_to_orig.push_back(cli);

        auto ptids = cl.point_ids();
        assert(ptids.size() == 3);
        for (auto& ptid : ptids) {
            max_index = std::max(ptid, max_index);
        }

        substor->surfaces.push_back(cl);

        /* the faces can be deduced, but we copy them because
         * we want to preserve boundary information */
        auto fcs = faces(msh, cl);
        for (auto fc : fcs) {
            auto bi = msh.boundary_info(fc);
            auto orig_num = offset(msh, fc);
            fwb.push_back( {fc, bi, orig_num} );
        }
    }

    /* mark the used indices, to take only the
     * necessary points */
    std::vector<std::optional<size_t>> used(max_index+1);
    for (auto& cl : substor->surfaces) {
        auto ptids = cl.point_ids();
        used[ptids[0]] = 1;
        used[ptids[1]] = 1;
        used[ptids[2]] = 1;
    }

    /* Compute the new indices */
    int ci = 0;
    for (auto& uopt : used) {
        if (uopt) {
            uopt = ci++;
        }
    }

    /* Copy only necessary points */
    substor->points.resize(ci);
    for (size_t i = 0; i < used.size(); i++) {
        const auto& uopt = used[i];
        if (uopt) {
            substor->points[uopt.value()] = stor->points[i];
        }
    }

    for (auto& cl : substor->surfaces) {
        auto ptids = cl.point_ids();
        assert(ptids.size() == 3);
        /* the optionals must be valid */
        auto np0 = used[ ptids[0] ].value();
        auto np1 = used[ ptids[1] ].value();
        auto np2 = used[ ptids[2] ].value();
        cl = cell_type({np0, np1, np2});
        //substor->edges.push_back( face_type{np0, np1} );
        //substor->edges.push_back( face_type{np1, np2} );
        //substor->edges.push_back( face_type{np0, np2} );
        substor->nodes.push_back( node_type{ np0 } );
        substor->nodes.push_back( node_type{ np1 } );
        substor->nodes.push_back( node_type{ np2 } );
    }

    /* Renumber faces */
    for (auto& [fc, bi, on] : fwb) {
        auto ptids = fc.point_ids();
        assert(ptids.size() == 2);
        /* the optionals must be valid */
        auto np0 = used[ ptids[0] ].value();
        auto np1 = used[ ptids[1] ].value();
        fc = face_type{np0, np1};
    }

    disk::priv::sort_uniq(fwb);
    substor->edges.resize(fwb.size());
    substor->boundary_info.resize(fwb.size());
    face_sub_to_orig.resize(fwb.size());

    for (size_t i = 0; i < fwb.size(); i++) {
        auto& [fc, bi, orig_num] = fwb[i];
        substor->edges[i] = fc;
        substor->boundary_info[i] = bi;
        face_sub_to_orig[i] = orig_num;
    }

    disk::priv::sort_uniq(substor->nodes);
    disk::priv::sort_uniq(substor->surfaces);
    substor->subdomain_info.resize(substor->surfaces.size(),
        disk::subdomain_descriptor(subdom_id) );
    detect_boundary(msh);
}

template<disk::mesh_2D Mesh>
void
diffusion_solver_refinement(const Mesh& cmsh, const solver_config& scfg)
{
    using mesh_type = Mesh;
    using scalar_type = typename mesh_type::coordinate_type;

    disk::simplicial_mesh<scalar_type, 2> fmsh;
    disk::submesh_via_gmsh(cmsh, fmsh, disk::average_diameter(cmsh)/4.0);

    std::vector<double> subdom_ids;
    for (auto& cl : fmsh) {
        auto di = fmsh.domain_info(cl);
        subdom_ids.push_back( di.tag() );
    }

    disk::silo_database silo;
    silo.create("c2f.silo");
    silo.add_mesh(cmsh, "coarse");
    silo.add_mesh(fmsh, "fine");
    silo.add_variable("fine", "domain_id", subdom_ids, disk::zonal_variable_t);

    
    for (auto& cl : cmsh) {
        std::cout << cl << std::endl;
        auto di = cmsh.domain_info(cl);
        auto fcs = faces(cmsh, cl);
        for (auto& fc : fcs) {
            std::cout << "  bndC: " << offset(cmsh, fc) << std::endl;
        }
    }

    for (size_t i = 0; i < cmsh.cells_size(); i++) {
        disk::simplicial_mesh<scalar_type, 2> smsh;
        std::vector<size_t> c_s2o, f_s2o;
        make_submesh_from_subdomain(fmsh, smsh, i+1, c_s2o, f_s2o);
        silo.add_mesh(smsh, "element_" + std::to_string(i));

        for (auto& cl : smsh) {
            auto fcs = faces(smsh, cl);
            for (auto& fc : fcs) {
                auto bi = smsh.boundary_info(fc);
                if (bi.is_boundary()) {
                    std::cout << "  bndS: " << bi.tag()-1 << std::endl;
                }
            }
        }
    }
    
}

template<typename T>
void
diffusion_solver_agglomeration(const disk::simplicial_mesh<T,2>& fmsh, const solver_config& scfg)
{
    using mesh_type = disk::simplicial_mesh<T,2>;
    using scalar_type = typename mesh_type::coordinate_type;

    disk::generic_mesh<T,2> cmsh;
    agglomerate_by_subdomain(fmsh, cmsh);

    std::vector<double> subdom_ids;
    for (auto& cl : fmsh) {
        auto di = fmsh.domain_info(cl);
        subdom_ids.push_back( di.tag() );
    }

    disk::silo_database silo;
    silo.create("f2c.silo");
    silo.add_mesh(cmsh, "coarse");
    silo.add_mesh(fmsh, "fine");
    silo.add_variable("fine", "domain_id", subdom_ids, disk::zonal_variable_t);
}

//template<typename Mesh>
template<typename T>
void
diffusion_solver(const disk::simplicial_mesh<T,2>& msh, const solver_config& scfg)
{
    using Mesh = disk::simplicial_mesh<T,2>;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using hho_space = hho_space<Mesh>;
    //using T = typename hho_space::scalar_type;
    using cell_basis_type = typename hho_space::cell_basis_type;
    using face_basis_type = typename hho_space::face_basis_type;

    bool mixed_order = (scfg.variant == hho_variant::mixed_order);

    auto face_degree = scfg.degree;
    auto cell_degree = mixed_order ? scfg.degree + 1 : scfg.degree;

    degree_info di(cell_degree, face_degree);

    disk::source_functor<Mesh> f;

    auto assm = make_assembler(msh, di);

    auto sizeF = face_basis_type::size_of_degree(di.face);

    using dm = disk::dynamic_matrix<T>;
    using dv = disk::dynamic_vector<T>;
    std::vector<std::pair<dm, dv>> local_contribs;
    local_contribs.reserve(msh.cells_size());

    /* MATRIX ASSEMBLY */
    std::cout << "Assembling global matrix..." << std::endl;
    for (auto& cl : msh)
    {
        dm lhs;
        if (mixed_order) {
            auto [R, A] = local_operator(msh, cl, di);
            auto S = local_stabilization_hdg(msh, cl, di);
            lhs = A+S;
        } else {
            auto [R, A] = local_operator(msh, cl, di);
            auto S = local_stabilization(msh, cl, di, R);
            lhs = A+S;
        }

        auto phiT = hho_space::cell_basis(msh, cl, di.cell);
        dv rhs = integrate(msh, cl, f, phiT);
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
        local_contribs.push_back({lhs, rhs});
    }
    assm.finalize();

    std::cout << "Unknowns: " << assm.LHS.rows() << " ";
    std::cout << "Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    /* PREPARE RAS */
    rasdd ras(msh, sizeF, scfg.overlap);
    ras.prepare( assm.LHS, assm.dirichlet_faces_flags() );
    if (scfg.ras_debug)
        ras.save_debug_data();

    /* PREPARE TWO-LEVEL */
    disk::generic_mesh<T, 2> coarsemsh;
    disk::sparse_matrix<T> proj;
    if (scfg.use_twolevel) {
        agglomerate_by_subdomain(msh, coarsemsh);
        auto cc2ff = make_cc2ff(msh, coarsemsh);
        proj = make_projectors(msh, coarsemsh, cc2ff.first, cc2ff.second, scfg.dgdegree, scfg.degree);
    }

    /* HELPER FUNCTION FOR POSTPRO */
    auto postpro = [&](const std::string& filename,
        const dv& sol, const dv& r, const dv& rr) {
        std::vector<T> u_data;
        std::vector<T> r_data;
        std::vector<T> rr_data;
        T error = 0.0;
        disk::solution_functor<Mesh> u_sol;
        for(size_t cell_i = 0; cell_i < msh.cells_size(); cell_i++)
        {
            auto cl = msh.cell_at(cell_i);   
            auto& [lhs, rhs] = local_contribs[cell_i];
            auto phiT = hho_space::cell_basis(msh, cl, di.cell);
            dv sol_ana = local_reduction(msh, cl, di, u_sol);
            auto locsolF = assm.take_local_solution(msh, cl, sol);
            dv locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            u_data.push_back(locsol(0));
            dv diff = locsol - sol_ana;
            error += diff.dot(lhs*diff);

            /* This does not make much sense, as the residual is not
             * in the solution space. It is however useful since it
             * allows to check that things are zero where they should
             * be zero. */
            auto locresF = assm.take_local_solution(msh, cl, r);
            dv zero = dv::Zero(rhs.size());
            dv locres = disk::hho::deschur(lhs, zero, locresF, phiT);
            r_data.push_back(locres(0));

            auto locresF2 = assm.take_local_solution(msh, cl, rr);
            dv locres2 = disk::hho::deschur(lhs, zero, locresF2, phiT);
            rr_data.push_back(locres2(0));
        }
        if (filename != "") {
            disk::silo_database silo_db;
            silo_db.create(filename);
            silo_db.add_mesh(msh, "mesh");
            silo_db.add_variable("mesh", "u", u_data, disk::zonal_variable_t);
            silo_db.add_variable("mesh", "residual", r_data, disk::zonal_variable_t);
            silo_db.add_variable("mesh", "residual2", rr_data, disk::zonal_variable_t);
        }
        return std::sqrt(error);
    };

    /* SOLVE */
    if (scfg.mode == ras_mode::iterate) {
        std::vector<iterdata> ids;
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        disk::sparse_matrix<T> Ac;
        mumps_solver<T> solverAc;
        if (scfg.use_twolevel) {
            Ac = proj.transpose() * (assm.LHS * proj);
            solverAc.factorize(Ac);
        }

        mumps_solver<T> solverLHS;
        solverLHS.factorize(assm.LHS);
        for (size_t niter = 0; niter < scfg.ras_maxiter; niter++) {
            std::cout << "RAS iteration " << niter+1 << std::endl;
            iterdata id;
            disk::dynamic_vector<T> r = assm.RHS - assm.LHS*sol;
            disk::dynamic_vector<T> rasr = ras(r);
            sol = sol + rasr;
            disk::dynamic_vector<T> e2 = disk::dynamic_vector<T>::Zero(r.size());
            if (scfg.use_twolevel) {
                r = assm.RHS - assm.LHS*sol;
                disk::dynamic_vector<T> rc = proj.transpose() * r;
                disk::dynamic_vector<T> ec = solverAc.solve(rc);
                e2 = proj * ec;
                sol = sol + e2;
            }

            disk::dynamic_vector<T> e = solverLHS.solve(r);
            
            //e2 = proj * (proj.transpose() * e);

            
            //disk::dynamic_vector<T> lc = proj.transpose() * assm.RHS;
            //disk::dynamic_vector<T> sol_uc = solverAc.solve( lc );
            
            //sol = proj * sol_uc;
            
            
            
            id.residual = r.norm();
            std::cout << "Postprocessing..." << std::endl;
            std::stringstream ss;
            
            if (scfg.fn_silo != "") {
                ss << scfg.fn_silo << "_iter" << niter << ".silo";
                std::cout << "  Saving solution to " << ss.str() << std::endl;
            }
            
            id.error = postpro(ss.str(), sol, e, e2);

            std::cout << "  A-norm error: " << id.error << ", residual norm: " << id.residual << std::endl;
            ids.push_back(id);
        }

        std::cout << "  Saving error history to " << scfg.fn_err_hist << std::endl;
        std::ofstream ofs(scfg.fn_err_hist);
        for (auto& id : ids)
            ofs << id << std::endl;
    }
    else {
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        bicgstab_io<T> bio;
        bio.verbose = true;
        sol = bicgstab(assm.LHS, assm.RHS, ras, bio);
        std::cout << "Postprocessing..." << std::endl;
        std::cout << "  Saving BiCGStab history to " << scfg.fn_bicg_hist << std::endl;
        std::ofstream ofs(scfg.fn_bicg_hist);
        for (auto& rr : bio.history)
            ofs << rr << std::endl;

        if (scfg.fn_silo != "") {
            std::cout << "  Saving solution to " << scfg.fn_silo << std::endl;
        }
        
        auto err = postpro(scfg.fn_silo, sol, bio.residual, bio.residual);
        std::cout << "  A-norm error: " << err << std::endl;
    }
}


void
run_on_internal_mesh(const solver_config& scfg)
{
    using T = double;

    if (scfg.meshtype == mesh_type::triangles) {
        disk::simplicial_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver_refinement(msh, scfg);
    }

    if (scfg.meshtype == mesh_type::cartesian) {
        disk::cartesian_mesh<T,2> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver_refinement(msh, scfg);
    }

    if (scfg.meshtype == mesh_type::hexas) {
        disk::generic_mesh<T,2> msh;
        auto mesher = make_fvca5_hex_mesher(msh);
        mesher.make_level(scfg.imesh_reflevels);

        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver_refinement(msh, scfg);
    }

    #if 0
    if (scfg.meshtype == mesh_type::tetras) {
        disk::simplicial_mesh<T,3> msh;
        auto mesher = make_simple_mesher(msh);
        for (size_t l = 0; l < scfg.imesh_reflevels; l++)
            mesher.refine();
        
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver_refinement(msh, scfg);
    }
    #endif
}


int main(int argc, char **argv)
{
    std::cout << "HHO-DD solver, revision " << GIT_REVISION << std::endl;
    
    rusage_monitor rm(false);

    solver_config scfg;
    scfg.overlap = 1;
    scfg.ras_maxiter = 20;
    scfg.degree = 0;
    scfg.ras_debug = false;
    scfg.variant = hho_variant::equal_order;
    scfg.mode = ras_mode::bicgstab;
    scfg.fn_bicg_hist = "bicgstab.txt";
    scfg.fn_err_hist = "error.txt";
    scfg.fn_silo = "";
    scfg.outdir = "";
    scfg.meshtype = mesh_type::undefined;
    scfg.imesh_partitions = 2;
    scfg.imesh_reflevels = 2;
    scfg.use_twolevel = false;
    scfg.dgdegree = 1;

    int opt;
    while ((opt = getopt(argc, argv, "b:d:e:i:k:m:o:p:q:r:s:DIMR2")) != -1) {
        switch (opt) {

        case 'b': /* BiCGStab history filename (bicgstab mode) */
            scfg.fn_bicg_hist = optarg;
            break;

        case 'd': /* output data directory */
            scfg.outdir = optarg;
            break;

        case 'e': /* Error history filename (iter mode) */
            scfg.fn_err_hist = optarg;
            break;

        case 'i': {
            int optval = std::max(0, atoi(optarg));
            scfg.ras_maxiter = optval;
            } break;
        
        case 'k': { /* hho method degree */
            int optval = std::max(0, atoi(optarg));
            scfg.degree = optval;
            } break;

        case 'm': { /* internal mesh type */
            if ( std::string(optarg) == "tri" )
                scfg.meshtype = mesh_type::triangles;
            else if ( std::string(optarg) == "quad" )
                scfg.meshtype = mesh_type::cartesian;
            else if ( std::string(optarg) == "hex" )
                scfg.meshtype = mesh_type::hexas;
            else if ( std::string(optarg) == "tet" )
                scfg.meshtype = mesh_type::tetras;
            else
                std::cout << "Warning: Wrong mesh type '" << optarg << "'" << std::endl;
            } break;

        case 'o': { /* number of overlap layers */
            int optval = std::max(0, atoi(optarg));
            scfg.overlap = optval;
            } break;

        case 'p': { /* partitions for internal mesh */
            int optval = std::max(0, atoi(optarg));
            scfg.imesh_partitions = optval;
            } break;

        case 'q': { /* coarse space degree */
            int optval = std::max(0, atoi(optarg));
            scfg.dgdegree = optval;
            } break;

        case 'r': { /* refinement levels for internal mesh */
            int optval = std::max(0, atoi(optarg));
            scfg.imesh_reflevels = optval;
            } break;
            
        case 's': /* silo filename or prefix */
            scfg.fn_silo = optarg;
            break;

        case 'D': /* Output subdomain debugging data */
            scfg.ras_debug = true;
            break;
        
        case 'I': /* Iterate instead of using BiCGStab */
            scfg.mode = ras_mode::iterate;
            break;

        case 'M': /* enable mixed order */
            scfg.variant = hho_variant::mixed_order;
            break;

        case 'R': /* Enable resource usage reporting */
            rm.enabled(true);
            break;

        case '2':
            scfg.use_twolevel = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;

    if ( (argc < 1) and (scfg.meshtype == mesh_type::undefined) ) {
        std::cout << "Usage: hhodd [options] [filename]\n\n";
        std::cout << " If a mesh is provided via [filename] is, the currently supported formats are:\n";
        std::cout << "   * GMSH .geo (not .msh): if a simplicial 2D mesh is wanted the GMSH script\n";
        std::cout << "     must have extension `.geo2s`, whereas if a 3D simplicial mesh is wanted\n";
        std::cout << "     the script must have `.geo3s` extension\n";
        std::cout << "   * FVCA5-style meshes in .typ1 format\n";
        std::cout << "   * FVCA6-style meshes in .msh format\n\n";
        std::cout << " Meshes can be also generated by the internal meshers, see -m below. In any\n";
        std::cout << " case the code assumes that the domain is the unit hypercube [0,1]^d.\n\n";
        std::cout << " Other options:\n";
        std::cout << "   -b <str> : BiCGStab history filename (default: bicgstab.txt)\n";
        std::cout << "   -d <str> : Directory for all the output files\n";
        std::cout << "   -e <str> : Error history filename (default: error.txt)\n";
        std::cout << "   -i <int> : RAS iterations (in iter mode)\n";
        std::cout << "   -k <int> : HHO degree\n";
        std::cout << "   -m <std> : Use internal meshers (tri, quad, hex or tet)\n";
        std::cout << "   -o <int> : Number of overlap layers\n";
        std::cout << "   -p <int> : Partitions for internal mesh\n";
        std::cout << "   -q <int> : DG degree\n";
        std::cout << "   -r <int> : Refinement levels for internal mesh\n";
        std::cout << "   -s <str> : Silo filename (bicgstab mode) or prefix (iter mode)\n";
        std::cout << "   -D       : Output subdomain debugging data\n";
        std::cout << "   -I       : Iterate instead of using BiCGStab\n";
        std::cout << "   -M       : Enable mixed-order HHO\n";
        std::cout << "   -R       : Report resource usage at exit\n";
        std::cout << "   -2       : Use two-level method\n";
        std::cout << "   -K <int> : Coarse space degree\n";
        return 1;
    }

    if (scfg.outdir != "") {
        std::filesystem::create_directory(scfg.outdir);
        scfg.fn_bicg_hist = scfg.outdir + "/" + scfg.fn_bicg_hist;
        scfg.fn_err_hist = scfg.outdir + "/" + scfg.fn_err_hist;
        if (scfg.fn_silo != "")
            scfg.fn_silo = scfg.outdir + "/" + scfg.fn_silo;
    }

    const char *mesh_filename = argv[0];
    if ( mesh_filename and (scfg.meshtype != mesh_type::undefined) ) {
        std::cout << "Warning: both '-m' and a filename were specified. Ignoring -m.\n";
    }

    using T = double;

    if (scfg.meshtype != mesh_type::undefined)
    {
        run_on_internal_mesh(scfg);
        return 0;
    }

    assert(argc > 0);

    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,2>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        rebase_subdomain_numbering(msh);
        diffusion_solver(msh, scfg);

        disk::silo_database db;
        db.create("plain_hho.silo");
        db.add_mesh(msh, "srcmesh");
        hho_diffusion_solver(msh, scfg.degree, db);
        return 0;
    }

#if 0
    /* GMSH */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial" << std::endl;

        using mesh_type = disk::simplicial_mesh<T,3>;
        mesh_type msh;
        disk::gmsh_geometry_loader< mesh_type > loader;
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        diffusion_solver(msh, scfg);
        return 0;
    }

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T,2> msh;
        disk::load_mesh_fvca5_2d<T>(mesh_filename, msh);
        //partition_unit_square_mesh(msh, scfg.imesh_partitions);
        disk::simplicial_mesh<T,2> fmsh;
        submesh_via_gmsh(msh, fmsh, 0.1);
        diffusion_solver(fmsh, scfg);
        return 0;
    }

    /* FVCA6 3D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::load_mesh_fvca6_3d<T>(mesh_filename, msh);
        partition_unit_square_mesh(msh, scfg.imesh_partitions);
        diffusion_solver(msh, scfg);
        return 0;
    }
#endif
    std::cout << "Unrecognized mesh format. Exiting." << std::endl;
    return 1;
}
