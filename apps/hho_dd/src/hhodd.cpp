/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <regex>
#include <set>

#include "diskpp/common/util.h"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho"
#include "rasdd.hpp"



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


} // namespace disk



enum class hho_variant {
    equal_order,
    mixed_order
};

enum class ras_mode {
    bicgstab,
    iterate
};

struct solver_config {
    size_t          overlap;
    size_t          ras_maxiter;
    size_t          degree;
    bool            ras_debug;
    hho_variant     variant;
    ras_mode        mode;
};


template<typename Mesh>
void
diffusion_solver(const Mesh& msh, const solver_config& scfg)
{
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;
    using cell_basis_type = typename hho_space<mesh_type>::cell_basis_type;
    using face_basis_type = typename hho_space<mesh_type>::face_basis_type;

    bool mixed_order = (scfg.variant == hho_variant::mixed_order);

    auto face_degree = scfg.degree;
    auto cell_degree = mixed_order ? scfg.degree + 1 : scfg.degree;

    degree_info di(cell_degree, face_degree);

    disk::source<Mesh> f;

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

        auto phiT = cell_basis_type(msh, cl, di.cell);
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

    /* HELPER FUNCTION FOR POSTPRO */
    auto postpro = [&](const std::string& filename, const dv& sol) {
        std::vector<T> u_data;
        T error = 0.0;
        disk::solution<Mesh> u_sol;
        for(size_t cell_i = 0; cell_i < msh.cells_size(); cell_i++)
        {
            auto cl = msh.cell_at(cell_i);   
            auto& [lhs, rhs] = local_contribs[cell_i];
            auto phiT = cell_basis_type(msh, cl, di.cell);
            dv sol_ana = local_reduction(msh, cl, di, u_sol);
            auto locsolF = assm.take_local_solution(msh, cl, sol);
            dv locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            u_data.push_back(locsol(0));
            dv diff = locsol - sol_ana;
            error += diff.dot(lhs*diff);
        }
        disk::silo_database silo_db;
        silo_db.create(filename);
        silo_db.add_mesh(msh, "mesh");
        silo_db.add_variable("mesh", "u", u_data, disk::zonal_variable_t);
        return std::sqrt(error);
    };

    /* SOLVE */
    if (scfg.mode == ras_mode::iterate) {
        std::vector<double> errors;
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        for (size_t niter = 0; niter < scfg.ras_maxiter; niter++) {
            sol = sol + ras(assm.RHS - assm.LHS*sol);
            
            std::cout << "Postprocessing..." << std::endl;
            std::stringstream ss;
            ss << "ras_iter_" << niter << ".silo";
            std::cout << "  Saving solution to " << ss.str() << std::endl;
            auto err = postpro(ss.str(), sol);
            errors.push_back(err);
            std::cout << "  A-norm error: " << err << std::endl;
        }

        std::cout << "  Saving error history to error.txt" << std::endl;
        std::ofstream ofs("error.txt");
        for (auto& e : errors)
            ofs << e << std::endl;
    }
    else {
        disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
        bicgstab_io<T> bio;
        bio.verbose = true;
        sol = bicgstab(assm.LHS, assm.RHS, ras, bio);
        std::cout << "Postprocessing..." << std::endl;
        std::cout << "  Saving BiCGStab history to bicgstab.txt" << std::endl;
        std::ofstream ofs("bicgstab.txt");
        for (auto& rr : bio.history)
            ofs << rr << std::endl;
        std::cout << "  Saving solution to ras.silo" << std::endl;
        auto err = postpro("ras.silo", sol);
        std::cout << "  A-norm error: " << err << std::endl;
    }
}

int main(int argc, char **argv)
{
    rusage_monitor rm;

    solver_config scfg;
    scfg.overlap = 1;
    scfg.ras_maxiter = 20;
    scfg.degree = 0;
    scfg.ras_debug = false;
    scfg.variant = hho_variant::equal_order;
    scfg.mode = ras_mode::bicgstab;

    int opt;
    while ((opt = getopt(argc, argv, "o:i:k:MID")) != -1) {
        switch (opt) {

        case 'o': { /* number of overlap layers */
            int optval = std::max(0, atoi(optarg));
            scfg.overlap = optval;
            } break;

        case 'i': {
            int optval = std::max(0, atoi(optarg));
            scfg.ras_maxiter = optval;
            } break;
        
        case 'k': { /* hho method degree */
            int optval = std::max(0, atoi(optarg));
            scfg.degree = optval;
            } break;
        
        case 'M': /* enable mixed order */
            scfg.variant = hho_variant::mixed_order;
            break;

        case 'I': /* Iterate instead of using BiCGStab */
            scfg.mode = ras_mode::iterate;
            break;

        case 'D': /* Output subdomain debugging data */
            scfg.mode = ras_mode::iterate;
            break;
        }
    }
    argc -= optind;
    argv += optind;

    if (argc < 1) {
        std::cout << "Usage: hhodd [options] <filename>\n\n";
        std::cout << " <filename> is the name of a GMSH script with `.geo2s`\n";
        std::cout << "   extension if the mesh is 2D or with `.geo3s` extension\n";
        std::cout << "   if the mesh is 3D.\n\nOther options:\n";
        std::cout << "   -o <int> : number of overlap layers\n";
        std::cout << "   -i <int> : maximum ras iterations (in iter mode)\n";
        std::cout << "   -k <int> : HHO degree\n";
        std::cout << "   -M       : enable mixed-order HHO\n";
        std::cout << "   -I       : iterate instead of using BiCGStab\n";
        std::cout << "   -D       : output subdomain debugging data\n";
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
        disk::make_interpartition_boundaries(msh);

        diffusion_solver(msh, scfg);
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
        disk::make_interpartition_boundaries(msh);

        diffusion_solver(msh, scfg);
        return 0;
    }

    std::cout << "Unrecognized mesh format. Exiting." << std::endl;
    return 1;
}
