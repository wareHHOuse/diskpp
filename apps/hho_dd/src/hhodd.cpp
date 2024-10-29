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


}

template<typename Mesh>
using sdmap_t = std::map<size_t, std::vector<typename Mesh::cell_type>>;






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

    rasdd ras(msh, sizeF, levels);
    ras.prepare( assm.LHS, assm.dirichlet_faces_flags() );

#ifdef ITERATE
    std::vector<double> errors;
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
    for (size_t niter = 0; niter < iters; niter++) {
        sol = sol + ras(assm.RHS - assm.LHS*sol);

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
        /*
        for (auto& [tag, cell_present] : Rj_cells_bytag)
        {
            std::vector<double> yesno(msh.cells_size());
            std::transform(cell_present.begin(), cell_present.end(),
                yesno.begin(), [](bool x) { return double(x); } );

            std::stringstream ss;
            ss << "domain" << tag;
            silo_db.add_variable("mesh", ss.str(), yesno, disk::zonal_variable_t);
        }
        */
    }

    std::ofstream ofs("error.txt");
    for (auto& e : errors)
        ofs << e << std::endl;
    
#else
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.size());
    sol = bicgstab(assm.LHS, assm.RHS, ras);

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
#endif
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
