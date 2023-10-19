#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"

#include "diskpp/methods/hho"

#include "mumps.hpp"

template<typename Mesh>
struct source;

template<disk::mesh_2D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        //if (pt.x() < 0.7)
        //    return (5.0*M_PI*M_PI/4.0)*std::sin(M_PI*pt.x()*0.5)*std::sin(M_PI*pt.y());
        //else
        //    return (10.0*M_PI*M_PI/4.0)*std::sin(M_PI*pt.x()*0.5)*std::sin(M_PI*pt.y());;
        
        if (pt.x() >= 1.0 and pt.x() <= 2.0 and pt.y() >= 1.0 and pt.y() <= 2.0)
            return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return 0.0;
    }
};


template<disk::mesh_3D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return 0.0;
    }
};


template<typename Mesh>
struct dirichlet_jump;

template<disk::mesh_2D Mesh>
struct dirichlet_jump<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return 0.0;
    }
};


template<disk::mesh_3D Mesh>
struct dirichlet_jump<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return 0.0;
    }
};

template<typename Mesh>
struct neumann_jump;

template<disk::mesh_2D Mesh>
struct neumann_jump<Mesh> {
    using point_type = typename Mesh::point_type;
    //using normal_type = static_vector<typename Mesh::coordinate_type, 2>;
    auto operator()(const point_type& pt) const {
        return 0.0;
    }
};


template<disk::mesh_3D Mesh>
struct neumann_jump<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return 0.0;
    }
};



template<typename Mesh>
void lt_solver(Mesh& msh, size_t degree)
{
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;

    using T = typename hho_space<Mesh>::scalar_type;
    using Space = hho_space<Mesh>;
    using cbt = typename hho_space<mesh_type>::cell_basis_type;
    using fbt = typename hho_space<mesh_type>::face_basis_type;

    degree_info di(degree);
    auto [szT_sb, szF_sb, szR_sb] = space_dimensions<Space>(di);
    /* Workaround for a Clang bug BEGIN*/
    auto szT = szT_sb;
    auto szF = szF_sb;
    auto szR = szR_sb;
    /* Workaround for a Clang bug END*/

    source<Mesh> f;
    dirichlet_jump<Mesh> gD;
    neumann_jump<Mesh> gN;

    auto assm = make_assembler(msh, di);

    timecounter tc;
    tc.tic();

    auto compute_local_contrib = [&](const mesh_type& msh, const cell_type& cl) {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);   

        auto subdomain_id = msh.domain_info(cl).tag();

        auto dc = 1.0;

        if (subdomain_id == 2)
            dc = 10;

        disk::dynamic_matrix<T> lhs = dc*(A+S);
    
        auto phiT = cbt(msh, cl, di.cell);

        size_t ofs = szT;
        auto fcs = faces(msh, cl);
        auto total_dofs = szT + fcs.size()*szF;
        
        disk::dynamic_vector<T> rhs = disk::dynamic_vector<T>::Zero(total_dofs);
        rhs.head(szT) = integrate(msh, cl, f, phiT);

        disk::dynamic_vector<T> gD_rhs = disk::dynamic_vector<T>::Zero(total_dofs);
        disk::dynamic_vector<T> gN_rhs = disk::dynamic_vector<T>::Zero(total_dofs);

        for (auto& fc : fcs) {

            auto boundary_id = msh.boundary_info(fc).tag();

            if (subdomain_id == 2 and (boundary_id == 8)) {
                auto phiF = fbt(msh, fc, di.face);
                gD_rhs.segment(ofs, szF) += L2_project(msh, fc, gD, phiF);
                gN_rhs.segment(ofs, szF) += integrate(msh, fc, gN, phiF);
            }
            ofs += szF;
        }

        rhs += -lhs*gD_rhs + gN_rhs;

        return std::tuple(lhs, rhs, phiT);
    };


    /* Assembly loop */
    for (auto& cl : msh)
    {
        auto [lhs, rhs, phiT] = compute_local_contrib(msh, cl);
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
    }
    assm.finalize();

    std::cout << "Assembly time: " << tc.toc() << std::endl;
    std::cout << "Unknowns: " << assm.LHS.rows() << " ";
    std::cout << "Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    /* Solve */
    tc.tic();
    disk::dynamic_vector<T> sol = mumps_lu(assm.LHS, assm.RHS);
    std::cout << "Solver time: " << tc.toc() << std::endl;

    /* Postprocess */
    std::vector<T> u_data;
    tc.tic();
    for (auto& cl : msh)
    {
        auto [lhs, rhs, phiT] = compute_local_contrib(msh, cl);
        auto locsolF = assm.take_local_solution(msh, cl, sol);
        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        u_data.push_back(locsol(0));
    }
    std::cout << "Postpro time: " << tc.toc() << std::endl;

    disk::silo_database silo_db;
    silo_db.create("diffusion.silo");
    silo_db.add_mesh(msh, "mesh");

    silo_db.add_variable("mesh", "u", u_data, disk::zonal_variable_t);

}

int main(int argc, char **argv)
{
    if (argc < 2)
        return 1;
    
    const char *mesh_filename = argv[1];
    size_t degree = 1;
    using T = double;

#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        lt_solver(msh, degree);
    }
#else
    std::cout << "GMSH support not compiled. Exiting." << std::endl;
    return 1;
#endif

    return 0;
}