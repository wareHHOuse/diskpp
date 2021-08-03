/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020,2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */

/* Implementation of "A mixed Hybrid High-Order formulation for
 * linear interior transmission elliptic problems"
 * by R. Bustinza & J. Munguia La-Cotera
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "solvers/mumps.hpp"

template<typename Mesh>
std::set<size_t>
subdomain_tags(const Mesh& msh)
{
    std::set<size_t> tags;

    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        tags.insert( di.tag() );
    }

    return tags;
}

template<typename Mesh>
size_t
num_interface_faces(Mesh& msh)
{
    size_t nif = 0;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        auto& fc = *itor;
        auto bi = msh.boundary_info(fc);
        if (bi.is_internal())
            nif++;
    }

    return nif;
}

template<typename Mesh>
typename Mesh::coordinate_type
rhs_term(const Mesh&, const typename Mesh::point_type& pt, size_t subdom)
{
    return 1.0;
}

template<typename Mesh>
void test(const Mesh& msh)
{
    auto subdom_tags = subdomain_tags(msh);
    if (subdom_tags.size() != 2)
    {
        std::cout << "The mesh must have exactly two subdomains" << std::endl;
        return;
    }

    size_t internal_tag = 1;

    size_t num_int_fcs = num_interface_faces(msh);
    std::cout << num_int_fcs << std::endl;

    using T = typename Mesh::coordinate_type;

    size_t degree = 1;
    disk::hho_degree_info hdi(degree);

    size_t fbs = disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension);
    size_t num_dofs = (msh.faces_size() + 2*num_int_fcs) * fbs;

    size_t ibase = msh.faces_size() * fbs;
    size_t iofs = 0;


    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);

        auto rhs_fun = [&](const typename Mesh::point_type& pt) -> auto {
            return rhs_term(msh, pt, di.tag());
        };

        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        disk::dynamic_matrix<T> A = gr.second + stab;
        auto sc     = make_scalar_static_condensation(msh, cl, hdi, A, rhs);

        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            for (size_t j = 0; j < fcs.size(); j++)
            {

            }
        }


            #if 0
            auto bi = msh.boundary_info(fc);
            if (bi.is_internal())
            {
                /* We are on the interface */
                if (di.tag() == internal_tag)
                {
                    /* We are in Ω₁ (the internal domain) */
                    auto fb = disk::make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                    disk::dynamic_matrix<T> mass = disk::make_mass_matrix(msh, fc, fb);
                    iofs += fbs;
                }
                else
                {
                    /* We are in Ω₂ (the external domain) */
                    auto fb = disk::make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                    disk::dynamic_matrix<T> mass = disk::make_mass_matrix(msh, fc, fb);
                }
            }
            else
            {

            }
            #endif
    }
}

int main(int argc, const char **argv)
{
    using T = double;

    if (argc != 2)
        return 1;
    
    const char *mesh_filename = argv[1];

    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh);

        disk::silo_database silo;
        silo.create("test.silo");
        silo.add_mesh(msh, "mesh");
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh);

        disk::silo_database silo;
        silo.create("test.silo");
        silo.add_mesh(msh, "mesh");
    }

    /* GMSH 3D generic */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D generic" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh);

        disk::silo_database silo;
        silo.create("test.silo");
        silo.add_mesh(msh, "mesh");
    }

    return 0;
}