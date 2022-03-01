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
void test(Mesh& msh, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    auto fun = [&](const typename Mesh::point_type& pt) {
        return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    };

    disk::hho_degree_info hdi(degree);

    T err_sq = 0.0;
    T err_proj = 0.0;
    for (auto& cl : msh)
    {
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto rb = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
        auto gr = make_scalar_hho_laplacian(msh, cl, hdi);

        disk::dynamic_vector<T> projH = disk::project_function(msh, cl, hdi, fun, 1);
        disk::dynamic_vector<T> projR = disk::dynamic_vector<T>::Zero(rb.size());
        disk::dynamic_vector<T> avefix = disk::dynamic_vector<T>::Zero(rb.size());
        projR.tail(rb.size()-1) = gr.first*projH;

        avefix.head(cb.size()) = projH.head(cb.size());
        avefix.tail(rb.size()-1) -= projR.tail(rb.size()-1);

        T avg = 0.0;
        auto qps = integrate(msh, cl, 2*hdi.reconstruction_degree()+1);
        for (auto& qp : qps)
        {
            auto rphi = rb.eval_functions(qp.point());
            for (size_t i = 0; i < rb.size(); i++)
                avg += qp.weight()*avefix(i)*rphi(i);
        }
        avg /= measure(msh,cl);


        projR(0) = avg;//projH(0); 

        for (auto& qp : qps)
        {
            auto cphi = cb.eval_functions(qp.point());
            auto diffT = projH.head(cb.size()).dot(cphi) - fun(qp.point());
            err_proj += qp.weight() * std::pow(diffT, 2.);

            auto rphi = rb.eval_functions(qp.point());
            auto diffR = projR.dot(rphi) - fun(qp.point());
            err_sq += qp.weight() * std::pow(diffR, 2);
        }
    }

    std::cout << "Projection error = " << std::sqrt(err_proj) << std::endl;
    std::cout << "Reconstruction error = " << std::sqrt(err_sq) << std::endl;
}

int main(int argc, const char **argv)
{
    using T = double;

    if (argc != 3)
        return 1;
    
    const char *mesh_filename = argv[1];
    size_t degree = atoi(argv[2]);

    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh, degree);
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh, degree);
    }

    /* GMSH 3D generic */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D generic" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        test(msh, degree);
    }

    return 0;
}
