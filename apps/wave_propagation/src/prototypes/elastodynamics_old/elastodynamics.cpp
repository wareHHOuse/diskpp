
//  Contributions by Omar Durán and Romain Mottier

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <list>
#include <getopt.h>


#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
using namespace Eigen;

#include "timecounter.h"
#include "methods/hho"
#include "geometry/geometry.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "output/silo.hpp"

// application common sources
#include "../common/display_settings.hpp"
#include "../common/fitted_geometry_builders.hpp"
#include "../common/linear_solver.hpp"
#include "../common/vec_analytic_functions.hpp"
#include "../common/preprocessor.hpp"
#include "../common/postprocessor.hpp"

// implicit RK schemes
#include "../common/dirk_hho_scheme.hpp"
#include "../common/dirk_butcher_tableau.hpp"

// explicit RK schemes
#include "../common/ssprk_hho_scheme.hpp"
#include "../common/ssprk_shu_osher_tableau.hpp"

void HeterogeneousGar6more2DIHHOSecondOrder(int argc, char **argv);

void HeterogeneousGar6more2DIHHOFirstOrder(int argc, char **argv);

void HeterogeneousGar6more2DIHHOFirstOrderTwoFields(int argc, char **argv);

void Gar6more2DIHHOSecondOrder(int argc, char **argv);

void Gar6more2DIHHOFirstOrder(int argc, char **argv);

void HeterogeneousIHHOFirstOrder(int argc, char **argv);

void HeterogeneousIHHOSecondOrder(int argc, char **argv);

void EHHOFirstOrder(int argc, char **argv);

void IHHOFirstOrderTwoFields(int argc, char **argv);

void IHHOFirstOrder(int argc, char **argv);

void IHHOSecondOrder(int argc, char **argv);

void HHOOneFieldConvergenceExample(int argc, char **argv);

void HHOTwoFieldsConvergenceExample(int argc, char **argv);

void HHOThreeFieldsConvergenceExample(int argc, char **argv);

int main(int argc, char **argv)
{

    HeterogeneousGar6more2DIHHOFirstOrderTwoFields(argc, argv);
    
//    HeterogeneousGar6more2DIHHOFirstOrder(argc, argv);
//    HeterogeneousGar6more2DIHHOSecondOrder(argc, argv);

//    Gar6more2DIHHOFirstOrder(argc, argv);
//    Gar6more2DIHHOSecondOrder(argc, argv);
    
//    HeterogeneousIHHOFirstOrder(argc, argv);
//    HeterogeneousIHHOSecondOrder(argc, argv);
    
//    EHHOFirstOrder(argc, argv);
//    IHHOFirstOrderTwoFields(argc, argv);
//    IHHOFirstOrder(argc, argv);
//    IHHOSecondOrder(argc, argv);
    
//    // Examples using main app objects for solving a linear elastic problem with optimal convergence rates
    // Primal HHO
//    HHOOneFieldConvergenceExample(argc, argv);
    // Dual HHO
//    HHOTwoFieldsConvergenceExample(argc, argv);
//    HHOThreeFieldsConvergenceExample(argc, argv);
    
    return 0;
}

void HHOOneFieldConvergenceExample(int argc, char **argv){
    
    using RealType = double;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    
    simulation_data sim_data = preprocessor::process_convergence_test_args(argc, argv);
    sim_data.print_simulation_data();

    // Manufactured exact solution
    bool quadratic_function_Q = false;
    bool Nonzero_Dirichlet_Q = false;
    RealType lambda = 1000.0;
    auto exact_vec_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        RealType x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType ux = (1 - x)*x*(1 - y)*y;
            RealType uy = (1 - x)*x*(1 - y)*y;
            return static_vector<RealType, 2>{ux, uy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType ux = - std::sin(M_PI * pt.x()) * std::cos(M_PI * pt.y());
                RealType uy = + std::cos(M_PI * pt.x()) * std::sin(M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }else{
                RealType ux = std::sin(2.0 * M_PI * pt.x()) * std::sin(2.0 * M_PI * pt.y());
                RealType uy = std::sin(3.0 * M_PI * pt.x()) * std::sin(3.0 * M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }

        }
        
    };

    auto exact_flux_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            RealType sxx = 2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            RealType sxy = (1 - x)*x*(1 - y) - (1 - x)*x*y + (1 - x)*(1 - y)*y - x*(1 - y)*y;
            RealType syy = 2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            sigma(0,0) = sxx;
            sigma(0,1) = sxy;
            sigma(1,0) = sxy;
            sigma(1,1) = syy;
            return sigma;
        }else{
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            if (Nonzero_Dirichlet_Q) {
                RealType sxx = - 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                RealType syy = + 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                sigma(0,0) = sxx;
                sigma(1,1) = syy;
                return sigma;
            }else{
                RealType sxx = 4*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType sxy = 2*M_PI*std::cos(2*M_PI*y)*std::sin(2*M_PI*x)
                            + 3*M_PI*std::cos(3*M_PI*x)*std::sin(3*M_PI*y);
                RealType syy = 6*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                sigma(0,0) = sxx;
                sigma(0,1) = sxy;
                sigma(1,0) = sxy;
                sigma(1,1) = syy;
                return sigma;
            }
        }

    };

    auto rhs_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType fx = 2*(1 + x*x + y*(-5 + 3*y) + x*(-3 + 4*y));
            RealType fy = 2*(1 + 3*x*x + (-3 + y)*y + x*(-5 + 4*y));
            return static_vector<RealType, 2>{-fx, -fy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType fx = + 2.0 * M_PI * M_PI * ( std::sin(M_PI * x) * std::cos( M_PI * y));
                RealType fy = - 2.0 * M_PI * M_PI * ( std::cos(M_PI * x) * std::sin( M_PI * y));
                return static_vector<RealType, 2>{-fx, -fy};
            }else{
                RealType fx = M_PI*M_PI*(9*(1 + lambda)*std::cos(3*M_PI*x)*std::cos(3*M_PI*y) - 4*(3 + lambda)*std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType fy = M_PI*M_PI*(4*(1 + lambda)*std::cos(2*M_PI*x)*std::cos(2*M_PI*y) - 9*(3 + lambda)*std::sin(3*M_PI*x)*std::sin(3*M_PI*y));
                return static_vector<RealType, 2>{-fx, -fy};
            }
        }
    };

    // simple material
    RealType rho = 1.0;
    RealType vp;
    if (Nonzero_Dirichlet_Q || quadratic_function_Q) {
        vp = std::sqrt(3.0);
    }else{
        vp = std::sqrt(2.0 + lambda);
    }
    RealType vs = 1.0;
    elastic_material_data<RealType> material(rho,vp,vs);
    
    std::ofstream error_file("steady_vector_error.txt");
    
    for(size_t k = 1; k <= sim_data.m_k_degree; k++){
        std::cout << bold << cyan << "Running an approximation with k : " << k << reset << std::endl;
        error_file << "Approximation with k : " << k << std::endl;
        for(size_t l = 0; l <= sim_data.m_n_divs; l++){
            
            // Building a cartesian mesh
            timecounter tc;
            tc.tic();
            RealType lx = 1.0;
            RealType ly = 1.0;
            size_t nx = 2;
            size_t ny = 2;
            mesh_type msh;
            
            cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
            mesh_builder.refine_mesh(l);
            mesh_builder.build_mesh();
            mesh_builder.move_to_mesh_storage(msh);
            std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
            
            // Creating HHO approximation spaces and corresponding linear operator
            size_t cell_k_degree = k;
            if(sim_data.m_hdg_stabilization_Q){
                cell_k_degree++;
            }
            disk::hho_degree_info hho_di(cell_k_degree,k);

            // Solving a scalar primal HHO problem
            boundary_type bnd(msh);
            bnd.addDirichletEverywhere(exact_vec_fun);
            tc.tic();
            auto assembler = elastodynamic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
            if(sim_data.m_hdg_stabilization_Q){
                assembler.set_hdg_stabilization();
            }
            assembler.load_material_data(msh,material);
            assembler.assemble(msh, rhs_fun);
            assembler.apply_bc(msh);
            tc.toc();
            std::cout << bold << cyan << "Assemble in : " << tc << " seconds" << reset << std::endl;
            
            // Solving LS
            
            Matrix<RealType, Dynamic, 1> x_dof;
            if (sim_data.m_sc_Q) {
                tc.tic();
                linear_solver<RealType> analysis(assembler.LHS,assembler.get_n_face_dof());
                analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations (SC) : " << analysis.n_equations() << reset << std::endl;
            }else{
                tc.tic();
                linear_solver<RealType> analysis(assembler.LHS);
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations : " << analysis.n_equations() << reset << std::endl;
            }
            
            // Computing errors
            postprocessor<mesh_type>::compute_errors_one_field_vectorial(msh, hho_di, assembler, x_dof, exact_vec_fun, exact_flux_fun,error_file);
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "steady_vector_k" + std::to_string(k) + "_";
                postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, l, msh, hho_di, x_dof, exact_vec_fun, false);
            }
        }
        error_file << std::endl << std::endl;
    }
    error_file.close();
}

void HHOTwoFieldsConvergenceExample(int argc, char **argv){
    
    using RealType = double;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    
    simulation_data sim_data = preprocessor::process_convergence_test_args(argc, argv);
    sim_data.print_simulation_data();

    // Manufactured exact solution
    bool quadratic_function_Q = false;
    bool Nonzero_Dirichlet_Q = false;
    RealType lambda = 1000.0;
    auto exact_vec_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        RealType x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType ux = (1 - x)*x*(1 - y)*y;
            RealType uy = (1 - x)*x*(1 - y)*y;
            return static_vector<RealType, 2>{ux, uy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType ux = - std::sin(M_PI * pt.x()) * std::cos(M_PI * pt.y());
                RealType uy = + std::cos(M_PI * pt.x()) * std::sin(M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }else{
                RealType ux = std::sin(2.0 * M_PI * pt.x()) * std::sin(2.0 * M_PI * pt.y());
                RealType uy = std::sin(3.0 * M_PI * pt.x()) * std::sin(3.0 * M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }

        }
        
    };

    auto exact_flux_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            RealType sxx = 2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            RealType sxy = (1 - x)*x*(1 - y) - (1 - x)*x*y + (1 - x)*(1 - y)*y - x*(1 - y)*y;
            RealType syy = 2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            sigma(0,0) = sxx;
            sigma(0,1) = sxy;
            sigma(1,0) = sxy;
            sigma(1,1) = syy;
            return sigma;
        }else{
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            if (Nonzero_Dirichlet_Q) {
                RealType sxx = - 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                RealType syy = + 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                sigma(0,0) = sxx;
                sigma(1,1) = syy;
                return sigma;
            }else{
                RealType sxx = 4*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType sxy = 2*M_PI*std::cos(2*M_PI*y)*std::sin(2*M_PI*x)
                            + 3*M_PI*std::cos(3*M_PI*x)*std::sin(3*M_PI*y);
                RealType syy = 6*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                sigma(0,0) = sxx;
                sigma(0,1) = sxy;
                sigma(1,0) = sxy;
                sigma(1,1) = syy;
                return sigma;
            }
        }

    };

    auto rhs_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType fx = 2*(1 + x*x + y*(-5 + 3*y) + x*(-3 + 4*y));
            RealType fy = 2*(1 + 3*x*x + (-3 + y)*y + x*(-5 + 4*y));
            return static_vector<RealType, 2>{-fx, -fy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType fx = + 2.0 * M_PI * M_PI * ( std::sin(M_PI * x) * std::cos( M_PI * y));
                RealType fy = - 2.0 * M_PI * M_PI * ( std::cos(M_PI * x) * std::sin( M_PI * y));
                return static_vector<RealType, 2>{-fx, -fy};
            }else{
                RealType fx = M_PI*M_PI*(9*(1 + lambda)*std::cos(3*M_PI*x)*std::cos(3*M_PI*y) - 4*(3 + lambda)*std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType fy = M_PI*M_PI*(4*(1 + lambda)*std::cos(2*M_PI*x)*std::cos(2*M_PI*y) - 9*(3 + lambda)*std::sin(3*M_PI*x)*std::sin(3*M_PI*y));
                return static_vector<RealType, 2>{-fx, -fy};
            }
        }
    };

    // simple material
    RealType rho = 1.0;
    RealType vp;
    if (Nonzero_Dirichlet_Q || quadratic_function_Q) {
        vp = std::sqrt(3.0);
    }else{
        vp = std::sqrt(2.0 + lambda);
    }
    RealType vs = 1.0;
    elastic_material_data<RealType> material(rho,vp,vs);
    
    std::ofstream error_file("steady_vector_mixed_two_fields_error.txt");
    
    for(size_t k = 1; k <= sim_data.m_k_degree; k++){
        std::cout << bold << cyan << "Running an approximation with k : " << k << reset << std::endl;
        error_file << "Approximation with k : " << k << std::endl;
        for(size_t l = 0; l <= sim_data.m_n_divs; l++){
            
            // Building a cartesian mesh
            timecounter tc;
            tc.tic();
            RealType lx = 1.0;
            RealType ly = 1.0;
            size_t nx = 2;
            size_t ny = 2;
            mesh_type msh;
            
            cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
            mesh_builder.refine_mesh(l);
            mesh_builder.build_mesh();
            mesh_builder.move_to_mesh_storage(msh);
            std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
            
            // Creating HHO approximation spaces and corresponding linear operator
            size_t cell_k_degree = k;
            if(sim_data.m_hdg_stabilization_Q){
                cell_k_degree++;
            }
            disk::hho_degree_info hho_di(cell_k_degree,k);

            // Solving a scalar primal HHO problem
            boundary_type bnd(msh);
            bnd.addDirichletEverywhere(exact_vec_fun);
            tc.tic();
            auto assembler = elastodynamic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);
            if(sim_data.m_hdg_stabilization_Q){
                assembler.set_hdg_stabilization();
            }
            if(sim_data.m_scaled_stabilization_Q){
                assembler.set_scaled_stabilization();
            }
            assembler.load_material_data(msh,material);
            assembler.assemble(msh, rhs_fun);
            assembler.assemble_mass(msh, false);
            assembler.apply_bc(msh);
            tc.toc();
            std::cout << bold << cyan << "Assemble in : " << tc << " seconds" << reset << std::endl;
            
            // Solving LS
            Matrix<RealType, Dynamic, 1> x_dof;
            if (sim_data.m_sc_Q) {
                tc.tic();
                SparseMatrix<RealType> Kg = assembler.LHS+assembler.MASS;
                linear_solver<RealType> analysis(Kg,assembler.get_n_face_dof());
                analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
                tc.toc();
                
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations (SC) : " << analysis.n_equations() << reset << std::endl;
            }else{
                tc.tic();
                SparseMatrix<RealType> Kg = assembler.LHS+assembler.MASS;
                linear_solver<RealType> analysis(Kg);
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations : " << analysis.n_equations() << reset << std::endl;
            }
            
            // Computing errors
            postprocessor<mesh_type>::compute_errors_two_fields_vectorial(msh, hho_di, x_dof, exact_vec_fun, exact_flux_fun, error_file);
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "steady_vector_mixed_two_fields_k" + std::to_string(k) + "_";
                postprocessor<mesh_type>::write_silo_two_fields_vectorial(silo_file_name, l, msh, hho_di, x_dof, exact_vec_fun, exact_flux_fun, false);
            }
        }
        error_file << std::endl << std::endl;
    }
    error_file.close();
    
    
}

void HHOThreeFieldsConvergenceExample(int argc, char **argv){
    
    using RealType = double;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    
    simulation_data sim_data = preprocessor::process_convergence_test_args(argc, argv);
    sim_data.print_simulation_data();

    // Manufactured exact solution
    bool quadratic_function_Q = false;
    bool Nonzero_Dirichlet_Q = false;
    RealType lambda = 1000.0;
    auto exact_vec_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        RealType x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType ux = (1 - x)*x*(1 - y)*y;
            RealType uy = (1 - x)*x*(1 - y)*y;
            return static_vector<RealType, 2>{ux, uy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType ux = - std::sin(M_PI * pt.x()) * std::cos(M_PI * pt.y());
                RealType uy = + std::cos(M_PI * pt.x()) * std::sin(M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }else{
                RealType ux = std::sin(2.0 * M_PI * pt.x()) * std::sin(2.0 * M_PI * pt.y());
                RealType uy = std::sin(3.0 * M_PI * pt.x()) * std::sin(3.0 * M_PI * pt.y());
                return static_vector<RealType, 2>{ux, uy};
            }

        }
        
    };

    auto exact_flux_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            RealType sxx = 2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            RealType sxy = (1 - x)*x*(1 - y) - (1 - x)*x*y + (1 - x)*(1 - y)*y - x*(1 - y)*y;
            RealType syy = 2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y + (2*(1 - x)*x*(1 - y) - 2*(1 - x)*x*y)/2. + (2*(1 - x)*(1 - y)*y - 2*x*(1 - y)*y)/2.;
            sigma(0,0) = sxx;
            sigma(0,1) = sxy;
            sigma(1,0) = sxy;
            sigma(1,1) = syy;
            return sigma;
        }else{
            static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
            if (Nonzero_Dirichlet_Q) {
                RealType sxx = - 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                RealType syy = + 2.0 * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
                sigma(0,0) = sxx;
                sigma(1,1) = syy;
                return sigma;
            }else{
                RealType sxx = 4*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType sxy = 2*M_PI*std::cos(2*M_PI*y)*std::sin(2*M_PI*x)
                            + 3*M_PI*std::cos(3*M_PI*x)*std::sin(3*M_PI*y);
                RealType syy = 6*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x)
                            + lambda*(3*M_PI*std::cos(3*M_PI*y)*std::sin(3*M_PI*x) + 2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y));
                sigma(0,0) = sxx;
                sigma(0,1) = sxy;
                sigma(1,0) = sxy;
                sigma(1,1) = syy;
                return sigma;
            }
        }

    };

    auto rhs_fun = [quadratic_function_Q,Nonzero_Dirichlet_Q,lambda](const typename mesh_type::point_type& pt) -> static_vector<RealType, 2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            RealType fx = 2*(1 + x*x + y*(-5 + 3*y) + x*(-3 + 4*y));
            RealType fy = 2*(1 + 3*x*x + (-3 + y)*y + x*(-5 + 4*y));
            return static_vector<RealType, 2>{-fx, -fy};
        }else{
            if (Nonzero_Dirichlet_Q) {
                RealType fx = + 2.0 * M_PI * M_PI * ( std::sin(M_PI * x) * std::cos( M_PI * y));
                RealType fy = - 2.0 * M_PI * M_PI * ( std::cos(M_PI * x) * std::sin( M_PI * y));
                return static_vector<RealType, 2>{-fx, -fy};
            }else{
                RealType fx = M_PI*M_PI*(9*(1 + lambda)*std::cos(3*M_PI*x)*std::cos(3*M_PI*y) - 4*(3 + lambda)*std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
                RealType fy = M_PI*M_PI*(4*(1 + lambda)*std::cos(2*M_PI*x)*std::cos(2*M_PI*y) - 9*(3 + lambda)*std::sin(3*M_PI*x)*std::sin(3*M_PI*y));
                return static_vector<RealType, 2>{-fx, -fy};
            }
        }
    };

    // simple material
    RealType rho = 1.0;
    RealType vp;
    if (Nonzero_Dirichlet_Q || quadratic_function_Q) {
        vp = std::sqrt(3.0);
    }else{
        vp = std::sqrt(2.0 + lambda);
    }
    RealType vs = 1.0;
    elastic_material_data<RealType> material(rho,vp,vs);
    
    std::ofstream error_file("steady_vector_mixed_error.txt");
    
    for(size_t k = 1; k <= sim_data.m_k_degree; k++){
        std::cout << bold << cyan << "Running an approximation with k : " << k << reset << std::endl;
        error_file << "Approximation with k : " << k << std::endl;
        for(size_t l = 0; l <= sim_data.m_n_divs; l++){
            
            // Building a cartesian mesh
            timecounter tc;
            tc.tic();
            RealType lx = 1.0;
            RealType ly = 1.0;
            size_t nx = 2;
            size_t ny = 2;
            mesh_type msh;
            
            cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
            mesh_builder.refine_mesh(l);
            mesh_builder.build_mesh();
            mesh_builder.move_to_mesh_storage(msh);
            std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
            
            // Creating HHO approximation spaces and corresponding linear operator
            size_t cell_k_degree = k;
            if(sim_data.m_hdg_stabilization_Q){
                cell_k_degree++;
            }
            disk::hho_degree_info hho_di(cell_k_degree,k);

            // Solving a scalar primal HHO problem
            boundary_type bnd(msh);
            bnd.addDirichletEverywhere(exact_vec_fun);
            tc.tic();
            auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
            if(sim_data.m_hdg_stabilization_Q){
                assembler.set_hdg_stabilization();
            }
            if(sim_data.m_scaled_stabilization_Q){
                assembler.set_scaled_stabilization();
            }
            assembler.load_material_data(msh,material);
            assembler.assemble(msh, rhs_fun);
            assembler.assemble_mass(msh, false);
            assembler.apply_bc(msh);
            tc.toc();
            std::cout << bold << cyan << "Assemble in : " << tc << " seconds" << reset << std::endl;
            
            // Solving LS
            Matrix<RealType, Dynamic, 1> x_dof;
            if (sim_data.m_sc_Q) {
                tc.tic();
                SparseMatrix<RealType> Kg = assembler.LHS+assembler.MASS;
                linear_solver<RealType> analysis(Kg,assembler.get_n_face_dof());
                analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations (SC) : " << analysis.n_equations() << reset << std::endl;
            }else{
                tc.tic();
                SparseMatrix<RealType> Kg = assembler.LHS+assembler.MASS;
                linear_solver<RealType> analysis(Kg);
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                std::cout << bold << cyan << "Number of equations : " << analysis.n_equations() << reset << std::endl;
            }
            
            // Computing errors
            postprocessor<mesh_type>::compute_errors_three_fields_vectorial(msh, hho_di, x_dof, exact_vec_fun, exact_flux_fun, error_file);
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "steady_vector_mixed_k" + std::to_string(k) + "_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, l, msh, hho_di, x_dof, exact_vec_fun, exact_flux_fun, false);
            }
        }
        error_file << std::endl << std::endl;
    }
    error_file.close();
    
    
}

void EHHOFirstOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = tf/nt;
    
    vec_analytic_functions functions;
    functions.set_function_type(vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
    RealType t = ti;
    auto exact_vel_fun      = functions.Evaluate_v(t);
    auto exact_flux_fun     = functions.Evaluate_sigma(t);
    auto rhs_fun            = functions.Evaluate_f(t);
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_vel_fun);
    tc.tic();
    auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
    assembler.load_material_data(msh);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, exact_vel_fun, exact_flux_fun);
    assembler.project_over_faces(msh, x_dof, exact_vel_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "vector_mixed_";
            postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_three_fields_explicit.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    int s = 3;
    Matrix<double, Dynamic, Dynamic> alpha;
    Matrix<double, Dynamic, Dynamic> beta;
    ssprk_shu_osher_tableau::ossprk_tables(s, alpha, beta);

    tc.tic();
    assembler.assemble(msh, rhs_fun);
    tc.toc();
    std::cout << bold << cyan << "Stiffness and rhs assembly completed: " << tc << " seconds" << reset << std::endl;
    size_t n_face_dof = assembler.get_n_face_dof();
    ssprk_hho_scheme<RealType> ssprk_an(assembler.LHS,assembler.RHS,assembler.MASS,n_face_dof);
    tc.toc();

    Matrix<double, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        {   // Updating rhs
            RealType t = dt*(it)+ti;
            auto rhs_fun            = functions.Evaluate_f(t);
            assembler.assemble_rhs(msh, rhs_fun);
            ssprk_an.SetFg(assembler.RHS);
        }
        RealType tn = dt*(it-1)+ti;
        tc.tic();
        {
            
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> ys = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s+1);
        
            Matrix<double, Dynamic, 1> yn, ysi, yj;
            ys.block(0, 0, n_dof, 1) = x_dof;
            for (int i = 0; i < s; i++) {
        
                ysi = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
                for (int j = 0; j <= i; j++) {
                    yn = ys.block(0, j, n_dof, 1);
                    ssprk_an.ssprk_weight(yn, yj, dt, alpha(i,j), beta(i,j));
                    ysi += yj;
                }
                ys.block(0, i+1, n_dof, 1) = ysi;
            }
        
            x_dof_n = ys.block(0, s, n_dof, 1);
        }
        tc.toc();
        std::cout << bold << cyan << "SSPRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        t = tn + dt;
        auto exact_vel_fun = functions.Evaluate_v(t);
        auto exact_flux_fun = functions.Evaluate_sigma(t);

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "vector_mixed_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
        }
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

        if(it == nt){
            // Computing errors
            postprocessor<mesh_type>::compute_errors_three_fields_vectorial(msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun);
        }
    }
    
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of SSPRKSS steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
}

void IHHOFirstOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;
    
    vec_analytic_functions functions;
    functions.set_function_type(vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
    RealType t = ti;
    auto exact_vel_fun      = functions.Evaluate_v(t);
    auto exact_flux_fun     = functions.Evaluate_sigma(t);
    auto rhs_fun            = functions.Evaluate_f(t);
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_vel_fun);
    tc.tic();
    auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
    assembler.load_material_data(msh);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, exact_vel_fun, exact_flux_fun);
    assembler.project_over_faces(msh, x_dof, exact_vel_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "vector_mixed_";
            postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_three_fields.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 1;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, rhs_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
                auto exact_vel_fun      = functions.Evaluate_v(t);
                auto rhs_fun            = functions.Evaluate_f(t);
                assembler.get_bc_conditions().updateDirichletFunction(exact_vel_fun, 0);
                assembler.assemble_rhs(msh, rhs_fun);
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "DIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        t = tn + dt;
        auto exact_vel_fun = functions.Evaluate_v(t);
        auto exact_flux_fun = functions.Evaluate_sigma(t);

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "vector_mixed_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
        }
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

        if(it == nt){
            // Computing errors
            postprocessor<mesh_type>::compute_errors_three_fields_vectorial(msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
    simulation_log << "Number of DIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

void IHHOFirstOrderTwoFields(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;
    
    vec_analytic_functions functions;
    functions.set_function_type(vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
    RealType t = ti;
    auto exact_vel_fun      = functions.Evaluate_v(t);
    auto exact_flux_fun     = functions.Evaluate_sigma(t);
    auto rhs_fun            = functions.Evaluate_f(t);
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_vel_fun);
    tc.tic();
    auto assembler = elastodynamic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);
    assembler.load_material_data(msh);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, exact_vel_fun, exact_flux_fun);
    assembler.project_over_faces(msh, x_dof, exact_vel_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "vector_mixed_two_fields";
            postprocessor<mesh_type>::write_silo_two_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_two_fields.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, rhs_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
                auto exact_vel_fun      = functions.Evaluate_v(t);
                auto rhs_fun            = functions.Evaluate_f(t);
                assembler.get_bc_conditions().updateDirichletFunction(exact_vel_fun, 0);
                assembler.assemble_rhs(msh, rhs_fun);
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "DIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        t = tn + dt;
        auto exact_vel_fun = functions.Evaluate_v(t);
        auto exact_flux_fun = functions.Evaluate_sigma(t);

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "vector_mixed__two_fields";
                postprocessor<mesh_type>::write_silo_two_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
        }
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

        if(it == nt){
            // Computing errors
            postprocessor<mesh_type>::compute_errors_three_fields_vectorial(msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
    simulation_log << "Number of DIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

void IHHOSecondOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);

    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;

    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.125;
    RealType dt     = (tf-ti)/nt;

    vec_analytic_functions functions;
    functions.set_function_type(vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
    RealType t = ti;
    auto exact_vec_fun      = functions.Evaluate_u(t);
    auto exact_vel_fun      = functions.Evaluate_v(t);
    auto exact_accel_fun    = functions.Evaluate_a(t);
    auto rhs_fun            = functions.Evaluate_f(t);
    auto exact_flux_fun     = functions.Evaluate_sigma(t);

    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_vec_fun);

    tc.tic();
    auto assembler = elastodynamic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    assembler.load_material_data(msh);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler created: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;

    // Projecting initial displacement, velocity and acceleration
    tc.tic();
    Matrix<RealType, Dynamic, 1> u_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, u_dof_n, exact_vec_fun);
    assembler.project_over_cells(msh, v_dof_n, exact_vel_fun);
    assembler.project_over_cells(msh, a_dof_n, exact_accel_fun);
    
    assembler.project_over_faces(msh, u_dof_n, exact_vec_fun);
    assembler.project_over_faces(msh, v_dof_n, exact_vel_fun);
    assembler.project_over_faces(msh, a_dof_n, exact_accel_fun);
    tc.toc();
    std::cout << bold << cyan << "Initialization completed: " << tc << " seconds" << reset << std::endl;
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "vec_";
        postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, u_dof_n, exact_vec_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_one_field.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, t, u_dof_n, v_dof_n, simulation_log);
    }
    
    bool standar_Q = true;
    // Newmark process
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.275;
        RealType gamma = 0.55;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 1.5;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, rhs_fun);
        SparseMatrix<RealType> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        linear_solver<RealType> analysis;
        if (sim_data.m_sc_Q) {
            analysis.set_Kg(assembler.LHS, assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
        }else{
            analysis.set_Kg(assembler.LHS);
        }
        analysis.set_direct_solver(true);
//        analysis.set_iterative_solver(true);
        analysis.factorize();
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            // Manufactured solution
            RealType t = dt*(it)+ti;
            auto exact_vec_fun      = functions.Evaluate_u(t);
            auto exact_vel_fun      = functions.Evaluate_v(t);
            auto exact_accel_fun    = functions.Evaluate_a(t);
            auto exact_flux_fun     = functions.Evaluate_sigma(t);
            auto rhs_fun            = functions.Evaluate_f(t);
            assembler.get_bc_conditions().updateDirichletFunction(exact_vec_fun, 0);
            assembler.assemble_rhs(msh, rhs_fun);
            
            // Compute intermediate state for scalar and rate
            u_dof_n = u_dof_n + dt*v_dof_n + 0.5*dt*dt*(1.0-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1.0-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*u_dof_n;

            assembler.RHS -= res;
            tc.toc();
            std::cout << bold << cyan << "Rhs assembly completed: " << tc << " seconds" << reset << std::endl;

            tc.tic();
            a_dof_np = analysis.solve(assembler.RHS); // new acceleration
            tc.toc();
            std::cout << bold << cyan << "Solution completed: " << tc << " seconds" << reset << std::endl;

            // update displacement, velocity and acceleration
            u_dof_n += beta*dt*dt*a_dof_np;
            v_dof_n += gamma*dt*a_dof_np;
            a_dof_n  = a_dof_np;
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "vec_";
                postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, u_dof_n, exact_vec_fun, false);
            }
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, t, u_dof_n, v_dof_n, simulation_log);
            }

            if(it == nt){
                postprocessor<mesh_type>::compute_errors_one_field_vectorial(msh, hho_di, assembler, u_dof_n, exact_vec_fun, exact_flux_fun, simulation_log);
            }

        }
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
    
}

void HeterogeneousIHHOSecondOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 3.0;
    RealType ly = 2.5;
    size_t nx = 3;
    size_t ny = 3;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-1.0, 0.0);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;

    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.75;
    RealType dt     = (tf-ti)/nt;
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = (2.0/3.0)+1.25;
            c = 10.0;
            lp = std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);

    tc.tic();
    auto assembler = elastodynamic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
         double x,y;
         x = pt.x();
         y = pt.y();
         std::vector<RealType> mat_data(3);
         RealType rho, vp, vs;
         rho = 1.0;
         if (y < 1.25) {
             vp = 1.0*std::sqrt(3.0);
             vs  = 1.0;
         }else{
             vp = std::sqrt(3.0);
             vs  = 1;
         }
         mat_data[0] = rho; // rho
         mat_data[1] = vp; // seismic compressional velocity vp
         mat_data[2] = vs; // seismic shear velocity vp
         return mat_data;
     };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler created: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;

    // Projecting initial displacement, velocity and acceleration
    tc.tic();
    Matrix<RealType, Dynamic, 1> u_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, u_dof_n, null_fun);
    assembler.project_over_faces(msh, u_dof_n, null_fun);
    
    assembler.project_over_cells(msh, v_dof_n, vec_fun);
    assembler.project_over_faces(msh, v_dof_n, vec_fun);
    
    assembler.project_over_cells(msh, a_dof_n, null_fun);
    assembler.project_over_faces(msh, a_dof_n, null_fun);
    tc.toc();
    std::cout << bold << cyan << "Initialization completed: " << tc << " seconds" << reset << std::endl;
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vec_";
        postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_one_field.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_one_field.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_one_field.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_one_field.csv");
    typename mesh_type::point_type s1_pt(0.5-2.0/3.0, 1.0/3.0);
    typename mesh_type::point_type s2_pt(0.5, 1.0/3.0);
    typename mesh_type::point_type s3_pt(0.5+2.0/3.0, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, ti, u_dof_n, v_dof_n, simulation_log);
    }
    
    linear_solver<RealType> analysis;
    bool standar_Q = true;
    // Newmark process
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.25;
        RealType gamma = 0.5;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 1.5;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, null_fun);
        SparseMatrix<RealType> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        if (sim_data.m_sc_Q) {
            tc.tic();
            analysis.set_Kg(assembler.LHS,assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
            tc.toc();
            std::cout << bold << cyan << "Equations condensed in : " << tc << " seconds" << reset << std::endl;
            
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
        
        }else{
            analysis.set_Kg(assembler.LHS);
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
            
        }
        
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            RealType t = dt*it+ti;
            tc.tic();
            // Compute intermediate state for scalar and rate
            u_dof_n = u_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*u_dof_n;

            assembler.RHS.setZero();
            assembler.RHS -= res;
            tc.toc();
            std::cout << bold << cyan << "Rhs assembly completed: " << tc << " seconds" << reset << std::endl;

            tc.tic();
            a_dof_np = analysis.solve(assembler.RHS); // new acceleration
            tc.toc();
            std::cout << bold << cyan << "Solution completed: " << tc << " seconds" << reset << std::endl;

            // update displacement, velocity and acceleration
            u_dof_n += beta*dt*dt*a_dof_np;
            v_dof_n += gamma*dt*a_dof_np;
            a_dof_n  = a_dof_np;

            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "inhomogeneous_vec_";
                postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
            }
            
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, t, u_dof_n, v_dof_n, simulation_log);
            }

        }
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
    
}

void HeterogeneousIHHOFirstOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;
    
    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
        return sigma;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 0.75;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (-4*std::sqrt(10.0/3.0)*(-1 + 1600.0*r*r))/(std::exp(800*r*r)*std::pow(M_PI,0.25));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(3);
        RealType rho, vp, vs;
        rho = 1.0;
        if (y < 0.5) {
            vp = 1.0*std::sqrt(3.0);
            vs  = 1.0;
        }else{
            vp = std::sqrt(3.0);
            vs  = 1;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
        mat_data[2] = vs; // seismic shear velocity vp
        return mat_data;
    };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, vec_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vec_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vector_mixed_";
            postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_three_fields.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, null_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    
    Matrix<double, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
                assembler.assemble_rhs(msh, null_fun);
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "DIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        RealType t = tn + dt;

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "inhomogeneous_vector_mixed_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
        }
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of DIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

void Gar6more2DIHHOSecondOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 2.0;
    RealType ly = 1.5;
    size_t nx = 3;
    size_t ny = 3;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-0.5, 0.0);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;

    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.75;
    RealType dt     = (tf-ti)/nt;
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 2.0/3.0;
            c = 10.0;
            lp = 1.0*std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);

    tc.tic();
    auto assembler = elastodynamic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
         double x,y;
         x = pt.x();
         y = pt.y();
         std::vector<RealType> mat_data(3);
         RealType rho, vp, vs;
         rho = 1.0;
         if (y < 0.5) {
             vp = 1.0*std::sqrt(3.0);
             vs  = 1.0;
         }else{
             vp = 1.0*std::sqrt(3.0);
             vs  = 1.0;
         }
         mat_data[0] = rho; // rho
         mat_data[1] = vp; // seismic compressional velocity vp
         mat_data[2] = vs; // seismic shear velocity vp
         return mat_data;
     };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler created: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;

    // Projecting initial displacement, velocity and acceleration
    tc.tic();
    Matrix<RealType, Dynamic, 1> u_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, u_dof_n, null_fun);
    assembler.project_over_faces(msh, u_dof_n, null_fun);
    
    assembler.project_over_cells(msh, v_dof_n, vec_fun);
    assembler.project_over_faces(msh, v_dof_n, vec_fun);
    
    assembler.project_over_cells(msh, a_dof_n, null_fun);
    assembler.project_over_faces(msh, a_dof_n, null_fun);
    tc.toc();
    std::cout << bold << cyan << "Initialization completed: " << tc << " seconds" << reset << std::endl;
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vec_";
        postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_one_field.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_one_field.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_one_field.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_one_field.csv");
    typename mesh_type::point_type s1_pt(0.5-2.0/3.0, 1.0/3.0);
    typename mesh_type::point_type s2_pt(0.5, 1.0/3.0);
    typename mesh_type::point_type s3_pt(0.5+2.0/3.0, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, ti, u_dof_n, v_dof_n, simulation_log);
    }
    
    linear_solver<RealType> analysis;
    bool standar_Q = true;
    // Newmark process
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.25;
        RealType gamma = 0.5;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 1.5;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, null_fun);
        SparseMatrix<RealType> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        if (sim_data.m_sc_Q) {
            tc.tic();
            analysis.set_Kg(assembler.LHS,assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
            tc.toc();
            std::cout << bold << cyan << "Equations condensed in : " << tc << " seconds" << reset << std::endl;
            
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
        
        }else{
            analysis.set_Kg(assembler.LHS);
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
            
        }
        
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            RealType t = dt*it+ti;
            tc.tic();
            // Compute intermediate state for scalar and rate
            u_dof_n = u_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*u_dof_n;

            assembler.RHS.setZero();
            assembler.RHS -= res;
            tc.toc();
            std::cout << bold << cyan << "Rhs assembly completed: " << tc << " seconds" << reset << std::endl;

            tc.tic();
            a_dof_np = analysis.solve(assembler.RHS); // new acceleration
            tc.toc();
            std::cout << bold << cyan << "Solution completed: " << tc << " seconds" << reset << std::endl;

            // update displacement, velocity and acceleration
            u_dof_n += beta*dt*dt*a_dof_np;
            v_dof_n += gamma*dt*a_dof_np;
            a_dof_n  = a_dof_np;

            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "inhomogeneous_vec_";
                postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
            }
            
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, t, u_dof_n, v_dof_n, simulation_log);
            }

        }
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
    
}

void Gar6more2DIHHOFirstOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 2.0;
    RealType ly = 1.5;
    size_t nx = 3;
    size_t ny = 3;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-0.5, 0.0);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.75;
    RealType dt     = (tf-ti)/nt;
    
    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
        return sigma;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 2.0/3.0;
            c = 10.0;
            lp = 1.0*std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(3);
        RealType rho, vp, vs;
        rho = 1.0;
        if (y < 0.5) {
            vp = 1.0*std::sqrt(3.0);
            vs  = 1.0;
        }else{
            vp = 1.0*std::sqrt(3.0);
            vs  = 1.0;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
        mat_data[2] = vs; // seismic shear velocity vp
        return mat_data;
    };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, vec_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vec_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vector_mixed_";
            postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_three_fields.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_three_fields.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_three_fields.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_three_fields.csv");
    typename mesh_type::point_type s1_pt(0.5-2.0/3.0, 1.0/3.0);
    typename mesh_type::point_type s2_pt(0.5, 1.0/3.0);
    typename mesh_type::point_type s3_pt(0.5+2.0/3.0, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, null_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    
    Matrix<double, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
//                assembler.assemble_rhs(msh, null_fun);
                assembler.RHS.setZero();
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "SDIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        RealType t = tn + dt;

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "inhomogeneous_vector_mixed_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
        }
        
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of SDIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

void HeterogeneousGar6more2DIHHOSecondOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();

    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 3.0;
    RealType ly = 3.0;
    size_t nx = 3;
    size_t ny = 3;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-1.5, -1.5);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;

    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.0;
            yc = 2.0/3.0;
            c = 10.0;
            lp = 2.0*std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);

    tc.tic();
    auto assembler = elastodynamic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
         double x,y;
         x = pt.x();
         y = pt.y();
         std::vector<RealType> mat_data(3);
         RealType rho, vp, vs;
         if (y > 0.0) {
             vp = 2.0*std::sqrt(3.0);
             vs  = 2.0;
             rho = 1.0;
         }else{
             vp = std::sqrt(3.0);
             vs  = 1.0;
             rho = 1.0;
         }
         mat_data[0] = rho; // rho
         mat_data[1] = vp; // seismic compressional velocity vp
         mat_data[2] = vs; // seismic shear velocity vp
         return mat_data;
     };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler created: " << tc << " seconds" << reset << std::endl;

    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;

    // Projecting initial displacement, velocity and acceleration
    tc.tic();
    Matrix<RealType, Dynamic, 1> u_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, u_dof_n, null_fun);
    assembler.project_over_faces(msh, u_dof_n, null_fun);
    
    assembler.project_over_cells(msh, v_dof_n, vec_fun);
    assembler.project_over_faces(msh, v_dof_n, vec_fun);
    
    assembler.project_over_cells(msh, a_dof_n, null_fun);
    assembler.project_over_faces(msh, a_dof_n, null_fun);
    tc.toc();
    std::cout << bold << cyan << "Initialization completed: " << tc << " seconds" << reset << std::endl;
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vec_";
        postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_one_field.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_one_field_h.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_one_field_h.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_one_field_h.csv");
    typename mesh_type::point_type s1_pt(-1.0/3.0, -1.0/3.0);
    typename mesh_type::point_type s2_pt( 0.0, -1.0/3.0);
    typename mesh_type::point_type s3_pt(+1.0/3.0, -1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_one_field(0, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, ti, u_dof_n, v_dof_n, simulation_log);
    }
    
    linear_solver<RealType> analysis;
    bool standar_Q = true;
    // Newmark process
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.25;
        RealType gamma = 0.5;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 1.5;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, null_fun);
        SparseMatrix<RealType> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        if (sim_data.m_sc_Q) {
            tc.tic();
            analysis.set_Kg(assembler.LHS,assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
            tc.toc();
            std::cout << bold << cyan << "Equations condensed in : " << tc << " seconds" << reset << std::endl;
            
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
        
        }else{
            analysis.set_Kg(assembler.LHS);
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
            
        }
        
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            RealType t = dt*it+ti;
            tc.tic();
            // Compute intermediate state for scalar and rate
            u_dof_n = u_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*u_dof_n;

            assembler.RHS.setZero();
            assembler.RHS -= res;
            tc.toc();
            std::cout << bold << cyan << "Rhs assembly completed: " << tc << " seconds" << reset << std::endl;

            tc.tic();
            a_dof_np = analysis.solve(assembler.RHS); // new acceleration
            tc.toc();
            std::cout << bold << cyan << "Solution completed: " << tc << " seconds" << reset << std::endl;

            // update displacement, velocity and acceleration
            u_dof_n += beta*dt*dt*a_dof_np;
            v_dof_n += gamma*dt*a_dof_np;
            a_dof_n  = a_dof_np;

            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "inhomogeneous_vec_";
                postprocessor<mesh_type>::write_silo_one_field_vectorial(silo_file_name, it, msh, hho_di, v_dof_n, vec_fun, false);
            }
            
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s1_pt_cell, msh, hho_di, v_dof_n, sensor_1_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s2_pt_cell, msh, hho_di, v_dof_n, sensor_2_log);
            postprocessor<mesh_type>::record_data_elastic_one_field(it, s3_pt_cell, msh, hho_di, v_dof_n, sensor_3_log);
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_elastic_energy_one_field(msh, hho_di, assembler, t, u_dof_n, v_dof_n, simulation_log);
            }

        }
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
    
}

void HeterogeneousGar6more2DIHHOFirstOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 3.0;
    RealType ly = 3.0;
    size_t nx = 3;
    size_t ny = 3;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-1.5, -1.5);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;
    
    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
        return sigma;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.0;
            yc = 2.0/3.0;
            c = 10.0;
            lp = 2.0*std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = elastodynamic_three_fields_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(3);
        RealType rho, vp, vs;
        if (y > 0.0) {
            vp = 2.0*std::sqrt(3.0);
            vs  = 2.0;
            rho = 1.0;
        }else{
            vp = std::sqrt(3.0);
            vs  = 1.0;
            rho = 1.0;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
        mat_data[2] = vs; // seismic shear velocity vp
        return mat_data;
    };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, vec_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vec_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vector_mixed_";
            postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_three_fields.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_three_fields_h.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_three_fields_h.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_three_fields_h.csv");
    typename mesh_type::point_type s1_pt(-1.0/3.0, 1.0/3.0);
    typename mesh_type::point_type s2_pt( 0.0, 1.0/3.0);
    typename mesh_type::point_type s3_pt(+1.0/3.0, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_three_fields(0, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, null_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    
    Matrix<double, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
                assembler.RHS.setZero();
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "SDIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        RealType t = tn + dt;

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "inhomogeneous_vector_mixed_";
                postprocessor<mesh_type>::write_silo_three_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
        }
        
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
        postprocessor<mesh_type>::record_data_elastic_three_fields(it, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of SDIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

void HeterogeneousGar6more2DIHHOFirstOrderTwoFields(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 2.0;
    RealType ly = 2.0;
    size_t nx = 2;
    size_t ny = 2;
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> boundary_type;
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-1.0, -1.25);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    tc.toc();
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.5;
    RealType dt     = (tf-ti)/nt;
    
    auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            static_vector<RealType, 2> f{0,0};
            return f;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
        return sigma;
    };
    
    auto vec_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
            RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
            x = pt.x();
            y = pt.y();
            xc = 0.0;
            yc = 0.25;
            c = 10.0;
            lp = 2.0*std::sqrt(3.0)/10.0;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
            vx = wave*(x-xc);
            vy = wave*(y-yc);
            static_vector<RealType, 2> v{vx,vy};
            return v;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = elastodynamic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(3);
        RealType rho, vp, vs;
        if (y > 0.0) {
            vp = std::sqrt(3.0);
            vs  = 1.0;
            rho = 1.0;
        }else{
            vp = std::sqrt(3.0);
            vs  = 1.0;
            rho = 1.0;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
        mat_data[2] = vs; // seismic shear velocity vp
        return mat_data;
    };
    
    assembler.load_material_data(msh,elastic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, vec_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vec_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_vector_mixed_two_fields_";
            postprocessor<mesh_type>::write_silo_two_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("elastodynamic_inhomogeneous_two_fields.txt");
    
    std::ofstream sensor_1_log("s1_elastodynamic_two_fields_h.csv");
    std::ofstream sensor_2_log("s2_elastodynamic_two_fields_h.csv");
    std::ofstream sensor_3_log("s3_elastodynamic_two_fields_h.csv");
    typename mesh_type::point_type s1_pt(-1.0/3.0, 1.0/3.0);
    typename mesh_type::point_type s2_pt( 0.0, 1.0/3.0);
    typename mesh_type::point_type s3_pt(+1.0/3.0, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_data_elastic_two_fields(0, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
    postprocessor<mesh_type>::record_data_elastic_two_fields(0, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
    postprocessor<mesh_type>::record_data_elastic_two_fields(0, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_elastic_energy_two_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;

    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;

    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }

    tc.tic();
    assembler.assemble(msh, null_fun);
    tc.toc();
    std::cout << bold << cyan << "First stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
//        dirk_an.setIterativeSolver();
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }
    
    Matrix<double, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;

        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);

            RealType t;
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {

                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }

                t = tn + c(i,0) * dt;
                assembler.RHS.setZero();
                dirk_an.SetFg(assembler.RHS);
                dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "SDIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        RealType t = tn + dt;

        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "inhomogeneous_vector_mixed_two_fields_";
                postprocessor<mesh_type>::write_silo_two_fields_vectorial(silo_file_name, it, msh, hho_di, x_dof, vec_fun, null_flux_fun, false);
        }
        
        postprocessor<mesh_type>::record_data_elastic_two_fields(it, s1_pt_cell, msh, hho_di, x_dof, sensor_1_log);
        postprocessor<mesh_type>::record_data_elastic_two_fields(it, s2_pt_cell, msh, hho_di, x_dof, sensor_2_log);
        postprocessor<mesh_type>::record_data_elastic_two_fields(it, s3_pt_cell, msh, hho_di, x_dof, sensor_3_log);
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_elastic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of SDIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

