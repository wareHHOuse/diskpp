  //
//  postprocessor.hpp
//  acoustics
//
//  Created by Omar Durán on 4/10/20.
//  Contributor: Romain Mottier

#pragma once
#ifndef postprocessor_hpp
#define postprocessor_hpp

#include <iomanip>
#include "../common/acoustic_one_field_assembler.hpp"
#include "../common/acoustic_two_fields_assembler.hpp"
#include "../common/elastodynamic_one_field_assembler.hpp"
#include "../common/elastodynamic_two_fields_assembler.hpp"
#include "../common/elastodynamic_three_fields_assembler.hpp"
#include "../common/elastoacoustic_two_fields_assembler.hpp"
#include "../common/elastoacoustic_four_fields_assembler.hpp"
#include "diskpp/bases/bases.hpp"

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif


template<typename Mesh>
class postprocessor {
    
public:
    
    // FIND CELLS & PICK CELLS & MESH QUALITY
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////

    static std::set<size_t> find_cells(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false) {
        
        using RealType = double;
        auto norm =  [](const typename Mesh::point_type& a, const typename Mesh::point_type& b ) -> RealType {
            RealType dx = (b.x() - a.x());
            RealType dy = (b.y() - a.y());
            RealType norm = std::sqrt(dx*dx + dy*dy);
            return norm;
        };
        
        // find minimum distance to the requested point
        size_t np = msh.points_size();
        std::vector<RealType> distances(np);
        
        size_t ip = 0;
        for (auto& point : msh.backend_storage()->points) {
            RealType dist = norm(pt,point);
            distances[ip] = dist;
            ip++;
        }
        
        size_t index = std::min_element(distances.begin(),distances.end()) - distances.begin();
        if(verbose_Q){
            RealType min_dist = *std::min_element(distances.begin(), distances.end());
            typename Mesh::point_type nearest_point = msh.backend_storage()->points.at(index);
        }
        
        std::set<size_t> cell_indexes;
        size_t cell_i = 0;
        for (auto& cell : msh) {
            auto points = cell.point_ids();
            size_t n_p = points.size();
            for (size_t l = 0; l < n_p; l++) {
                auto pt_id = points[l];
                if (index == pt_id)
                    cell_indexes.insert(cell_i);
            }
            cell_i++;
        }
        
        if(verbose_Q){
            // std::cout << "Detected cells indexes : " << std::endl;
            for(auto index : cell_indexes){
                // std::cout << index << std::endl;
            }
        }
        
        return cell_indexes;
    }
    
    static std::set<size_t> find_elastic_cells(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false) {
        
        using RealType = double;
        auto norm = [](const typename Mesh::point_type& a, const typename Mesh::point_type& b) -> RealType {
            RealType dx = (b.x() - a.x());
            RealType dy = (b.y() - a.y());
            return std::sqrt(dx * dx + dy * dy);
        };
        
        // Trouver le point du maillage le plus proche
        size_t np = msh.points_size();
        std::vector<RealType> distances(np);
        
        size_t ip = 0;
        for (auto& point : msh.backend_storage()->points) {
            distances[ip] = norm(pt, point);
            ip++;
        }
        
        size_t index = std::min_element(distances.begin(), distances.end()) - distances.begin();
    
        if (verbose_Q) {
            RealType min_dist = *std::min_element(distances.begin(), distances.end());
            typename Mesh::point_type nearest_point = msh.backend_storage()->points.at(index);
        }
        
        std::set<size_t> cell_indexes;
        size_t cell_i = 0;
        for (auto& cell : msh) {
            auto points = cell.point_ids();
            size_t n_p = points.size();
            
            for (size_t l = 0; l < n_p; l++) {
                if (index == points[l]) {
                    auto bar = barycenter(msh, cell);
                    if (bar.y() <= 0.0) 
                        cell_indexes.insert(cell_i);
                }
            }
            cell_i++;
        }   
        return cell_indexes;
    }

    static std::set<size_t> find_acoustic_cells(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false) {
        
        using RealType = double;
        auto norm = [](const typename Mesh::point_type& a, const typename Mesh::point_type& b) -> RealType {
            RealType dx = (b.x() - a.x());
            RealType dy = (b.y() - a.y());
            return std::sqrt(dx * dx + dy * dy);
        };
        
        // Trouver le point du maillage le plus proche
        size_t np = msh.points_size();
        std::vector<RealType> distances(np);
        
        size_t ip = 0;
        for (auto& point : msh.backend_storage()->points) {
            distances[ip] = norm(pt, point);
            ip++;
        }
        
        size_t index = std::min_element(distances.begin(), distances.end()) - distances.begin();
        
        if (verbose_Q) {
            RealType min_dist = *std::min_element(distances.begin(), distances.end());
            typename Mesh::point_type nearest_point = msh.backend_storage()->points.at(index);
        }
        
        std::set<size_t> cell_indexes;
        size_t cell_i = 0;
        for (auto& cell : msh) {
            auto points = cell.point_ids();
            size_t n_p = points.size();
            
            for (size_t l = 0; l < n_p; l++) {
                if (index == points[l]) {
                    auto bar = barycenter(msh, cell);
                    if (bar.y() >= 0) 
                        cell_indexes.insert(cell_i);
                }
            }
            cell_i++;
        }
        
        return cell_indexes;
    }
    
    static void mesh_quality(Mesh & msh, elastoacoustic_four_fields_assembler<Mesh> & assembler,std::ostream & mesh_file = std::cout){
    
        using RealType = double;
        auto dim = Mesh::dimension;
        
        RealType o_to_ten          = 0.0;
        RealType ten_to_twenty     = 0.0;
        RealType twenty_to_thirty  = 0.0;
        RealType thirty_to_forty   = 0.0;
        RealType forty_to_fifty    = 0.0;
        RealType fifty_to_sixty    = 0.0;
        RealType sixty_to_seventy  = 0.0;
        RealType seventy_to_eighty = 0.0;
        RealType eighty_to_ninety  = 0.0;
        RealType ninety_to_hundred = 0.0;
        
        // Determination h 
        RealType h = 10.0;
        auto storage = msh.backend_storage();
        for (auto& e_chunk : assembler.get_e_material_data()) {
            auto& cell = storage->surfaces[e_chunk.first];
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }
        }
        for (auto& a_chunk : assembler.get_a_material_data()) {
            auto& cell = storage->surfaces[a_chunk.first];
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }
        }
        
        // Determination of small edges 
        RealType nb_faces = 0.0;
        for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
            const auto face = *face_it;
            auto pts = points(msh, face);
            RealType longueur = (pts[1] - pts[0]).to_vector().norm();
            if (longueur <= 1.e-1*h) {
                o_to_ten++;
            } else if (longueur <= 2e-1*h) {
                ten_to_twenty++;
            } else if (longueur <= 3e-1*h) {
                twenty_to_thirty++;
            } else if (longueur <= 4e-1*h) {
                thirty_to_forty++;
            } else if (longueur <= 5e-1*h) {
                forty_to_fifty++;
            } else if (longueur <= 6e-1*h) {
                fifty_to_sixty++;
            } else if (longueur <= 7e-1*h) {
                sixty_to_seventy++;
            } else if (longueur <= 8e-1*h) {
                seventy_to_eighty++;
            } else if (longueur <= 9e-1*h) {
                eighty_to_ninety++;
            } else {
                ninety_to_hundred++;
            } 
            nb_faces++;
        }
        // Number of faces 
        mesh_file << "Number of faces: " << nb_faces << std::endl<< std::endl;
        
        mesh_file << "      ratio <= 10% : " << "   " << o_to_ten << std::endl;
        mesh_file << "10% < ratio <= 20% : " << "   " << ten_to_twenty << std::endl;
        mesh_file << "20% < ratio <= 30% : " << "   " << twenty_to_thirty << std::endl;
        mesh_file << "30% < ratio <= 40% : " << "   " << thirty_to_forty << std::endl;
        mesh_file << "40% < ratio <= 50% : " << "   " << forty_to_fifty << std::endl;
        mesh_file << "50% < ratio <= 60% : " << "   " << fifty_to_sixty << std::endl;
        mesh_file << "60% < ratio <= 70% : " << "   " << sixty_to_seventy << std::endl;
        mesh_file << "70% < ratio <= 80% : " << "   " << seventy_to_eighty << std::endl;
        mesh_file << "80% < ratio <= 90% : " << "   " << eighty_to_ninety << std::endl;
        mesh_file << "90% < ratio <= 100%: " << "   " << ninety_to_hundred << std::endl;
        
        RealType prct_o_to_ten          = o_to_ten/nb_faces;
        RealType prct_ten_to_twenty     = ten_to_twenty/nb_faces;
        RealType prct_twenty_to_thirty  = twenty_to_thirty/nb_faces;
        RealType prct_thirty_to_forty   = thirty_to_forty/nb_faces;
        RealType prct_forty_to_fifty    = forty_to_fifty/nb_faces;
        RealType prct_fifty_to_sixty    = fifty_to_sixty/nb_faces;
        RealType prct_sixty_to_seventy  = sixty_to_seventy/nb_faces;
        RealType prct_seventy_to_eighty = seventy_to_eighty/nb_faces;
        RealType prct_eighty_to_ninety  = eighty_to_ninety/nb_faces;
        RealType prct_ninety_to_hundred = ninety_to_hundred/nb_faces;
        
        // Percentage of faces 
        mesh_file << "Percentage: " << std::endl;
        mesh_file << "      ratio <= 10% : " << "   " << prct_o_to_ten << std::endl;
        mesh_file << "10% < ratio <= 20% : " << "   " << prct_ten_to_twenty << std::endl;
        mesh_file << "20% < ratio <= 30% : " << "   " << prct_twenty_to_thirty << std::endl;
        mesh_file << "30% < ratio <= 40% : " << "   " << prct_thirty_to_forty << std::endl;
        mesh_file << "40% < ratio <= 50% : " << "   " << prct_forty_to_fifty << std::endl;
        mesh_file << "50% < ratio <= 60% : " << "   " << prct_fifty_to_sixty << std::endl;
        mesh_file << "60% < ratio <= 70% : " << "   " << prct_sixty_to_seventy << std::endl;
        mesh_file << "70% < ratio <= 80% : " << "   " << prct_seventy_to_eighty << std::endl;
        mesh_file << "80% < ratio <= 90% : " << "   " << prct_eighty_to_ninety << std::endl;
        mesh_file << "90% < ratio <= 100%: " << "   " << prct_ninety_to_hundred << std::endl;
        RealType total = prct_o_to_ten + prct_ten_to_twenty + prct_twenty_to_thirty + prct_thirty_to_forty + prct_forty_to_fifty + prct_fifty_to_sixty + prct_sixty_to_seventy + prct_seventy_to_eighty + prct_eighty_to_ninety + prct_ninety_to_hundred;
        mesh_file << "              TOTAL: " << "   " << total << std::endl << std::endl;
        
    }

    static size_t pick_cell(typename Mesh::point_type & pt, Mesh & msh, std::set<size_t> & cell_indexes, bool verbose_Q = false){ 
        
        using RealType = double;
        
        auto triangle_member_Q = [] (typename Mesh::point_type & p, typename Mesh::point_type & p0, typename Mesh::point_type & p1, typename Mesh::point_type & p2) {
            RealType dx = p.x()-p2.x();
            RealType dy = p.y()-p2.y();
            RealType dx21 = p2.x()-p1.x();
            RealType dy12 = p1.y()-p2.y();
            RealType d = dy12*(p0.x()-p2.x()) + dx21*(p0.y()-p2.y());
            RealType s = dy12*dx + dx21*dy;
            RealType t = (p2.y()-p0.y())*dx + (p0.x()-p2.x())*dy;
            if (d < 0.0) {
                return s<=0.0 && t<=0.0 && s+t>=d;
            }
            return s>=0 && t>=0 && s+t<=d;
        };
        
        size_t n_cells = cell_indexes.size();
        if (n_cells == 1) {
            size_t first_index = *cell_indexes.begin();
            return first_index;
        }
        bool is_member_Q = false;
        for(auto index : cell_indexes){
            auto& cell = msh.backend_storage()->surfaces[index];
            auto bar = barycenter(msh, cell);
            auto points = cell.point_ids();
            size_t n_p = points.size();
            
            // building teselation
            std::vector<std::vector<typename Mesh::point_type>> triangles(n_p);
            for (size_t l = 0; l < n_p; l++) {
                std::vector<typename Mesh::point_type> chunk(3);
                if( l == n_p - 1){
                    chunk[0] = msh.backend_storage()->points.at(points[l]);
                    chunk[1] = msh.backend_storage()->points.at(points[0]);
                    chunk[2] = bar;
                }
                else {
                    chunk[0] = msh.backend_storage()->points.at(points[l]);
                    chunk[1] = msh.backend_storage()->points.at(points[l+1]);
                    chunk[2] = bar;
                }
                triangles[l] = chunk;
            }
            // check whether the point is memeber of any triangle
            for (auto triangle : triangles) {
                is_member_Q = triangle_member_Q(pt,triangle[0],triangle[1],triangle[2]);
                if (is_member_Q) {
                    // std::cout << "Detected cell index = " << index << std::endl;
                    return index;
                }
            }
            
        }
        
        if(!is_member_Q){
            if(verbose_Q){
                std::cout << "Point is not member of cells set. Returning cell_indexes[0] " << std::endl;
            }
            size_t first_index = *cell_indexes.begin();
            return first_index;
        }
        
        return -1;
    }
    
    // COMPUTE ERRORS UNCOUPLED PROBLEM
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////

    /// Compute L2 and H1 errors for one field approximation
    static void compute_errors_one_field(Mesh & msh, disk::hho_degree_info & hho_di, acoustic_one_field_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<std::vector<double>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
        
        RealType scalar_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        size_t cell_i = 0;
        RealType h = 10.0;
        for (auto& cell : msh)
        {
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof, 0, cell_dof, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, scal_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - scalar_cell_dof;
                scalar_l2_error += diff.dot(mass*diff);

            }

            // flux evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto gr = make_scalar_hho_laplacian(msh, cell, hho_di);
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_dof_data(msh, cell, x_dof);
                Matrix<RealType, Dynamic, 1> recdofs = gr.first * all_dofs;

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    auto t_dphi = rec_basis.eval_gradients( point_pair.point() );
                    Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();

                    for (size_t i = 1; i < t_dphi.rows(); i++){
                        grad_uh = grad_uh + recdofs(i-1)*t_dphi.block(i, 0, 1, 2);
                    }

                    Matrix<RealType, 1, 2> grad_u_exact = Matrix<RealType, 1, 2>::Zero();
                    grad_u_exact(0,0) =  flux_fun(point_pair.point())[0];
                    grad_u_exact(0,1) =  flux_fun(point_pair.point())[1];
                    flux_l2_error += omega * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);

                }
            }

            cell_i++;
        }
        tc.toc();
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();


    }
    
    /// Compute L2 and H1 errors for two fields approximation
    static void compute_errors_two_fields(Mesh & msh, disk::hho_degree_info & hho_di, acoustic_two_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<std::vector<double>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout, bool recons_error_Q = false){

        timecounter tc;
        tc.tic();

        using RealType = double;
        size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_dof = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t cell_dof = n_scal_dof + n_vec_dof;
        
        RealType scalar_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        size_t cell_i = 0;
        RealType h = 10.0;
        for (auto& cell : msh)
        {
             RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            // scalar evaluation
            if(recons_error_Q){
                auto rec_cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto gr = make_scalar_hho_laplacian(msh, cell, hho_di);
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_dof_data(msh, cell, x_dof);
                size_t n_cell_dof = all_dofs.rows() - n_vec_dof;

                Matrix<RealType, Dynamic, 1> cell_dofs = all_dofs.block(n_vec_dof, 0, n_cell_dof, 1);
                Matrix<RealType, Dynamic, 1> rec_scalar_cell_dof = gr.first * cell_dofs;
                
                size_t n_rbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension);
                Matrix<RealType, Dynamic, 1> rec_scalar_dof =  Matrix<RealType, Dynamic, 1>::Zero(n_rbs, 1);
                
                rec_scalar_dof(0, 0) = cell_dofs(0, 0); // The constant is the same.
                rec_scalar_dof.block(1, 0, n_rbs-1, 1) = rec_scalar_cell_dof;
                
                
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, rec_cell_basis, hho_di.reconstruction_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, rec_cell_basis, scal_fun, 1);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - rec_scalar_dof;
                scalar_l2_error += diff.dot(mass*diff);
            }else{
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof+n_vec_dof, 0, n_scal_dof, 1);
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, scal_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                
                Matrix<RealType, Dynamic, 1> diff = real_dofs - scalar_cell_dof;
                scalar_l2_error += diff.dot(mass*diff);
            }

            // flux evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                Matrix<RealType, Dynamic, 1> flux_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_vec_dof, 1);
                
                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    auto t_dphi = cell_basis.eval_gradients( point_pair.point() );
                    
                    Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
                    for (size_t i = 1; i < t_dphi.rows(); i++){
                      grad_uh = grad_uh + flux_cell_dof(i-1)*t_dphi.block(i, 0, 1, 2);
                    }

                    Matrix<RealType, 1, 2> grad_u_exact = Matrix<RealType, 1, 2>::Zero();
                    grad_u_exact(0,0) =  flux_fun(point_pair.point())[0];
                    grad_u_exact(0,1) =  flux_fun(point_pair.point())[1];
                    flux_l2_error += omega * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);

                }
            }

            cell_i++;
        }
        tc.toc();
        
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();
        
    }
    
    /// Compute L2 and H1 errors for one field vectorial approximation
    static void compute_errors_one_field_vectorial(Mesh & msh, disk::hho_degree_info & hho_di, elastodynamic_one_field_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        
        RealType vector_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        RealType h = 10.0;
        size_t cell_ind = 0;
        for (auto& cell : msh)
        {
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, vec_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - vec_cell_dof;
                vector_l2_error += diff.dot(mass*diff);

            }

            elastic_material_data<RealType> material = assembler.get_material_data()[cell_ind];
            RealType mu = material.rho()*material.vs()*material.vs();
            RealType lambda = material.rho()*material.vp()*material.vp() - 2.0*mu;
            
            // flux evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_dof_data(msh, cell, x_dof);
                
                auto           sgr = make_vector_hho_symmetric_laplacian(msh, cell, hho_di);
                disk::dynamic_vector<RealType> GTu = sgr.first * all_dofs;

                auto           dr   = make_hho_divergence_reconstruction(msh, cell, hho_di);
                disk::dynamic_vector<RealType> divu = dr.first * all_dofs;

                auto cbas_v = disk::make_vector_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto cbas_s = disk::make_scalar_monomial_basis(msh, cell, hho_di.face_degree());

                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    
                    auto t_dphi = rec_basis.eval_gradients( point_pair.point() );
                    auto gphi   = cbas_v.eval_sgradients(point_pair.point());
                    auto epsilon = disk::eval(GTu, gphi, dim);
                    auto divphi   = cbas_s.eval_functions(point_pair.point());
                    auto trace_epsilon = disk::eval(divu, divphi);
                    auto sigma = 2.0 * mu * epsilon + lambda * trace_epsilon * disk::static_matrix<RealType, 2, 2>::Identity();
                    auto flux_diff = (flux_fun(point_pair.point()) - sigma).eval();
                    flux_l2_error += omega * flux_diff.squaredNorm();

                }
            }

            cell_ind++;
        }
        tc.toc();
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(vector_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();


    }
    
    /// Compute L2 and H1 errors for two fields vectorial approximation
    static void compute_errors_two_fields_vectorial(Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        size_t cell_dof = n_ten_cbs + n_vec_cbs;
        
        RealType vector_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        size_t cell_i = 0;
        RealType h = 10.0;
        for (auto& cell : msh)
        {
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof.block(cell_i*cell_dof+n_ten_cbs, 0, n_vec_cbs, 1);

            // vector evaluation
            {
                auto cell_basis = disk::make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, vec_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - vec_cell_dof;
                vector_l2_error += diff.dot(mass*diff);

            }

            
            // tensor evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                
                auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    
                    auto t_ten_phi = ten_basis.eval_functions( point_pair.point() );
                    assert(t_ten_phi.size() == ten_basis.size());
                    auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                    
                    auto flux_diff = (flux_fun(point_pair.point()) - sigma_h).eval();
                    flux_l2_error += omega * flux_diff.squaredNorm();

                }
                
            }

            cell_i++;
        }
        tc.toc();
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(vector_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();

    }
    
    /// Compute L2 and H1 errors for three fields vectorial approximation
    static void compute_errors_three_fields_vectorial(Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
        size_t n_sca_cbs = disk::scalar_basis_size(hho_di.face_degree(), dim);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        size_t cell_dof = n_ten_cbs + n_sca_cbs + n_vec_cbs;
        
        RealType vector_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        size_t cell_i = 0;
        RealType h = 10.0;
        for (auto& cell : msh)
        {
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof.block(cell_i*cell_dof+n_ten_cbs+n_sca_cbs, 0, n_vec_cbs, 1);

            // vector evaluation
            {
                auto cell_basis = disk::make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, vec_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - vec_cell_dof;
                vector_l2_error += diff.dot(mass*diff);

            }

            // tensor evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                
                auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);

                auto sca_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.face_degree());
                Matrix<RealType, Dynamic, 1> sigma_v_x_cell_dof = x_dof.block(cell_i*cell_dof+n_ten_cbs, 0, n_sca_cbs, 1);

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    
                    auto t_ten_phi = ten_basis.eval_functions( point_pair.point() );
                    assert(t_ten_phi.size() == ten_basis.size());
                    auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                    
                    auto t_sca_phi = sca_basis.eval_functions( point_pair.point() );
                    assert(t_sca_phi.size() == sca_basis.size());
                    auto sigma_v_h = disk::eval(sigma_v_x_cell_dof, t_sca_phi);
                    
                    sigma_h  += sigma_v_h * disk::static_matrix<RealType, 2, 2>::Identity();
                    
                    auto flux_diff = (flux_fun(point_pair.point()) - sigma_h).eval();
                    flux_l2_error += omega * flux_diff.squaredNorm();

                }
                
            }

            cell_i++;
        }
        tc.toc();
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(vector_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();

    }

    // COMPUTE ENERGY UNCOUPLED PROBLEM
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////

    /// Compute the discrete acoustic energy for one field approximation
    static double compute_acoustic_energy_one_field(Mesh & msh, disk::hho_degree_info & hho_di, acoustic_one_field_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & p_dof, Matrix<double, Dynamic, 1> & v_dof, std::ostream & energy_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
        
        std::vector<RealType> energy_vec(msh.cells_size());
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [&msh,&assembler,&energy_vec,&p_dof,&v_dof,&cell_dof] (size_t & cell_ind){
                            auto& cell = msh.backend_storage()->surfaces[cell_ind];
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> cell_alpha_dof_n_v = v_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                            Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_alpha_dof_n_v;
                            Matrix<RealType, 1, 1> term_1 = cell_alpha_dof_n_v.transpose() * cell_mass_tested;
                            
                            energy_vec[cell_ind] = term_1(0,0);
                
                            Matrix<RealType, Dynamic, Dynamic> laplacian_loc = assembler.laplacian_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> cell_p_dofs = assembler.gather_dof_data(msh, cell, p_dof);
                            Matrix<RealType, Dynamic, 1> cell_stiff_tested = laplacian_loc * cell_p_dofs;
                            Matrix<RealType, 1, 1> term_2 = cell_p_dofs.transpose() * cell_stiff_tested;
                
                            energy_vec[cell_ind] += term_2(0,0);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> cell_alpha_dof_n_v = v_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_alpha_dof_n_v;
                Matrix<RealType, 1, 1> term_1 = cell_alpha_dof_n_v.transpose() * cell_mass_tested;
                
                energy_vec[cell_ind] = term_1(0,0);
    
                Matrix<RealType, Dynamic, Dynamic> laplacian_loc = assembler.laplacian_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> cell_p_dofs = assembler.gather_dof_data(msh, cell, p_dof);
                Matrix<RealType, Dynamic, 1> cell_stiff_tested = laplacian_loc * cell_p_dofs;
                Matrix<RealType, 1, 1> term_2 = cell_p_dofs.transpose() * cell_stiff_tested;
    
                energy_vec[cell_ind] += term_2(0,0);
        
            }
        #endif
    
        RealType energy_h = std::accumulate(energy_vec.begin(), energy_vec.end(),0.0);
        energy_h *= 0.5;
        
        tc.toc();
        std::cout << bold << cyan << "Energy completed: " << tc << " seconds" << reset << std::endl;
        energy_file << time << "   " << std::setprecision(16) << energy_h << std::endl;
        return energy_h;
    }
    
    /// Compute the discrete acoustic energy for one field approximation
    static double compute_acoustic_energy_two_fields(Mesh & msh, disk::hho_degree_info & hho_di, acoustic_two_fields_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & x_dof, std::ostream & energy_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_cbs = n_scal_cbs + n_vec_cbs;
        
        std::vector<RealType> energy_vec(msh.cells_size());
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [&msh,&assembler,&energy_vec,&x_dof,&n_cbs] (size_t & cell_ind){
                            auto& cell = msh.backend_storage()->surfaces[cell_ind];
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> cell_dof = x_dof.block(cell_ind*n_cbs, 0, n_cbs, 1);
                            Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_dof;
                            Matrix<RealType, 1, 1> term = cell_dof.transpose() * cell_mass_tested;
                            energy_vec[cell_ind] = term(0,0);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> cell_dof = x_dof.block(cell_ind*n_cbs, 0, n_cbs, 1);
                Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_dof;
                Matrix<RealType, 1, 1> term = cell_dof.transpose() * cell_mass_tested;
            
                energy_vec[cell_ind] = term(0,0);
        
            }
        #endif
    
        RealType energy_h = std::accumulate(energy_vec.begin(), energy_vec.end(),0.0);
        energy_h *= 0.5;
        
        tc.toc();
        std::cout << bold << cyan << "Energy completed: " << tc << " seconds" << reset << std::endl;
        energy_file << time << "   " << std::setprecision(16) << energy_h << std::endl;
        
        return energy_h;
    }
    
    /// Compute the discrete elastic energy for one field approximation
    static double compute_elastic_energy_one_field(Mesh & msh, disk::hho_degree_info & hho_di, elastodynamic_one_field_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & u_dof, Matrix<double, Dynamic, 1> & v_dof, std::ostream & energy_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        
        std::vector<RealType> energy_vec(msh.cells_size());
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [&msh,&assembler,&energy_vec,&u_dof,&v_dof,&cell_dof] (size_t & cell_ind){
                            auto& cell = msh.backend_storage()->surfaces[cell_ind];
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> cell_alpha_dof_n_v = v_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                            Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_alpha_dof_n_v;
                            Matrix<RealType, 1, 1> term_1 = cell_alpha_dof_n_v.transpose() * cell_mass_tested;
                            
                            energy_vec[cell_ind] = term_1(0,0);

                            Matrix<RealType, Dynamic, Dynamic> laplacian_loc = assembler.laplacian_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> cell_p_dofs = assembler.gather_dof_data(msh, cell, u_dof);
                            Matrix<RealType, Dynamic, 1> cell_stiff_tested = laplacian_loc * cell_p_dofs;
                            Matrix<RealType, 1, 1> term_2 = cell_p_dofs.transpose() * cell_stiff_tested;

                            energy_vec[cell_ind] += term_2(0,0);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> cell_alpha_dof_n_v = v_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                Matrix<RealType, Dynamic, 1> cell_mass_tested = mass_matrix * cell_alpha_dof_n_v;
                Matrix<RealType, 1, 1> term_1 = cell_alpha_dof_n_v.transpose() * cell_mass_tested;
                
                energy_vec[cell_ind] = term_1(0,0);

                Matrix<RealType, Dynamic, Dynamic> laplacian_loc = assembler.laplacian_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> cell_u_dofs = assembler.gather_dof_data(msh, cell, u_dof);
                Matrix<RealType, Dynamic, 1> cell_stiff_tested = laplacian_loc * cell_u_dofs;
                Matrix<RealType, 1, 1> term_2 = cell_u_dofs.transpose() * cell_stiff_tested;

                energy_vec[cell_ind] += term_2(0,0);
        
            }
        #endif
        
        RealType energy_h = std::accumulate(energy_vec.begin(), energy_vec.end(),0.0);
        energy_h *= 0.5;
        
        tc.toc();
        std::cout << bold << cyan << "Energy completed: " << tc << " seconds" << reset << std::endl;
        energy_file << time << "   " << std::setprecision(16) << energy_h << std::endl;
        return energy_h;
    }
    
    /// Compute the discrete elastic energy for two fields approximation
    static double compute_elastic_energy_two_fields(Mesh & msh, disk::hho_degree_info & hho_di, elastodynamic_two_fields_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & x_dof, std::ostream & energy_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        
        std::vector<RealType> energy_vec(msh.cells_size());
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [&msh,&assembler,&energy_vec,&x_dof,&n_ten_cbs,&n_vec_cbs] (size_t & cell_ind){
                            auto& cell = msh.backend_storage()->surfaces[cell_ind];
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_dof_data(msh, cell, x_dof);
                            
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs);
                            Matrix<RealType, Dynamic, 1> v_dof = x_dof_loc.block(n_ten_cbs, 0, n_vec_cbs, 1);
                            Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * v_dof;
                            Matrix<RealType, 1, 1> term_1 = v_dof.transpose() * v_mass_tested;
                            energy_vec[cell_ind] = term_1(0,0);
                            
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix_stress = mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs);
                            Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs, 1);
                            Matrix<RealType, Dynamic, 1> epsilon_mass = mass_matrix_stress * sigma_dof;
                            Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
                            energy_vec[cell_ind] += term_2(0,0);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_dof_data(msh, cell, x_dof);
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs);
                Matrix<RealType, Dynamic, 1> v_dof = x_dof_loc.block(n_ten_cbs, 0, n_vec_cbs, 1);
                Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * v_dof;
                Matrix<RealType, 1, 1> term_1 = v_dof.transpose() * v_mass_tested;
                energy_vec[cell_ind] = term_1(0,0);
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix_stress = mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs);
                Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs, 1);
                Matrix<RealType, Dynamic, 1> epsilon_mass = mass_matrix_stress * sigma_dof;
                Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
                energy_vec[cell_ind] += term_2(0,0);
        
            }
        #endif
    
        RealType energy_h = std::accumulate(energy_vec.begin(), energy_vec.end(),0.0);
        energy_h *= 0.5;
        
        tc.toc();
        std::cout << bold << cyan << "Energy completed: " << tc << " seconds" << reset << std::endl;
        energy_file << time << "   " << std::setprecision(16) << energy_h << std::endl;
        return energy_h;
    }
    
    /// Compute the discrete elastic energy for three field approximation
    static double compute_elastic_energy_three_fields(Mesh & msh, disk::hho_degree_info & hho_di, elastodynamic_three_fields_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & x_dof, std::ostream & energy_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_sca_cbs = disk::scalar_basis_size(hho_di.face_degree(), Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        
        std::vector<RealType> energy_vec(msh.cells_size());
        #ifdef HAVE_INTEL_TBB
                size_t n_cells = msh.cells_size();
                tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                    [&msh,&assembler,&energy_vec,&x_dof,&n_ten_cbs,&n_sca_cbs,&n_vec_cbs] (size_t & cell_ind){
                            auto& cell = msh.backend_storage()->surfaces[cell_ind];
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_dof_data(msh, cell, x_dof);
                            
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs + n_sca_cbs, n_ten_cbs + n_sca_cbs, n_vec_cbs, n_vec_cbs);
                            Matrix<RealType, Dynamic, 1> v_dof = x_dof_loc.block(n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
                            Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * v_dof;
                            Matrix<RealType, 1, 1> term_1 = v_dof.transpose() * v_mass_tested;
                            energy_vec[cell_ind] = term_1(0,0);
                            
                            Matrix<RealType, Dynamic, Dynamic> mass_matrix_stress = mass_matrix.block(0, 0, n_ten_cbs + n_sca_cbs, n_ten_cbs + n_sca_cbs);
                            Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs + n_sca_cbs, 1);
                            Matrix<RealType, Dynamic, 1> epsilon_mass = mass_matrix_stress * sigma_dof;
                            Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
                            energy_vec[cell_ind] += term_2(0,0);
                }
            );
        #else
            for (size_t cell_ind = 0; cell_ind < msh.cells_size(); cell_ind++)
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.mass_operator(cell_ind, msh, cell);
                Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_dof_data(msh, cell, x_dof);
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs + n_sca_cbs, n_ten_cbs + n_sca_cbs, n_vec_cbs, n_vec_cbs);
                Matrix<RealType, Dynamic, 1> v_dof = x_dof_loc.block(n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
                Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * v_dof;
                Matrix<RealType, 1, 1> term_1 = v_dof.transpose() * v_mass_tested;
                energy_vec[cell_ind] = term_1(0,0);
                
                Matrix<RealType, Dynamic, Dynamic> mass_matrix_stress = mass_matrix.block(0, 0, n_ten_cbs + n_sca_cbs, n_ten_cbs + n_sca_cbs);
                Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs + n_sca_cbs, 1);
                Matrix<RealType, Dynamic, 1> epsilon_mass = mass_matrix_stress * sigma_dof;
                Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
                energy_vec[cell_ind] += term_2(0,0);
        
            }
        #endif
    
        RealType energy_h = std::accumulate(energy_vec.begin(), energy_vec.end(),0.0);
        energy_h *= 0.5;
        
        tc.toc();
        std::cout << bold << cyan << "Energy completed: " << tc << " seconds" << reset << std::endl;
        energy_file << time << "   " << std::setprecision(16) << energy_h << std::endl;
        return energy_h;
    }
    
    // SENSORS UNCOUPLED PROBELM
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////
    
    /// Record data at provided point for one field approximation
    static void record_data_acoustic_one_field(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);

        RealType vh = 0.0;

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vh\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // scalar evaluation
            {
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                auto t_phi = cell_basis.eval_functions( pt );
                vh = scalar_cell_dof.dot( t_phi );
            }
            
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh << std::endl;
        seismogram_file.flush();

    }
    
    /// Record velocity data at provided point for one field approximation
    static void record_velocity_data_acoustic_one_field(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, acoustic_one_field_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);

        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // reconstructed velocity evaluation
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto gr = make_scalar_hho_laplacian(msh, cell, hho_di);
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_dof_data(msh, cell, x_dof);
                Matrix<RealType, Dynamic, 1> recdofs = gr.first * all_dofs;
                
                auto t_dphi = rec_basis.eval_gradients( pt );
                Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();

                for (size_t i = 1; i < t_dphi.rows(); i++){
                    grad_uh = grad_uh + recdofs(i-1)*t_dphi.block(i, 0, 1, 2);
                }

                RealType rho = assembler.get_material_data()[cell_ind].rho();
                vh = (1.0/rho)*(grad_uh);
            }
            
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }
    
    /// Record data at provided point for two fields approximation
    static void record_data_acoustic_two_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_dof = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t cell_dof = n_scal_dof + n_vec_dof;

        RealType vh = 0.0;

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vh\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // scalar evaluation
            {
                
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_ind*cell_dof+n_vec_dof, 0, n_scal_dof, 1);
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                auto t_phi = cell_basis.eval_functions( pt );
                vh = scalar_cell_dof.dot( t_phi );
            }
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh << std::endl;
        seismogram_file.flush();

    }
    
    /// Record data at provided point for one field vectorial approximation
    static void record_data_elastic_one_field(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);

        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // vector evaluation
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_ind*cell_dof, 0, cell_dof, 1);
                auto t_phi = cell_basis.eval_functions( pt );
                assert(t_phi.rows() == cell_basis.size());
                vh = disk::eval(vec_x_cell_dof, t_phi);
            }
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }
    
    /// Record data at provided point for two fields vectorial approximation
    static void record_data_elastic_two_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        size_t cell_dof = n_ten_cbs + n_vec_cbs;

        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // vector evaluation
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_ind*cell_dof + n_ten_cbs , 0, n_vec_cbs, 1);
                auto t_phi = cell_basis.eval_functions( pt );
                assert(t_phi.rows() == cell_basis.size());
                vh = disk::eval(vec_x_cell_dof, t_phi);
            }
            
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }
    
    /// Record data at provided point for three fields vectorial approximation
    static void record_data_elastic_three_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
          size_t n_sca_cbs = disk::scalar_basis_size(hho_di.face_degree(), dim);
          size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
          size_t cell_dof = n_ten_cbs + n_sca_cbs + n_vec_cbs;

        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            // vector evaluation
            {
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_ind*cell_dof + n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
                auto t_phi = cell_basis.eval_functions( pt );
                assert(t_phi.rows() == cell_basis.size());
                vh = disk::eval(vec_x_cell_dof, t_phi);
            }
            
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }
    
    // WRITE SILO UNCOUPLED PROBLEM
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////

    // Write a silo file for one field approximation
    static void write_silo_one_field(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof,
    std::function<double(const typename Mesh::point_type& )> scal_fun, bool cell_centered_Q){
    
            timecounter tc;
            tc.tic();
            
            auto dim = Mesh::dimension;
            auto num_cells = msh.cells_size();
            auto num_points = msh.points_size();
            using RealType = double;
            std::vector<RealType> exact_u, approx_u;
            size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
            
            if (cell_centered_Q) {
                exact_u.reserve( num_cells );
                approx_u.reserve( num_cells );
    
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto bar = barycenter(msh, cell);
                    exact_u.push_back( scal_fun(bar) );
                    
                    // scalar evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof, 0, cell_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.push_back(uh);
                    }
                    cell_i++;
                }
    
            }else{
    
                exact_u.reserve( num_points );
                approx_u.reserve( num_points );
    
                // scan for selected cells, common cells are discardable
                std::map<size_t, size_t> point_to_cell;
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto points = cell.point_ids();
                    size_t n_p = points.size();
                    for (size_t l = 0; l < n_p; l++)
                    {
                        auto pt_id = points[l];
                        point_to_cell[pt_id] = cell_i;
                    }
                    cell_i++;
                }
    
                for (auto& pt_id : point_to_cell)
                {
                    auto bar = *std::next(msh.points_begin(), pt_id.first);
                    exact_u.push_back( scal_fun(bar) );
    
                    cell_i = pt_id.second;
                    auto cell = *std::next(msh.cells_begin(), cell_i);
                    // scalar evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof, 0, cell_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.push_back(uh);
                    }
                }
    
            }
    
            disk::silo_database silo;
            silo_file_name += std::to_string(it) + ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            if (cell_centered_Q) {
                disk::silo_zonal_variable<double> v_silo("v", exact_u);
                disk::silo_zonal_variable<double> vh_silo("vh", approx_u);
                silo.add_variable("mesh", v_silo);
                silo.add_variable("mesh", vh_silo);
            }else{
                disk::silo_nodal_variable<double> v_silo("v", exact_u);
                disk::silo_nodal_variable<double> vh_silo("vh", approx_u);
                silo.add_variable("mesh", v_silo);
                silo.add_variable("mesh", vh_silo);
            }
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file for two fields approximation
    static void write_silo_two_fields(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof,
    std::function<double(const typename Mesh::point_type& )> scal_fun,  std::function<std::vector<double>(const typename Mesh::point_type& )> flux_fun, bool cell_centered_Q){
    
            timecounter tc;
            tc.tic();
            
            auto num_cells = msh.cells_size();
            auto num_points = msh.points_size();
            using RealType = double;
            std::vector<RealType> exact_u, approx_u;
            std::vector<RealType> exact_dux, exact_duy, approx_dux, approx_duy;
            
            size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
            size_t n_vec_dof = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
            size_t cell_dof = n_scal_dof + n_vec_dof;
            
            if (cell_centered_Q) {
                exact_u.reserve( num_cells );
                approx_u.reserve( num_cells );
                exact_dux.reserve( num_cells );
                exact_duy.reserve( num_cells );
                approx_dux.reserve( num_cells );
                approx_duy.reserve( num_cells );
    
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto bar = barycenter(msh, cell);
                    exact_u.push_back( scal_fun(bar) );
                    exact_dux.push_back( flux_fun(bar)[0] );
                    exact_duy.push_back( flux_fun(bar)[1] );
                    
                    // scalar evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof+n_vec_dof, 0, n_scal_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.push_back(uh);
                    }
                    
                    // flux evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                        Matrix<RealType, Dynamic, 1> flux_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_vec_dof, 1);
                        auto t_dphi = cell_basis.eval_gradients( bar );
                        
                        Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
                        for (size_t i = 1; i < t_dphi.rows(); i++){
                          grad_uh = grad_uh + flux_cell_dof(i-1)*t_dphi.block(i, 0, 1, 2);
                        }
    
                        approx_dux.push_back(grad_uh(0,0));
                        approx_duy.push_back(grad_uh(0,1));
                    }
                    
                    cell_i++;
                }
    
            }else{
    
                exact_u.reserve( num_points );
                approx_u.reserve( num_points );
                exact_dux.reserve( num_points );
                exact_duy.reserve( num_points );
                approx_dux.reserve( num_points );
                approx_duy.reserve( num_points );
    
                // scan for selected cells, common cells are discardable
                std::map<size_t, size_t> point_to_cell;
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto points = cell.point_ids();
                    size_t n_p = points.size();
                    for (size_t l = 0; l < n_p; l++)
                    {
                        auto pt_id = points[l];
                        point_to_cell[pt_id] = cell_i;
                    }
                    cell_i++;
                }
    
                for (auto& pt_id : point_to_cell)
                {
                    auto bar = *std::next(msh.points_begin(), pt_id.first);
                    cell_i = pt_id.second;
                    auto cell = *std::next(msh.cells_begin(), cell_i);
                    
                    exact_u.push_back( scal_fun(bar) );
                    exact_dux.push_back( flux_fun(bar)[0] );
                    exact_duy.push_back( flux_fun(bar)[1] );
                    
                    // scalar evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_i*cell_dof+n_vec_dof, 0, n_scal_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.push_back(uh);
                    }
                    
                    // flux evaluation
                    {
                        auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                        Matrix<RealType, Dynamic, 1> flux_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_vec_dof, 1);
                        auto t_dphi = cell_basis.eval_gradients( bar );
                        
                        Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
                        for (size_t i = 1; i < t_dphi.rows(); i++){
                          grad_uh = grad_uh + flux_cell_dof(i-1)*t_dphi.block(i, 0, 1, 2);
                        }
    
                        approx_dux.push_back(grad_uh(0,0));
                        approx_duy.push_back(grad_uh(0,1));
                    }
                }
    
            }
    
            disk::silo_database silo;
            silo_file_name += std::to_string(it) + ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            if (cell_centered_Q) {
                disk::silo_zonal_variable<double> v_silo("v", exact_u);
                disk::silo_zonal_variable<double> qx_silo("qx", exact_dux);
                disk::silo_zonal_variable<double> qy_silo("qy", exact_duy);
                disk::silo_zonal_variable<double> vh_silo("vh", approx_u);
                disk::silo_zonal_variable<double> qhx_silo("qhx", approx_dux);
                disk::silo_zonal_variable<double> qhy_silo("qhy", approx_duy);
                silo.add_variable("mesh", v_silo);
                silo.add_variable("mesh", qx_silo);
                silo.add_variable("mesh", qy_silo);
                silo.add_variable("mesh", vh_silo);
                silo.add_variable("mesh", qhx_silo);
                silo.add_variable("mesh", qhy_silo);
            }else{
                disk::silo_nodal_variable<double> v_silo("v", exact_u);
                disk::silo_nodal_variable<double> qx_silo("qx", exact_dux);
                disk::silo_nodal_variable<double> qy_silo("qy", exact_duy);
                disk::silo_nodal_variable<double> vh_silo("vh", approx_u);
                disk::silo_nodal_variable<double> qhx_silo("qhx", approx_dux);
                disk::silo_nodal_variable<double> qhy_silo("qhy", approx_duy);
                silo.add_variable("mesh", v_silo);
                silo.add_variable("mesh", qx_silo);
                silo.add_variable("mesh", qy_silo);
                silo.add_variable("mesh", vh_silo);
                silo.add_variable("mesh", qhx_silo);
                silo.add_variable("mesh", qhy_silo);
            }
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file for one field vectorial approximation
    static void write_silo_one_field_vectorial(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof,
    std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, bool cell_centered_Q){
    
            timecounter tc;
            tc.tic();
            
            auto dim = Mesh::dimension;
            auto num_cells = msh.cells_size();
            auto num_points = msh.points_size();
            using RealType = double;
            std::vector<RealType> exact_ux, exact_uy, approx_ux, approx_uy;
            size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
            
            if (cell_centered_Q) {
                exact_ux.reserve( num_cells );
                exact_uy.reserve( num_cells );
                approx_ux.reserve( num_cells );
                approx_uy.reserve( num_cells );
    
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto bar = barycenter(msh, cell);
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, cell_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                    cell_i++;
                }
    
            }else{
    
                exact_ux.reserve( num_points );
                exact_uy.reserve( num_points );
                approx_ux.reserve( num_points );
                approx_uy.reserve( num_points );
                
                // scan for selected cells, common cells are discardable
                std::map<size_t, size_t> point_to_cell;
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto points = cell.point_ids();
                    size_t n_p = points.size();
                    for (size_t l = 0; l < n_p; l++)
                    {
                        auto pt_id = points[l];
                        point_to_cell[pt_id] = cell_i;
                    }
                    cell_i++;
                }
    
                for (auto& pt_id : point_to_cell)
                {
                    auto bar = *std::next(msh.points_begin(), pt_id.first);
                    cell_i = pt_id.second;
                    auto cell = *std::next(msh.cells_begin(), cell_i);
                    
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, cell_dof, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                }
    
            }
    
            disk::silo_database silo;
            silo_file_name += std::to_string(it) + ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            if (cell_centered_Q) {
                disk::silo_zonal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_zonal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_zonal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_zonal_variable<double> vhy_silo("vhy", approx_uy);
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
            }else{
                disk::silo_nodal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_nodal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_nodal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_nodal_variable<double> vhy_silo("vhy", approx_uy);
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
            }
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file for two fields vectorial approximation
    static void write_silo_two_fields_vectorial(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof,std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> flux_fun, bool cell_centered_Q){
    
            timecounter tc;
            tc.tic();
            
            auto dim = Mesh::dimension;
            auto num_cells = msh.cells_size();
            auto num_points = msh.points_size();
            using RealType = double;
            std::vector<RealType> exact_ux, exact_uy, approx_ux, approx_uy;
            std::vector<RealType> exact_sxx, exact_sxy, exact_syy;
            std::vector<RealType> approx_sxx, approx_sxy, approx_syy;
            size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
            size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
            size_t cell_dof = n_ten_cbs + n_vec_cbs;
            
            if (cell_centered_Q) {
                exact_ux.reserve( num_cells );
                exact_uy.reserve( num_cells );
                approx_ux.reserve( num_cells );
                approx_uy.reserve( num_cells );
                
                exact_sxx.reserve( num_cells );
                exact_sxy.reserve( num_cells );
                exact_syy.reserve( num_cells );
                approx_sxx.reserve( num_cells );
                approx_sxy.reserve( num_cells );
                approx_syy.reserve( num_cells );
    
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto bar = barycenter(msh, cell);
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    exact_sxx.push_back( flux_fun(bar)(0,0) );
                    exact_sxy.push_back( flux_fun(bar)(0,1) );
                    exact_syy.push_back( flux_fun(bar)(1,1) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof + n_ten_cbs, 0, n_vec_cbs, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                    
                    // tensor evaluation
                    {
                        auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                        Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);
    
                        auto t_ten_phi = ten_basis.eval_functions( bar );
                        assert(t_ten_phi.size() == ten_basis.size());
                        auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                        
                        approx_sxx.push_back(sigma_h(0,0));
                        approx_sxy.push_back(sigma_h(0,1));
                        approx_syy.push_back(sigma_h(1,1));
                    }
                    
                    cell_i++;
                }
    
            }else{
    
                exact_ux.reserve( num_points );
                exact_uy.reserve( num_points );
                approx_ux.reserve( num_points );
                approx_uy.reserve( num_points );
                
                exact_sxx.reserve( num_points );
                exact_sxy.reserve( num_points );
                exact_syy.reserve( num_points );
                approx_sxx.reserve( num_points );
                approx_sxy.reserve( num_points );
                approx_syy.reserve( num_points );
                
                // scan for selected cells, common cells are discardable
                std::map<size_t, size_t> point_to_cell;
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto points = cell.point_ids();
                    size_t n_p = points.size();
                    for (size_t l = 0; l < n_p; l++)
                    {
                        auto pt_id = points[l];
                        point_to_cell[pt_id] = cell_i;
                    }
                    cell_i++;
                }
    
                for (auto& pt_id : point_to_cell)
                {
                    auto bar = *std::next(msh.points_begin(), pt_id.first);
                    cell_i = pt_id.second;
                    auto cell = *std::next(msh.cells_begin(), cell_i);
                    
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    exact_sxx.push_back( flux_fun(bar)(0,0) );
                    exact_sxy.push_back( flux_fun(bar)(0,1) );
                    exact_syy.push_back( flux_fun(bar)(1,1) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof + n_ten_cbs, 0, n_vec_cbs, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                    
                    // tensor evaluation
                    {
                        auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                        Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);
    
                        auto t_ten_phi = ten_basis.eval_functions( bar );
                        assert(t_ten_phi.size() == ten_basis.size());
                        auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                        
                        approx_sxx.push_back(sigma_h(0,0));
                        approx_sxy.push_back(sigma_h(0,1));
                        approx_syy.push_back(sigma_h(1,1));
                    }
                    
                }
    
            }
    
            disk::silo_database silo;
            silo_file_name += std::to_string(it) + ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            if (cell_centered_Q) {
                disk::silo_zonal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_zonal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_zonal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_zonal_variable<double> vhy_silo("vhy", approx_uy);
                
                disk::silo_nodal_variable<double> sxx_silo("sxx", exact_sxx);
                disk::silo_nodal_variable<double> sxy_silo("sxy", exact_sxy);
                disk::silo_nodal_variable<double> syy_silo("syy", exact_syy);
                disk::silo_nodal_variable<double> shxx_silo("shxx", approx_sxx);
                disk::silo_nodal_variable<double> shxy_silo("shxy", approx_sxy);
                disk::silo_nodal_variable<double> shyy_silo("shyy", approx_syy);
                
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
                silo.add_variable("mesh", sxx_silo);
                silo.add_variable("mesh", sxy_silo);
                silo.add_variable("mesh", syy_silo);
                silo.add_variable("mesh", shxx_silo);
                silo.add_variable("mesh", shxy_silo);
                silo.add_variable("mesh", shyy_silo);
            }else{
                disk::silo_nodal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_nodal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_nodal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_nodal_variable<double> vhy_silo("vhy", approx_uy);
                
                disk::silo_nodal_variable<double> sxx_silo("sxx", exact_sxx);
                disk::silo_nodal_variable<double> sxy_silo("sxy", exact_sxy);
                disk::silo_nodal_variable<double> syy_silo("syy", exact_syy);
                disk::silo_nodal_variable<double> shxx_silo("shxx", approx_sxx);
                disk::silo_nodal_variable<double> shxy_silo("shxy", approx_sxy);
                disk::silo_nodal_variable<double> shyy_silo("shyy", approx_syy);
                
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
                silo.add_variable("mesh", sxx_silo);
                silo.add_variable("mesh", sxy_silo);
                silo.add_variable("mesh", syy_silo);
                silo.add_variable("mesh", shxx_silo);
                silo.add_variable("mesh", shxy_silo);
                silo.add_variable("mesh", shyy_silo);
                
            }
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file for three fields vectorial approximation
    static void write_silo_three_fields_vectorial(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof,std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> flux_fun, bool cell_centered_Q){
    
            timecounter tc;
            tc.tic();
            
            auto dim = Mesh::dimension;
            auto num_cells = msh.cells_size();
            auto num_points = msh.points_size();
            using RealType = double;
            std::vector<RealType> exact_ux, exact_uy, approx_ux, approx_uy;
            std::vector<RealType> exact_sxx, exact_sxy, exact_syy;
            std::vector<RealType> approx_sxx, approx_sxy, approx_syy;
            size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), dim, dim);
            size_t n_sca_cbs = disk::scalar_basis_size(hho_di.face_degree(), dim);
            size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
            size_t cell_dof = n_ten_cbs + n_sca_cbs + n_vec_cbs;
            
            if (cell_centered_Q) {
                exact_ux.reserve( num_cells );
                exact_uy.reserve( num_cells );
                approx_ux.reserve( num_cells );
                approx_uy.reserve( num_cells );
                
                exact_sxx.reserve( num_cells );
                exact_sxy.reserve( num_cells );
                exact_syy.reserve( num_cells );
                approx_sxx.reserve( num_cells );
                approx_sxy.reserve( num_cells );
                approx_syy.reserve( num_cells );
    
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto bar = barycenter(msh, cell);
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    exact_sxx.push_back( flux_fun(bar)(0,0) );
                    exact_sxy.push_back( flux_fun(bar)(0,1) );
                    exact_syy.push_back( flux_fun(bar)(1,1) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof + n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                    
                    // tensor evaluation
                    {
                        auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                        Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);
    
                        auto sca_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.face_degree());
                        Matrix<RealType, Dynamic, 1> sigma_v_x_cell_dof = x_dof.block(cell_i*cell_dof+n_ten_cbs, 0, n_sca_cbs, 1);
    
                        auto t_ten_phi = ten_basis.eval_functions( bar );
                        assert(t_ten_phi.size() == ten_basis.size());
                        auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                        
                        auto t_sca_phi = sca_basis.eval_functions( bar );
                        assert(t_sca_phi.size() == sca_basis.size());
                        auto sigma_v_h = disk::eval(sigma_v_x_cell_dof, t_sca_phi);
                    
                        sigma_h  += sigma_v_h * disk::static_matrix<RealType, 2, 2>::Identity();
    
                        approx_sxx.push_back(sigma_h(0,0));
                        approx_sxy.push_back(sigma_h(0,1));
                        approx_syy.push_back(sigma_h(1,1));
                    }
                    
                    cell_i++;
                }
    
            }else{
    
                exact_ux.reserve( num_points );
                exact_uy.reserve( num_points );
                approx_ux.reserve( num_points );
                approx_uy.reserve( num_points );
                
                exact_sxx.reserve( num_points );
                exact_sxy.reserve( num_points );
                exact_syy.reserve( num_points );
                approx_sxx.reserve( num_points );
                approx_sxy.reserve( num_points );
                approx_syy.reserve( num_points );
                
                // scan for selected cells, common cells are discardable
                std::map<size_t, size_t> point_to_cell;
                size_t cell_i = 0;
                for (auto& cell : msh)
                {
                    auto points = cell.point_ids();
                    size_t n_p = points.size();
                    for (size_t l = 0; l < n_p; l++)
                    {
                        auto pt_id = points[l];
                        point_to_cell[pt_id] = cell_i;
                    }
                    cell_i++;
                }
    
                for (auto& pt_id : point_to_cell)
                {
                    auto bar = *std::next(msh.points_begin(), pt_id.first);
                    cell_i = pt_id.second;
                    auto cell = *std::next(msh.cells_begin(), cell_i);
                    
                    exact_ux.push_back( vec_fun(bar)(0,0) );
                    exact_uy.push_back( vec_fun(bar)(1,0) );
                    
                    exact_sxx.push_back( flux_fun(bar)(0,0) );
                    exact_sxy.push_back( flux_fun(bar)(0,1) );
                    exact_syy.push_back( flux_fun(bar)(1,1) );
                    
                    // vector evaluation
                    {
                        auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                        Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(cell_i*cell_dof + n_ten_cbs + n_sca_cbs, 0, n_vec_cbs, 1);
                        auto t_phi = cell_basis.eval_functions( bar );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.push_back(uh(0,0));
                        approx_uy.push_back(uh(1,0));
                    }
                    
                    // tensor evaluation
                    {
                        auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                        Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(cell_i*cell_dof, 0, n_ten_cbs, 1);
    
                        auto sca_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.face_degree());
                        Matrix<RealType, Dynamic, 1> sigma_v_x_cell_dof = x_dof.block(cell_i*cell_dof+n_ten_cbs, 0, n_sca_cbs, 1);
    
                        auto t_ten_phi = ten_basis.eval_functions( bar );
                        assert(t_ten_phi.size() == ten_basis.size());
                        auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                        
                        auto t_sca_phi = sca_basis.eval_functions( bar );
                        assert(t_sca_phi.size() == sca_basis.size());
                        auto sigma_v_h = disk::eval(sigma_v_x_cell_dof, t_sca_phi);
                    
                        sigma_h  += sigma_v_h * disk::static_matrix<RealType, 2, 2>::Identity();
    
                        approx_sxx.push_back(sigma_h(0,0));
                        approx_sxy.push_back(sigma_h(0,1));
                        approx_syy.push_back(sigma_h(1,1));
                    }
                    
                }
    
            }
    
            disk::silo_database silo;
            silo_file_name += std::to_string(it) + ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            if (cell_centered_Q) {
                disk::silo_zonal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_zonal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_zonal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_zonal_variable<double> vhy_silo("vhy", approx_uy);
                
                disk::silo_nodal_variable<double> sxx_silo("sxx", exact_sxx);
                disk::silo_nodal_variable<double> sxy_silo("sxy", exact_sxy);
                disk::silo_nodal_variable<double> syy_silo("syy", exact_syy);
                disk::silo_nodal_variable<double> shxx_silo("shxx", approx_sxx);
                disk::silo_nodal_variable<double> shxy_silo("shxy", approx_sxy);
                disk::silo_nodal_variable<double> shyy_silo("shyy", approx_syy);
                
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
                silo.add_variable("mesh", sxx_silo);
                silo.add_variable("mesh", sxy_silo);
                silo.add_variable("mesh", syy_silo);
                silo.add_variable("mesh", shxx_silo);
                silo.add_variable("mesh", shxy_silo);
                silo.add_variable("mesh", shyy_silo);
            }else{
                disk::silo_nodal_variable<double> vx_silo("vx", exact_ux);
                disk::silo_nodal_variable<double> vy_silo("vy", exact_uy);
                disk::silo_nodal_variable<double> vhx_silo("vhx", approx_ux);
                disk::silo_nodal_variable<double> vhy_silo("vhy", approx_uy);
                
                disk::silo_nodal_variable<double> sxx_silo("sxx", exact_sxx);
                disk::silo_nodal_variable<double> sxy_silo("sxy", exact_sxy);
                disk::silo_nodal_variable<double> syy_silo("syy", exact_syy);
                disk::silo_nodal_variable<double> shxx_silo("shxx", approx_sxx);
                disk::silo_nodal_variable<double> shxy_silo("shxy", approx_sxy);
                disk::silo_nodal_variable<double> shyy_silo("shyy", approx_syy);
                
                silo.add_variable("mesh", vx_silo);
                silo.add_variable("mesh", vy_silo);
                silo.add_variable("mesh", vhx_silo);
                silo.add_variable("mesh", vhy_silo);
                silo.add_variable("mesh", sxx_silo);
                silo.add_variable("mesh", sxy_silo);
                silo.add_variable("mesh", syy_silo);
                silo.add_variable("mesh", shxx_silo);
                silo.add_variable("mesh", shxy_silo);
                silo.add_variable("mesh", shyy_silo);
                
            }
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file with acoustic properties as zonal variables
    static void write_silo_acoustic_property_map(std::string silo_file_name, Mesh & msh, std::vector< acoustic_material_data<double> > & material){
    
            timecounter tc;
            tc.tic();
            
            auto num_cells = msh.cells_size();
            std::vector<double> rho_data, vp_data;
            vp_data.reserve( num_cells );
            rho_data.reserve( num_cells );
            
            for (size_t cell_id = 0; cell_id < num_cells; cell_id++)
            {
                acoustic_material_data<double> acoustic_data = material[cell_id];
                rho_data.push_back( acoustic_data.rho() );
                vp_data.push_back( acoustic_data.vp() );
            }
    
            disk::silo_database silo;
            silo_file_name += ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            disk::silo_zonal_variable<double> rho_silo("rho", rho_data);
            disk::silo_zonal_variable<double> vp_silo("vp", vp_data);
            silo.add_variable("mesh", rho_silo);
            silo.add_variable("mesh", vp_silo);
    
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Properties file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // Write a silo file with elastic properties as zonal variables
    static void write_silo_elastic_property_map(std::string silo_file_name, Mesh & msh, std::vector< elastic_material_data<double> > & material){
            
            timecounter tc;
            tc.tic();
            
            auto num_cells = msh.cells_size();
            std::vector<double> rho_data, vp_data, vs_data;
            rho_data.reserve( num_cells );
            vp_data.reserve( num_cells );
            vs_data.reserve( num_cells );
            
            for (size_t cell_id = 0; cell_id < num_cells; cell_id++) {
                elastic_material_data<double> elastic_data = material[cell_id];
                rho_data.push_back( elastic_data.rho() );
                vp_data.push_back( elastic_data.vp() );
                vs_data.push_back( elastic_data.vs() );
            }
            
            disk::silo_database silo;
            silo_file_name += ".silo";
            silo.create(silo_file_name.c_str());
            silo.add_mesh(msh, "mesh");
            disk::silo_zonal_variable<double> rho_silo("rho", rho_data);
            disk::silo_zonal_variable<double> vp_silo("vp", vp_data);
            disk::silo_zonal_variable<double> vs_silo("vs", vs_data);
            silo.add_variable("mesh", rho_silo);
            silo.add_variable("mesh", vp_silo);
            silo.add_variable("mesh", vs_silo);
            
            silo.close();
            tc.toc();
            std::cout << std::endl;
            std::cout << bold << cyan << "Properties file rendered in : " << tc << " seconds" << reset << std::endl;
        }
        
    // COUPLING POSTPROCESS UNCOUPLED PROBLEM
    ////////////////////////////////////////////////////////////////////////////////      
    ////////////////////////////////////////////////////////////////////////////////

    /// Compute errors
    static void compute_errors_two_fields_elastoacoustic(Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_two_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> sigma_fun,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        RealType scalar_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        RealType vector_l2_error = 0.0;
        RealType sigma_l2_error = 0.0;
        RealType h = 10.0;
        size_t e_cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        size_t a_cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
        auto storage = msh.backend_storage();
        
        size_t e_cell_ind = 0;
        for (auto& e_chunk : assembler.get_e_material_data()) {
            
            auto& cell = storage->surfaces[e_chunk.first];
            
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof.block(e_cell_ind*e_cell_dof, 0, e_cell_dof, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, vec_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - vec_cell_dof;
                vector_l2_error += diff.dot(mass*diff);

            }

            elastic_material_data<RealType> material = e_chunk.second;
            RealType mu = material.rho()*material.vs()*material.vs();
            RealType lambda = material.rho()*material.vp()*material.vp() - 2.0*mu;
            
            // sigma evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_e_dof_data(e_cell_ind, msh, cell, x_dof);
                
                auto           sgr = make_vector_hho_symmetric_laplacian(msh, cell, hho_di);
                disk::dynamic_vector<RealType> GTu = sgr.first * all_dofs;

                auto           dr   = make_hho_divergence_reconstruction(msh, cell, hho_di);
                disk::dynamic_vector<RealType> divu = dr.first * all_dofs;

                auto cbas_v = disk::make_vector_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto cbas_s = disk::make_scalar_monomial_basis(msh, cell, hho_di.face_degree());

                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    
                    auto t_dphi = rec_basis.eval_gradients( point_pair.point() );
                    auto gphi   = cbas_v.eval_sgradients(point_pair.point());
                    auto epsilon = disk::eval(GTu, gphi, dim);
                    auto divphi   = cbas_s.eval_functions(point_pair.point());
                    auto trace_epsilon = disk::eval(divu, divphi);
                    auto sigma = 2.0 * mu * epsilon + lambda * trace_epsilon * disk::static_matrix<RealType, 2, 2>::Identity();
                    auto flux_diff = (sigma_fun(point_pair.point()) - sigma).eval();
                    sigma_l2_error += omega * flux_diff.squaredNorm();

                }
            }

            e_cell_ind++;
        }
        
        size_t n_elastic_cell_dof = assembler.get_e_material_data().size() * e_cell_dof;
        size_t a_cell_ind = 0;
        for (auto& a_chunk : assembler.get_a_material_data()) {
            
            auto& cell = storage->surfaces[a_chunk.first];
            
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof+n_elastic_cell_dof, 0, a_cell_dof, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, scal_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - scalar_cell_dof;
                scalar_l2_error += diff.dot(mass*diff);

            }

            // flux evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto gr = make_scalar_hho_laplacian(msh, cell, hho_di);
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_a_dof_data(a_cell_ind, msh, cell, x_dof);
                Matrix<RealType, Dynamic, 1> recdofs = gr.first * all_dofs;

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    auto t_dphi = rec_basis.eval_gradients( point_pair.point() );
                    Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();

                    for (size_t i = 1; i < t_dphi.rows(); i++){
                        grad_uh = grad_uh + recdofs(i-1)*t_dphi.block(i, 0, 1, 2);
                    }

                    Matrix<RealType, 1, 2> grad_u_exact = flux_fun(point_pair.point());
                    flux_l2_error += omega * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);

                }
            }
            a_cell_ind++;
        }
        
        tc.toc();
        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << h << std::endl;
        error_file << "Elastic region : " << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(vector_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(sigma_l2_error) << std::endl;
        error_file << "Acoustic region : " << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();
    }
    
    static void compute_errors_four_fields_elastoacoustic(Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> sigma_fun,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        RealType scalar_l2_error = 0.0;
        RealType flux_l2_error = 0.0;
        RealType vector_l2_error = 0.0;
        RealType sigma_l2_error = 0.0;
        RealType h = 10.0;

        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t e_cell_dof = n_ten_cbs + n_vec_cbs;

        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_cell_dof = n_vel_scal_cbs + n_scal_cbs;
        
        auto storage = msh.backend_storage();
        size_t e_cell_ind = 0;
        for (auto& e_chunk : assembler.get_e_material_data()) {
            auto& cell = storage->surfaces[e_chunk.first];
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }
            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof.block(e_cell_ind*e_cell_dof+n_ten_cbs, 0, n_vec_cbs, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, vec_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - vec_cell_dof;
                vector_l2_error += diff.dot(mass*diff);

            }

            elastic_material_data<RealType> material = e_chunk.second;
            RealType mu = material.rho()*material.vs()*material.vs();
            RealType lambda = material.rho()*material.vp()*material.vp() - 2.0*mu;
            
            
            // tensor evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                
                auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
                Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof, 0, n_ten_cbs, 1);

                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    
                    auto t_ten_phi = ten_basis.eval_functions( point_pair.point() );
                    assert(t_ten_phi.size() == ten_basis.size());
                    auto sigma_h = disk::eval(ten_x_cell_dof, t_ten_phi);
                    
                    auto flux_diff = (sigma_fun(point_pair.point()) - sigma_h).eval();
                    sigma_l2_error += omega * flux_diff.squaredNorm();

                }
                
            }
            
            e_cell_ind++;
        }
        
        size_t n_elastic_cell_dof = assembler.get_e_material_data().size() * e_cell_dof;
        size_t a_cell_ind = 0;
        for (auto& a_chunk : assembler.get_a_material_data()) {
            
            auto& cell = storage->surfaces[a_chunk.first];
            
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof+n_elastic_cell_dof+n_vel_scal_cbs, 0, n_scal_cbs, 1);

            // scalar evaluation
            {
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, scal_fun);
                Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                Matrix<RealType, Dynamic, 1> diff = real_dofs - scalar_cell_dof;
                scalar_l2_error += diff.dot(mass*diff);

            }

            // flux evaluation
            {
                auto int_rule = integrate(msh, cell, 2*(hho_di.cell_degree()+1));
                auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                Matrix<RealType, Dynamic, 1> flux_cell_dof = x_dof.block(a_cell_ind*a_cell_dof+n_elastic_cell_dof, 0, n_vel_scal_cbs, 1);
                
                // Error integrals
                for (auto & point_pair : int_rule) {

                    RealType omega = point_pair.weight();
                    auto t_dphi = cell_basis.eval_gradients( point_pair.point() );
                    
                    Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
                    for (size_t i = 1; i < t_dphi.rows(); i++){
                      grad_uh = grad_uh + flux_cell_dof(i-1)*t_dphi.block(i, 0, 1, 2);
                    }

                    Matrix<RealType, 1, 2> grad_u_exact = Matrix<RealType, 1, 2>::Zero();
                    grad_u_exact(0,0) =  flux_fun(point_pair.point())[0];
                    grad_u_exact(0,1) =  flux_fun(point_pair.point())[1];
                    flux_l2_error += omega * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);

                }
            }
            a_cell_ind++;
        }
        
        tc.toc();
        // std::cout << bold << red << "   ERROR COMPLETED: " << tc << " seconds" << reset << std::endl;
        error_file << "Characteristic h size = " << h << std::endl;
        error_file << "Elastic region : " << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(vector_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(sigma_l2_error) << std::endl;
        error_file << "Acoustic region : " << std::endl;
        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();
    }
    
    static void compute_errors_four_fields_elastoacoustic_energy_norm(Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> vec_fun, std::function<disk::static_matrix<double, 2, 2>(const typename Mesh::point_type& )> sigma_fun,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<disk::static_vector<double, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        RealType dg_acoutic_l2_error = 0.0;
        RealType dg_elastic_l2_error = 0.0;
        RealType h = 10.0;

        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.face_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t e_cell_dof = n_ten_cbs + n_vec_cbs;

        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_cell_dof = n_vel_scal_cbs + n_scal_cbs;
        
        auto storage = msh.backend_storage();
        
        // Elastic dG unknowns 
        size_t e_cell_ind = 0;
        for (auto& e_chunk : assembler.get_e_material_data()) {

            auto& cell = storage->surfaces[e_chunk.first];

            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }

            { // tensor evaluation
            Matrix<RealType, Dynamic, 1> ten_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof, 0, n_ten_cbs, 1);
            auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.face_degree());
            Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, ten_basis, hho_di.face_degree());
            Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, ten_basis, sigma_fun);
            Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
            Matrix<RealType, Dynamic, 1> diff = real_dofs - ten_x_cell_dof;
            dg_elastic_l2_error += diff.dot(mass*diff);
            }
            
            e_cell_ind++;
        }
        
        size_t n_elastic_cell_dof = assembler.get_e_material_data().size() * e_cell_dof;

        // Acoustic dG unknowns 
        size_t a_cell_ind = 0;
        for (auto& a_chunk : assembler.get_a_material_data()) {
            
            auto& cell = storage->surfaces[a_chunk.first];
            
            RealType h_l = diameter(msh, cell);
            if (h_l < h) {
                h = h_l;
            }
            
            
            { // Vectorial evaluation
            Matrix<RealType, Dynamic, 1> flux_cell_dof = x_dof.block(a_cell_ind*a_cell_dof+n_elastic_cell_dof, 0, n_vel_scal_cbs, 1);
            Matrix<RealType, Dynamic, 1> real_dofs = flux_cell_dof;
      
            auto recdeg = hho_di.reconstruction_degree();
            auto rec_basis = make_scalar_monomial_basis(msh, cell, recdeg);
            auto rbs = disk::scalar_basis_size(recdeg, Mesh::dimension);
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_q_full  = make_stiffness_matrix(msh, cell, rec_basis);
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_q = Matrix<RealType, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
            mass_matrix_q = mass_matrix_q_full.block(1, 1, rbs-1, rbs-1);

            Matrix<RealType, Dynamic, 1> rhs = Matrix<RealType, Dynamic, 1>::Zero(rbs-1);
            Matrix<RealType, 1, 2> f_vec = Matrix<RealType, Dynamic, Dynamic>::Zero(1, 2);
            const auto qps = integrate(msh, cell, 2*recdeg);
            for (auto& qp : qps) {
                auto dphi = rec_basis.eval_gradients(qp.point());
                f_vec(0,0) = flux_fun(qp.point())[0];
                f_vec(0,1) = flux_fun(qp.point())[1];
                for (size_t i = 0; i < rbs-1; i++){
                    Matrix<RealType, 2, 1> phi_i = dphi.block(i+1, 0, 1, 2).transpose();
                    rhs(i,0) = rhs(i,0) + (qp.weight() * f_vec*phi_i)(0,0);
                }
            }
            real_dofs = mass_matrix_q.llt().solve(rhs);

            // std::cout << "flux_cell_dof = " << flux_cell_dof.size() << std::endl;
            // auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
            // std::cout << "cell_basis = " << cell_basis.size() << std::endl;
            // Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cell, cell_basis, hho_di.face_degree());
            // std::cout << "mass = " << mass_matrix_q.size() << std::endl;
            // Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cell, cell_basis, flux_fun);
            // std::cout << "rhs = " << rhs.size() << std::endl;
            // Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
            // std::cout << "real_dofs = " << real_dofs.size() << std::endl;
            Matrix<RealType, Dynamic, 1> diff = real_dofs - flux_cell_dof;
            // std::cout << "diff = " << diff.size() << std::endl;
            dg_acoutic_l2_error += diff.dot(mass_matrix_q*diff);
            // std::cout << std::endl;

        }

            a_cell_ind++;
        }
        
        tc.toc();
        error_file << "L2 errors of dG unknowns: " << std::endl;
        error_file << "Acoustic region : " << std::setprecision(16) << std::sqrt(dg_acoutic_l2_error) << std::endl;
        error_file << "Elastic  region : " << std::setprecision(16) << std::sqrt(dg_elastic_l2_error) << std::endl;
        error_file << std::endl;
        error_file.flush();
    }

    /// Compute energy 
    static double compute_elasto_acoustic_energy_four_field(Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & x_dof, std::ostream & energy_file = std::cout) {

        timecounter tc;
        tc.tic();
        
        using RealType = double;
        auto dim = Mesh::dimension;
        
        // Elastic basis functions
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), Mesh::dimension, Mesh::dimension);   
        size_t e_cell_dof = n_ten_cbs + n_vec_cbs;
        
        // Acoustic basis functions 
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension) - 1;  
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_cell_dof = n_vel_scal_cbs + n_scal_cbs;
        
        
        auto storage = msh.backend_storage();
        std::vector<RealType> elastic_energy_vec(msh.cells_size()); 
        std::vector<RealType> acoustic_energy_vec(msh.cells_size());    
        
        // Computation of the elastic energy 
        size_t e_cell_ind = 0;
        
        //Loop on elastic elements 
        for (auto& e_chunk : assembler.get_e_material_data()) {
            auto& cell = storage->surfaces[e_chunk.first];
            
            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.e_mass_operator(e_chunk.second, msh, cell);
            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_e_dof_data(e_chunk.first, msh, cell, x_dof);
            
            // Evaluation of the L2 norm of the velocity (vectorial component)
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs);  
            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof_loc.block(n_ten_cbs, 0, n_vec_cbs, 1);        
            Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * vec_cell_dof;
            Matrix<RealType, 1, 1> term_1 = vec_cell_dof.transpose() * v_mass_tested;
            elastic_energy_vec[e_cell_ind] = term_1(0,0);
            
            // Evaluation of the L2 norm of the strain tensor 
            Matrix<RealType, Dynamic, Dynamic> mass_sigma = mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs);
            Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs, 1);
            Matrix<RealType, Dynamic, 1> epsilon_mass = mass_sigma * sigma_dof;
            Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
            elastic_energy_vec[e_cell_ind] += term_2(0,0);
    
            e_cell_ind++;
            
        }
        
        // // Computation of the acoustic energy
        size_t a_cell_ind = 0;
        
        // Loop on acoustic elements 
        for (auto& a_chunk : assembler.get_a_material_data()) {
            auto& cell = storage->surfaces[a_chunk.first];
            
            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.a_mass_operator(a_chunk.second, msh, cell);
            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_a_dof_data(a_chunk.first, msh, cell, x_dof);
            
            // Evaluation of the L2 norm of the pressure (scalar component) 
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_p = mass_matrix.block(n_vel_scal_cbs, n_vel_scal_cbs, n_scal_cbs, n_scal_cbs); 
            Matrix<RealType, Dynamic, 1> p_cell_dof = x_dof_loc.block(n_vel_scal_cbs, 0, n_scal_cbs, 1);
            Matrix<RealType, Dynamic, 1> p_cell_mass_tested = mass_matrix_p * p_cell_dof;
            Matrix<RealType, 1, 1> term_3 = p_cell_dof.transpose() * p_cell_mass_tested;
            acoustic_energy_vec[a_cell_ind] = term_3(0,0);
            
            // Evaluation of the L2 norm of the acoustic velocity (vectorial component)
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(0, 0, n_vel_scal_cbs, n_vel_scal_cbs); 
            Matrix<RealType, Dynamic, 1> v_cell_dof = x_dof_loc.block(0, 0, n_vel_scal_cbs, 1);
            Matrix<RealType, Dynamic, 1> v_cell_mass_tested = mass_matrix_v * v_cell_dof;
            Matrix<RealType, 1, 1> term_4 = v_cell_dof.transpose() * v_cell_mass_tested;
            acoustic_energy_vec[a_cell_ind] += term_4(0,0);
            
            a_cell_ind++;
        }
        
        RealType elastic_energy_h = std::accumulate(elastic_energy_vec.begin(), elastic_energy_vec.end(), 0.0);
        RealType acoustic_energy_h = std::accumulate(acoustic_energy_vec.begin(), acoustic_energy_vec.end(), 0.0);
        elastic_energy_h  *= 0.5;
        acoustic_energy_h *= 0.5;
        RealType total_energy = elastic_energy_h + acoustic_energy_h;
        tc.toc();
        
        // std::cout << bold << cyan << "      Energy completed: " << tc << " seconds" << reset;
        energy_file << time << "         " << std::setprecision(16) 
        << elastic_energy_h  << "         " 
        << acoustic_energy_h << "         " << total_energy << std::endl;
        
        return total_energy;
    }
    
    static double compute_elasto_acoustic_energy_four_field_bassin(polygon_2d_mesh_reader<double> mesh_builder, Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, double & time, Matrix<double, Dynamic, 1> & x_dof, std::ostream & energy_file = std::cout) {

        timecounter tc;
        tc.tic();
        
        using RealType = double;
        auto dim = Mesh::dimension;
        
        // Elastic basis functions
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(),Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), Mesh::dimension, Mesh::dimension);   
        size_t e_cell_dof = n_ten_cbs + n_vec_cbs;
        
        // Acoustic basis functions 
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension) - 1;  
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_cell_dof = n_vel_scal_cbs + n_scal_cbs;
    
        auto storage = msh.backend_storage();
        std::vector<RealType> elastic_energy_vec(msh.cells_size()); 
        std::vector<RealType> socle_energy_vec(msh.cells_size()); 
        std::vector<RealType> sediment_energy_vec(msh.cells_size()); 
        std::vector<RealType> acoustic_energy_vec(msh.cells_size());    
        
        // Computation of the elastic energy 
        size_t e_cell_ind = 0;
        
        //Loop on elastic elements 
        for (auto& e_chunk : assembler.get_e_material_data()) {
            auto& cell = storage->surfaces[e_chunk.first];
            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.e_mass_operator(e_chunk.second, msh, cell);
            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_e_dof_data(e_chunk.first, msh, cell, x_dof);
            
            // Evaluation of the L2 norm of the velocity (vectorial component)
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(n_ten_cbs, n_ten_cbs, n_vec_cbs, n_vec_cbs);  
            Matrix<RealType, Dynamic, 1> vec_cell_dof = x_dof_loc.block(n_ten_cbs, 0, n_vec_cbs, 1);        
            Matrix<RealType, Dynamic, 1> v_mass_tested = mass_matrix_v * vec_cell_dof;
            Matrix<RealType, 1, 1> term_1 = vec_cell_dof.transpose() * v_mass_tested;
            elastic_energy_vec[e_cell_ind] = term_1(0,0);
            
            // Evaluation of the L2 norm of the strain tensor 
            Matrix<RealType, Dynamic, Dynamic> mass_sigma = mass_matrix.block(0, 0, n_ten_cbs, n_ten_cbs);
            Matrix<RealType, Dynamic, 1> sigma_dof = x_dof_loc.block(0, 0, n_ten_cbs, 1);
            Matrix<RealType, Dynamic, 1> epsilon_mass = mass_sigma * sigma_dof;
            Matrix<RealType, 1, 1> term_2 = sigma_dof.transpose() * epsilon_mass;
            if (mesh_builder.polygons[e_chunk.first].m_material == 1) {
                socle_energy_vec[e_cell_ind] = term_1(0,0) + term_2(0,0);
            }
            if (mesh_builder.polygons[e_chunk.first].m_material == 2) {
                sediment_energy_vec[e_cell_ind] = term_1(0,0) + term_2(0,0);
            }
            elastic_energy_vec[e_cell_ind] += term_2(0,0);
            
            e_cell_ind++;
            
        }
        
        // // Computation of the acoustic energy
        size_t a_cell_ind = 0;
        
        // Loop on acoustic elements 
        for (auto& a_chunk : assembler.get_a_material_data()) {
            auto& cell = storage->surfaces[a_chunk.first];
            
            Matrix<RealType, Dynamic, Dynamic> mass_matrix = assembler.a_mass_operator(a_chunk.second, msh, cell);
            Matrix<RealType, Dynamic, 1> x_dof_loc = assembler.gather_a_dof_data(a_chunk.first, msh, cell, x_dof);
            
            // Evaluation of the L2 norm of the pressure (scalar component) 
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_p = mass_matrix.block(n_vel_scal_cbs, n_vel_scal_cbs, n_scal_cbs, n_scal_cbs); 
            Matrix<RealType, Dynamic, 1> p_cell_dof = x_dof_loc.block(n_vel_scal_cbs, 0, n_scal_cbs, 1);
            Matrix<RealType, Dynamic, 1> p_cell_mass_tested = mass_matrix_p * p_cell_dof;
            Matrix<RealType, 1, 1> term_3 = p_cell_dof.transpose() * p_cell_mass_tested;
            acoustic_energy_vec[a_cell_ind] = term_3(0,0);
            
            // Evaluation of the L2 norm of the acoustic velocity (vectorial component)
            Matrix<RealType, Dynamic, Dynamic> mass_matrix_v = mass_matrix.block(0, 0, n_vel_scal_cbs, n_vel_scal_cbs); 
            Matrix<RealType, Dynamic, 1> v_cell_dof = x_dof_loc.block(0, 0, n_vel_scal_cbs, 1);
            Matrix<RealType, Dynamic, 1> v_cell_mass_tested = mass_matrix_v * v_cell_dof;
            Matrix<RealType, 1, 1> term_4 = v_cell_dof.transpose() * v_cell_mass_tested;
            acoustic_energy_vec[a_cell_ind] += term_4(0,0);
            
            a_cell_ind++;
        }
        
        RealType socle_energy_h    = std::accumulate(socle_energy_vec.begin(), socle_energy_vec.end(), 0.0);
        RealType sediment_energy_h = std::accumulate(sediment_energy_vec.begin(), sediment_energy_vec.end(), 0.0);
        RealType elastic_energy_h  = std::accumulate(elastic_energy_vec.begin(), elastic_energy_vec.end(), 0.0);
        RealType acoustic_energy_h = std::accumulate(acoustic_energy_vec.begin(), acoustic_energy_vec.end(), 0.0);
        socle_energy_h    *= 0.5;
        sediment_energy_h *= 0.5;
        elastic_energy_h  *= 0.5;
        acoustic_energy_h *= 0.5;
        RealType total_energy = elastic_energy_h + acoustic_energy_h;
        tc.toc();
        
        // std::cout << bold << cyan << "      Energy completed: " << tc << " seconds" << reset;
        energy_file << time << "         " << std::setprecision(16) 
        << socle_energy_h  << "         " 
        << sediment_energy_h  << "         " 
        << elastic_energy_h  << "         " 
        << acoustic_energy_h << "         " << total_energy << std::endl;
        
        return total_energy;
        
    }
  
    // Write silo file 
    static void write_silo_two_fields_elastoacoustic(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::map<size_t,elastic_material_data<double>> & e_material, std::map<size_t,acoustic_material_data<double>> & a_material, bool cell_centered_Q){

        timecounter tc;
        tc.tic();
        
        auto dim = Mesh::dimension;
        auto num_cells = msh.cells_size();
        auto num_points = msh.points_size();
        using RealType = double;
        std::vector<RealType> approx_ux, approx_uy;
        std::vector<RealType> approx_u;
        size_t e_cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
        size_t a_cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
        auto storage = msh.backend_storage();
        
        if (cell_centered_Q) {
            approx_ux.resize( num_cells );
            approx_uy.resize( num_cells );
            approx_u.resize( num_cells );
            
            size_t e_cell_ind = 0;
            for (auto& e_chunk : e_material) {
                
                auto& cell = storage->surfaces[e_chunk.first];
                auto bar = barycenter(msh, cell);
                approx_u.at(e_chunk.first) =( 0.0/ 0.0);
                
                // vector evaluation
                {
                    auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                    Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof, 0, e_cell_dof, 1);
                    auto t_phi = cell_basis.eval_functions( bar );
                    assert(t_phi.rows() == cell_basis.size());
                    auto uh = disk::eval(vec_x_cell_dof, t_phi);
                    approx_ux.at(e_chunk.first) = (uh(0,0));
                    approx_uy.at(e_chunk.first) = (uh(1,0));
                }
                e_cell_ind++;
            }
            
            size_t n_elastic_cell_dof = e_material.size() * e_cell_dof;
            size_t a_cell_ind = 0;
            for (auto& a_chunk : a_material) {
                
                auto& cell = storage->surfaces[a_chunk.first];
                auto bar = barycenter(msh, cell);
                approx_ux.at(a_chunk.first) = ( 0.0/0.0 );
                approx_uy.at(a_chunk.first) = ( 0.0/0.0 );
                
                // scalar evaluation
                {
                    auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                    Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof + n_elastic_cell_dof, 0, a_cell_dof, 1);
                    auto t_phi = cell_basis.eval_functions( bar );
                    RealType uh = scalar_cell_dof.dot( t_phi );
                    approx_u.at(a_chunk.first) = (uh);
                }
                a_cell_ind++;
            }


        }else{
            
            // Filling with nan (It is weird but useful in Paraview)
            approx_ux.resize( num_points , 0.0/ 0.0);
            approx_uy.resize( num_points , 0.0/ 0.0);
            approx_u.resize( num_points , 0.0/ 0.0);
            
            std::map<size_t,size_t> e_cell_index;
            std::map<size_t,size_t> a_cell_index;
            
            // elastic data
            size_t e_cell_ind = 0;
            for (auto chunk : e_material) {
                e_cell_index.insert(std::make_pair(chunk.first,e_cell_ind));
                e_cell_ind++;
            }

            // acoustic data
            size_t a_cell_ind = 0;
            for (auto chunk : a_material) {
                a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
                a_cell_ind++;
            }
            
            for (auto& e_chunk : e_material) {
                
                auto& cell = storage->surfaces[e_chunk.first];
                size_t e_cell_ind = e_cell_index[e_chunk.first];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                 Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof, 0, e_cell_dof, 1);
                
                auto points = cell.point_ids();
                size_t n_p = points.size();
                for (size_t l = 0; l < n_p; l++)
                {
                    auto pt_id = points[l];
                    auto pt_coord = *std::next(msh.points_begin(), pt_id);
                    // vector evaluation
                    {
                        auto t_phi = cell_basis.eval_functions( pt_coord );
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        approx_ux.at(pt_id) = uh(0,0);
                        approx_uy.at(pt_id) = uh(1,0);
                    }
                }
            }
            
            size_t n_elastic_cell_dof = e_material.size() * e_cell_dof;
            for (auto& a_chunk : a_material) {
                
                auto& cell = storage->surfaces[a_chunk.first];
                size_t a_cell_ind = a_cell_index[a_chunk.first];
                auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof + n_elastic_cell_dof, 0, a_cell_dof, 1);
                
                auto points = cell.point_ids();
                size_t n_p = points.size();
                for (size_t l = 0; l < n_p; l++)
                {
                    auto pt_id = points[l];
                    auto pt_coord = *std::next(msh.points_begin(), pt_id);
                    // scalar evaluation
                    {
                        auto t_phi = cell_basis.eval_functions( pt_coord );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.at(pt_id) = uh;
                    }
                }
            }
            
        }

        disk::silo_database silo;
        silo_file_name += std::to_string(it) + ".silo";
        silo.create(silo_file_name.c_str());
        silo.add_mesh(msh, "mesh");
        if (cell_centered_Q) {
            disk::silo_zonal_variable<double> vhx_silo("vhx", approx_ux);
            disk::silo_zonal_variable<double> vhy_silo("vhy", approx_uy);
            disk::silo_zonal_variable<double> vh_silo("vh", approx_u);
            silo.add_variable("mesh", vhx_silo);
            silo.add_variable("mesh", vhy_silo);
            silo.add_variable("mesh", vh_silo);
        }else{
            disk::silo_nodal_variable<double> vhx_silo("vhx", approx_ux);
            disk::silo_nodal_variable<double> vhy_silo("vhy", approx_uy);
            disk::silo_nodal_variable<double> vh_silo("vh", approx_u);
            silo.add_variable("mesh", vhx_silo);
            silo.add_variable("mesh", vhy_silo);
            silo.add_variable("mesh", vh_silo);
        }

        silo.close();
        tc.toc();
        std::cout << std::endl;
        std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
    }
    
    static void write_silo_four_fields_elastoacoustic(std::string silo_file_name, size_t it, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & x_dof, std::map<size_t,elastic_material_data<double>> & e_material, std::map<size_t,acoustic_material_data<double>> & a_material, bool cell_centered_Q) {

        timecounter tc;
        tc.tic();
    
        auto dim = Mesh::dimension;
        auto num_cells  = msh.cells_size(); 
        auto num_points = msh.points_size();
        
        using RealType  = double;
        std::vector<RealType> approx_ux, approx_uy, approx_uy_bis;
        std::vector<RealType> approx_u;
        
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(), Mesh::dimension, Mesh::dimension);
        size_t e_cell_dof = n_ten_cbs + n_vec_cbs;
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension) -1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_cell_dof = n_vel_scal_cbs + n_scal_cbs;
        
        auto storage = msh.backend_storage();
        
        if (cell_centered_Q) {
            approx_ux.resize(num_cells);
            approx_uy.resize(num_cells);
            approx_uy_bis.resize(num_cells);
            approx_u.resize (num_cells);
            
            size_t e_cell_ind = 0;
            for (auto& e_chunk : e_material) {
	
                auto& cell = storage->surfaces[e_chunk.first];
                auto bar = barycenter(msh, cell);
                approx_u.at(e_chunk.first) = (0.0/0.0);
                
                // vector evaluation
                {
                    auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                    Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof 
                    + n_ten_cbs, 0, n_vec_cbs, 1);
                    auto t_phi = cell_basis.eval_functions( bar );
                    assert(t_phi.rows() == cell_basis.size());
                    auto uh = disk::eval(vec_x_cell_dof, t_phi);
                    approx_ux.at(e_chunk.first) = (uh(0,0));
                    approx_uy.at(e_chunk.first) = (uh(1,0));
                    approx_uy_bis.at(e_chunk.first) = (uh(1,0));
                }
                e_cell_ind++;
            }
            
            size_t n_elastic_cell_dof = e_material.size() * e_cell_dof;
            size_t a_cell_ind = 0;
            for (auto& a_chunk : a_material) {
	
                auto& cell = storage->surfaces[a_chunk.first];
                auto bar = barycenter(msh, cell);
                approx_ux.at(a_chunk.first) = (0.0/0.0);
                approx_uy.at(a_chunk.first) = (0.0/0.0);
                approx_uy_bis.at(a_chunk.first) = (0.0/0.0);
                
                // scalar evaluation
                {
                    auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                    Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof +
                    n_vel_scal_cbs +
                    n_elastic_cell_dof, 0,
                    n_scal_cbs, 1);
                    auto t_phi = cell_basis.eval_functions( bar );
                    RealType uh = scalar_cell_dof.dot( t_phi );
                    approx_u.at(a_chunk.first) = (uh);
                }
                a_cell_ind++;
            }
        }
        
        else {
      
            // Filling with nan (It is weird but useful in Paraview)
            approx_ux.resize(num_points,0.0/0.0);
            approx_uy.resize(num_points,0.0/0.0);
            approx_uy_bis.resize(num_points,0.0/0.0);
            approx_u.resize(num_points,0.0/0.0);
            
            std::map<size_t,size_t> e_cell_index;
            std::map<size_t,size_t> a_cell_index;
            size_t e_cell_ind = 0;
            for (auto chunk : e_material) {
                e_cell_index.insert(std::make_pair(chunk.first,e_cell_ind));
                e_cell_ind++;
            }
            
            // acoustic data
            size_t a_cell_ind = 0;
            for (auto chunk : a_material) {
                a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
                a_cell_ind++;
            }
            
            for (auto& e_chunk : e_material) {
                
                auto& cell = storage->surfaces[e_chunk.first];
                size_t e_cell_ind = e_cell_index[e_chunk.first];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> vec_x_cell_dof = x_dof.block(e_cell_ind*e_cell_dof
                + n_ten_cbs, 0, n_vec_cbs, 1);
                
                auto points = cell.point_ids();
                size_t n_p  = points.size();
                // problème sur cette boucle 
                for (size_t l = 0; l < n_p; l++) {
                    auto pt_id = points[l];
                    auto pt_coord = *std::next(msh.points_begin(), pt_id);
                    
                    // vector evaluation
                    {
                        auto t_phi = cell_basis.eval_functions(pt_coord);
                        assert(t_phi.rows() == cell_basis.size());
                        auto uh = disk::eval(vec_x_cell_dof, t_phi);
                        
                        approx_ux.at(pt_id) = uh(0,0);
                        approx_uy.at(pt_id) = uh(1,0);
                        approx_uy_bis.at(pt_id) = std::sqrt( std::abs(uh(0,0))*std::abs(uh(0,0))
                        + std::abs(uh(1,0))*std::abs(uh(1,0)) );
                    }
                }
            }
            
            size_t n_elastic_cell_dof = e_material.size() * e_cell_dof;
            for (auto& a_chunk : a_material) {	
                auto& cell = storage->surfaces[a_chunk.first];
                size_t a_cell_ind = a_cell_index[a_chunk.first];
                auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(a_cell_ind*a_cell_dof +
                n_vel_scal_cbs + n_elastic_cell_dof,
                0, n_scal_cbs, 1);
                
                auto points = cell.point_ids();
                size_t n_p = points.size();
                // std::cout << n_p <<std::cout;
                for (size_t l = 0; l < n_p; l++) {
                    auto pt_id = points[l];
                    auto pt_coord = *std::next(msh.points_begin(), pt_id);
                    // scalar evaluation
                    {
                        auto t_phi = cell_basis.eval_functions( pt_coord );
                        RealType uh = scalar_cell_dof.dot( t_phi );
                        approx_u.at(pt_id) = uh;
                    }
                }
            }
            
        }
    
        disk::silo_database silo;
        silo_file_name += std::to_string(it) + ".silo";
        silo.create(silo_file_name.c_str());
        silo.add_mesh(msh, "mesh");
        if (cell_centered_Q) {
            disk::silo_zonal_variable<double> vhx_silo("vhx", approx_ux);
            disk::silo_zonal_variable<double> vhy_silo("vhy", approx_uy);
            disk::silo_zonal_variable<double> vhy_bis_silo("norm_v", approx_uy_bis);
            disk::silo_zonal_variable<double> vh_silo("vh", approx_u);
            silo.add_variable("mesh", vhx_silo);
            silo.add_variable("mesh", vhy_silo);
            silo.add_variable("mesh", vhy_bis_silo);
            silo.add_variable("mesh", vh_silo);
        }
        else {
            disk::silo_nodal_variable<double> vhx_silo("vhx", approx_ux);
            disk::silo_nodal_variable<double> vhy_silo("vhy", approx_uy);
            disk::silo_nodal_variable<double> vhy_bis_silo("norm_v", approx_uy_bis);
            disk::silo_nodal_variable<double> vh_silo("vh", approx_u);
            silo.add_variable("mesh", vhx_silo);
            silo.add_variable("mesh", vhy_silo);
            silo.add_variable("mesh", vhy_bis_silo);
            silo.add_variable("mesh", vh_silo);
        }
        
        silo.close();
        tc.toc();
        // std::cout << "      ";
        // std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset
        //           << std::endl;
    
    }
  
    /// Record acoustic data 
    static void record_acoustic_data_elasto_acoustic_four_fields( size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index,  Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, bool e_side_Q, std::ostream & seismogram_file = std::cout) {
       
        using RealType = double;
        auto dim = Mesh::dimension;
        
        // Basis 
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_n_cbs = n_vel_scal_cbs + n_scal_cbs;
        size_t a_n_cell_dof = assembler.get_a_n_cells_dof();
        
        RealType vh = 0.0;
        
        typename Mesh::point_type pt = pt_cell_index.first;
        
        // FIND CELL
        if (pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_acoustic_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vh\"" << std::endl;
        }
    
        size_t cell_ind = pt_cell_index.second;
        auto& cell = msh.backend_storage()->surfaces[cell_ind];
        auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
        Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_a_dof_data(cell_ind, msh, cell, x_dof);
        Matrix<RealType, Dynamic, 1> scalar_cell_dof = all_dofs.block(n_vel_scal_cbs,0,n_scal_cbs,1);
        auto t_phi = cell_basis.eval_functions( pt );
        vh = scalar_cell_dof.dot( t_phi );
        
        seismogram_file << it << "," << std::setprecision(16) <<  vh << std::endl;
        seismogram_file.flush();
        
    }
    
    /// Record elastic data 
    static void record_elastic_data_elasto_acoustic_four_fields( size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index,  Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, bool e_side_Q, std::ostream & seismogram_file = std::cout) {

        timecounter tc;
        tc.tic();
        
        using RealType = double;
        auto dim = Mesh::dimension;
        
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t e_n_cbs = n_ten_cbs + n_vec_cbs;
        size_t e_n_cell_dof = assembler.get_e_n_cells_dof();
        
        Matrix<double, Dynamic, Dynamic> eh = Matrix<double, Dynamic, Dynamic>::Zero(2, 2);
        
        typename Mesh::point_type pt = pt_cell_index.first;    
        if (pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_elastic_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"Sxx\"" << "," << "\"Sxy\"" << "," << "\"Syy\"" << std::endl;
        }
        size_t cell_ind = pt_cell_index.second;
        auto& cell = msh.backend_storage()->surfaces[cell_ind];
        
        // stress tensor evaluation
        size_t ten_dofs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        auto ten_basis = make_sym_matrix_monomial_basis(msh, cell, hho_di.grad_degree());
        Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_e_dof_data(cell_ind, msh, cell, x_dof);
        Matrix<RealType, Dynamic, 1> ten_x_cell_dof = all_dofs.block(0, 0, ten_dofs, 1);
        auto t_phi = ten_basis.eval_functions(pt);
        eh = disk::eval(ten_x_cell_dof, t_phi);
        
        tc.toc();
        seismogram_file << it << "," << std::setprecision(16) <<  eh(0,0) << "," << std::setprecision(16) <<  eh(0,1) << "," << std::setprecision(16) << eh(1,1) << std::endl;
        seismogram_file.flush();
        
    }
    
    /// Record velocity 
    static void record_velocity_data_elasto_acoustic_two_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_two_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & u_dof, Matrix<double, Dynamic, 1> & v_dof, bool e_side_Q, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;

        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        if(pt_cell_index.second == -1){
            std::set<size_t> cell_indexes = find_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
        }

        {
            size_t cell_ind = pt_cell_index.second;
            
            
            if(e_side_Q){// vector evaluation
                
                size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_e_dof_data(cell_ind, msh, cell, v_dof);
                Matrix<RealType, Dynamic, 1> vec_x_cell_dof = all_dofs.block(0, 0, cell_dof, 1);
                auto t_phi = cell_basis.eval_functions( pt );
                assert(t_phi.rows() == cell_basis.size());
                vh = disk::eval(vec_x_cell_dof, t_phi);
            }
            else
            {// reconstructed velocity evaluation

                size_t cell_dof = disk::scalar_basis_size(hho_di.cell_degree(), dim);
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto rec_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
                auto gr = make_scalar_hho_laplacian(msh, cell, hho_di);
                Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_a_dof_data(cell_ind, msh, cell, u_dof);
                Matrix<RealType, Dynamic, 1> recdofs = gr.first * all_dofs;
                
                auto t_dphi = rec_basis.eval_gradients( pt );
                Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();

                for (size_t i = 1; i < t_dphi.rows(); i++){
                    grad_uh = grad_uh + recdofs(i-1)*t_dphi.block(i, 0, 1, 2);
                }
                acoustic_material_data<RealType> a_mat = assembler.get_a_material_data().find(cell_ind)->second;
                RealType rho = a_mat.rho();
                vh = (1.0/rho)*(grad_uh);
            }
            
        }
        tc.toc();
        std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }
    
    static void record_velocity_data_elasto_acoustic_four_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, bool e_side_Q, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t e_n_cbs = n_ten_cbs + n_vec_cbs;
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_n_cbs = n_vel_scal_cbs + n_scal_cbs;
        size_t e_n_cell_dof = assembler.get_e_n_cells_dof();
        size_t a_n_cell_dof = assembler.get_a_n_cells_dof();
        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;

        if (e_side_Q) {
            // FIND CELL
            if (pt_cell_index.second == -1){
                std::set<size_t> cell_indexes = find_elastic_cells(pt, msh, true);
                size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
                assert(cell_index != -1);
                pt_cell_index.second = cell_index;
                seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
            }
            // RECORD DATA
            size_t cell_ind = pt_cell_index.second;
            size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
            auto& cell = msh.backend_storage()->surfaces[cell_ind];
            auto cell_basis = make_vector_monomial_basis(msh, cell, hho_di.cell_degree());
            Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_e_dof_data(cell_ind, msh, cell, x_dof);
            Matrix<RealType, Dynamic, 1> vec_x_cell_dof = all_dofs.block(0 + n_ten_cbs, 0, cell_dof, 1);
            auto t_phi = cell_basis.eval_functions(pt);
            assert(t_phi.rows() == cell_basis.size());
            vh = disk::eval(vec_x_cell_dof, t_phi);
        }
        else {
            // FIND CELL
            if (pt_cell_index.second == -1){
                std::set<size_t> cell_indexes = find_acoustic_cells(pt, msh, true);
                size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
                assert(cell_index != -1);
                pt_cell_index.second = cell_index;
                seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
            }
            // RECORD DATA
            size_t cell_ind = pt_cell_index.second;
            size_t cell_dof = disk::vector_basis_size(hho_di.cell_degree(), dim, dim);
            auto& cell = msh.backend_storage()->surfaces[cell_ind];
            auto cell_basis = make_scalar_monomial_basis(msh, cell, hho_di.reconstruction_degree());
            Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_a_dof_data(cell_ind, msh, cell, x_dof);
            Matrix<RealType, Dynamic, 1> flux_cell_dof = all_dofs.block(0, 0, cell_dof, 1);
            auto t_dphi = cell_basis.eval_gradients( pt );
            Matrix<RealType, 1, 2> flux_uh = Matrix<RealType, 1, 2>::Zero();
            for (size_t i = 1; i < t_dphi.rows(); i++)
                flux_uh = flux_uh + flux_cell_dof(i-1)*t_dphi.block(i, 0, 1, 2);
            vh = flux_uh;
        }
        
        tc.toc();
        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();
    }
    
    /// Record interface data
    static void record_acoustic_pressure_coupling_data_elasto_acoustic_four_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, bool e_side_Q, std::set<size_t> interface_face_indexes, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        size_t n_vel_scal_cbs = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t n_scal_cbs = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t a_n_cbs = n_vel_scal_cbs + n_scal_cbs;
        size_t a_n_fbs = disk::scalar_basis_size(hho_di.face_degree(), Mesh::dimension - 1);
        RealType vh = 0.0;

        typename Mesh::point_type pt = pt_cell_index.first;
        typename Mesh::face_type fb;

        // FIND FACE
        if (pt_cell_index.second == -1) {
            seismogram_file << "\"Time\"" << "," << "\"vh\"" << std::endl;
            std::set<size_t> cell_indexes = find_acoustic_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
        }
        auto& cell = msh.backend_storage()->surfaces[pt_cell_index.second];
        auto cell_faces = faces(msh, cell);
        for (auto face : cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) 
                fb = face;
        }
        
        auto fc_id = msh.lookup(fb);
        auto face_basis = make_scalar_monomial_basis(msh, fb, hho_di.face_degree());
        Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_a_dof_data(pt_cell_index.second, msh, cell, x_dof);
        Matrix<RealType, Dynamic, 1> face_dofs = all_dofs.block(a_n_cbs, 0, a_n_fbs, 1);
        auto t_phi = face_basis.eval_functions(pt);
        vh = face_dofs.dot(t_phi);

        seismogram_file << it << "," << std::setprecision(16) <<  vh << std::endl;
        seismogram_file.flush();

    }

    static void record_elastic_velocity_coupling_data_elasto_acoustic_four_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, elastoacoustic_four_fields_assembler<Mesh> & assembler, Matrix<double, Dynamic, 1> & x_dof, bool e_side_Q, std::set<size_t> interface_face_indexes, std::ostream & seismogram_file = std::cout){

        timecounter tc;
        tc.tic();

        using RealType = double;
        auto dim = Mesh::dimension;
        
        size_t n_ten_cbs = disk::sym_matrix_basis_size(hho_di.grad_degree(), Mesh::dimension, Mesh::dimension);
        size_t n_vec_cbs = disk::vector_basis_size(hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
        size_t e_n_cbs = n_ten_cbs + n_vec_cbs;
        size_t e_n_cell_dof = assembler.get_e_n_cells_dof();
        size_t e_n_fbs = disk::vector_basis_size(hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
        Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

        typename Mesh::point_type pt = pt_cell_index.first;
        typename Mesh::face_type fb;

        // FIND FACE
        if (pt_cell_index.second == -1) {
            seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
            std::set<size_t> cell_indexes = find_elastic_cells(pt, msh, true);
            size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
            assert(cell_index != -1);
            pt_cell_index.second = cell_index;
        }
        auto& cell = msh.backend_storage()->surfaces[pt_cell_index.second];
        auto cell_faces = faces(msh, cell);
        for (auto face : cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) 
                fb = face;
        }
        
        // RECORD DATA
        auto fc_id = msh.lookup(fb);
        auto face_basis = make_vector_monomial_basis(msh, fb, hho_di.face_degree());
        Matrix<RealType, Dynamic, 1> all_dofs = assembler.gather_e_dof_data(pt_cell_index.second, msh, cell, x_dof);
        Matrix<RealType, Dynamic, 1> face_dofs = all_dofs.block(e_n_cbs, 0, e_n_fbs, 1);
        auto t_phi = face_basis.eval_functions(pt);
        vh = disk::eval(face_dofs, t_phi);

        seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
        seismogram_file.flush();

    }

};



#endif /* postprocessor_hpp */

