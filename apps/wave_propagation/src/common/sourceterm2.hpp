
#pragma once
#ifndef source_term_hpp
#define source_term_hpp

#include <iomanip>
#include "../common/acoustic_one_field_assembler.hpp"
#include "../common/acoustic_two_fields_assembler.hpp"
#include "../common/elastodynamic_one_field_assembler.hpp"
#include "../common/elastodynamic_two_fields_assembler.hpp"
#include "../common/elastodynamic_three_fields_assembler.hpp"
#include "../common/elastoacoustic_two_fields_assembler.hpp"
#include "../common/elastoacoustic_four_fields_assembler.hpp"

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class source_term {
    
public:
    
    using RealType = double;
    
    // find distances to the requested point    
    static std::vector<RealType> make_distances(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false) {

        auto norm =  [](const typename Mesh::point_type& a, const typename Mesh::point_type& b ) -> RealType {
            RealType dx = (b.x() - a.x());
            RealType dy = (b.y() - a.y());
            RealType norm = std::sqrt(dx*dx + dy*dy);
            return norm;
        };

        size_t np = msh.points_size();
        std::vector<RealType> distances(np);

        size_t ip = 0;
        for (auto& point : msh.backend_storage()->points)
        {
            RealType dist = norm(pt,point);
            distances[ip] = dist;
            ip++;
        }

        return distances;
    }

    // find minimum distances to the requested point 
    static RealType min_distance(std::vector<RealType> distances, bool verbose_Q = false) {
        RealType min_dist; 
        if(verbose_Q){
            min_dist = *std::min_element(distances.begin(), distances.end());
        }
        return min_dist;
    }

    /// Find the cells associated to the requested point
    static std::set<size_t> find_cells_source(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false){

        std::vector<RealType> distances = make_distances(pt, msh, true);
        size_t index = std::min_element(distances.begin(),distances.end()) - distances.begin();
        if(verbose_Q){
            RealType min_dist = *std::min_element(distances.begin(), distances.end());
            typename Mesh::point_type nearest_point = msh.backend_storage()->points.at(index);
            std::cout << "Nearest point detected : " << std::endl;
            std::cout << "  x =  " << nearest_point.x() << std::endl;
            std::cout << "  y =  " << nearest_point.y() << std::endl;
            std::cout << "Distance = " << min_dist << std::endl;
            std::cout << "Global index = " << index << std::endl;
        }

        std::set<size_t> cell_indexes;

        size_t cell_i = 0;
        for (auto& cell : msh)
        {
            auto points = cell.point_ids();
            size_t n_p = points.size();
            for (size_t l = 0; l < n_p; l++)
            {
                auto pt_id = points[l];
                if(index == pt_id){
                    cell_indexes.insert(cell_i);
                }
            }
            cell_i++;
        }
        
        if(verbose_Q){
            std::cout << "Detected cells indexes : " << std::endl;
            for(auto index : cell_indexes){
                std::cout << index << std::endl;
            }
        }

        return cell_indexes;
    
    }
    
    /// Pick the cell that contains the requested point
    static std::vector<size_t> pick_cell_source(typename Mesh::point_type & pt, Mesh & msh, std::set<size_t> & cell_indexes, bool verbose_Q = false){
        
        std::vector<RealType> distances = make_distances(pt, msh, true);

        if (min_distance(distances, true) == 0.0) {
        std::vector<size_t> cells_sources(4);
            int i = 0;
            for (auto cells : cell_indexes) {
                cells_sources[i] = cells;
                i++;
            }
            std::cout << "Pas de problème quand on est sur un point du maillage" << std::endl;
            std::cout << "Quand on est au milieu d'une maille à implémenter" << std::endl;
            return cells_sources;
        }
  
        else {
        std::vector<size_t> cells_sources(1);
            auto triangle_member_Q = [] (typename Mesh::point_type & p, typename Mesh::point_type & p0, typename Mesh::point_type & p1, typename Mesh::point_type & p2)
            {
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
        std::cout << "Problème laaaaaaa" << std::endl;     
            size_t n_cells = cell_indexes.size();
            if (n_cells == 1) {
                size_t first_index = *cell_indexes.begin();
                cells_sources[0] = first_index;
                return cells_sources;
            }
            bool is_member_Q = false;
            for(auto index : cell_indexes){
                auto& cell = msh.backend_storage()->surfaces[index];
                auto bar = barycenter(msh, cell);
                auto points = cell.point_ids();
                size_t n_p = points.size();
                
                // building teselation
                std::vector<std::vector<typename Mesh::point_type>> triangles(n_p);
                for (size_t l = 0; l < n_p; l++)
                {

                    std::vector<typename Mesh::point_type> chunk(3);
                    if( l == n_p - 1){
                        chunk[0] = msh.backend_storage()->points.at(points[l]);
                        chunk[1] = msh.backend_storage()->points.at(points[0]);
                        chunk[2] = bar;
                    }else{
                        chunk[0] = msh.backend_storage()->points.at(points[l]);
                        chunk[1] = msh.backend_storage()->points.at(points[l+1]);
                        chunk[2] = bar;
                    }
                    triangles[l] = chunk;
                }
                
                // check whether the point is memeber of any triangle
                int j = 0;
                for (auto triangle : triangles) {
                    is_member_Q = triangle_member_Q(pt,triangle[0],triangle[1],triangle[2]);
                    if (is_member_Q) {
                        std::cout << "Detected cell index = " << index << std::endl;
                        cells_sources[j] = index;
                        return cells_sources;
                    }
                    j++;
                }

            }
            
            if(!is_member_Q){
                if(verbose_Q){
                    std::cout << "Point is not member of cells set. Returning cell_indexes[0] " << std::endl;
                }
                size_t first_index = *cell_indexes.begin();
                cells_sources[0] = first_index;
                return cells_sources;
            }
            
        return cells_sources;

        }
    }
    
    // Punctual Source
    static Matrix<RealType, Dynamic, 1> punctual_source(size_t it, double dt, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, disk::hho_degree_info & hho_di, Matrix<RealType, Dynamic, 1> & x_dof){

        auto dim = Mesh::dimension;
        size_t n_scal_dof = disk::scalar_basis_size(hho_di.cell_degree(), Mesh::dimension);
        size_t n_vec_dof = disk::scalar_basis_size(hho_di.reconstruction_degree(), Mesh::dimension)-1;
        size_t cell_dof = n_scal_dof + n_vec_dof;

        typename Mesh::point_type pt = pt_cell_index.first;

        std::vector<size_t> cells_index;
        std::set<size_t> cell_indexes;

        if(pt_cell_index.second == -1){
            cell_indexes = find_cells_source(pt, msh, true);
            cells_index = pick_cell_source(pt, msh, cell_indexes, true);
            std::cout << cells_index.size() << std::endl;
            assert(cells_index[0] != -1);
            pt_cell_index.second = cells_index[0];
        }

        if (cells_index.size() == 4) {
            std::cout << "source sur un point du maillage" << std::endl;
        }

        else {
            size_t cell_ind = pt_cell_index.second;
            // scalar evaluation
            {               

            std::cout << "source sur un point du maillage" << std::endl;
                Matrix<RealType, Dynamic, 1> scalar_cell_dof = x_dof.block(cell_ind*cell_dof+n_vec_dof, 0, n_scal_dof, 1);
                auto& cell = msh.backend_storage()->surfaces[cell_ind];
                auto cell_basis = disk::make_scalar_monomial_basis(msh, cell, hho_di.cell_degree());
                auto t_phi = cell_basis.eval_functions( pt );
                std::cout << cell_ind << std::endl;

                for (int i = scalar_cell_dof.size()-1; i < scalar_cell_dof.size(); i++) {
                    scalar_cell_dof(i) = scalar_cell_dof(i) - t_phi(i)*ricker_fluid(it,dt);
                }
                                
                x_dof.block(cell_ind*cell_dof+n_vec_dof, 0, n_scal_dof, 1) = scalar_cell_dof;
            }
        } 

        return x_dof;

    }

    static RealType ricker_fluid(size_t it, RealType dt){

        RealType tau, freq, Ricker, Ricker_fl, alpha, sigma, sigmabis, sigma2, sigma3, tn;

        tn    = it * dt;
        tau   = 0.4;
        freq  = 5.0; 
        alpha = - M_PI*M_PI * freq*freq;
        
        if ( tn < 2.5*tau ) {
            sigma  = alpha*(tn-tau)*(tn-tau);
            Ricker = 2.0 * alpha * (1+2*sigma) * exp(sigma);
            sigmabis  = M_PI*freq*(tn-tau);
            sigma2 = sigmabis*sigmabis;
            sigma3 = M_PI*freq;
            Ricker_fl = -4.0*sigmabis*sigma3*exp(-sigma2)-2.0*sigmabis*sigma3*Ricker;
        }
        else Ricker_fl = 0.0; 

        std::cout << Ricker_fl << std::endl;

        return Ricker_fl;
    }

};


#endif /* postprocessor_hpp */
