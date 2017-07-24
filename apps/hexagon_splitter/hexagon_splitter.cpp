/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>
#include <fstream>

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "common/eigen.hpp"


std::ostream&
operator<<(std::ostream& os, const std::vector<size_t>& vec)
{
    for (auto& e : vec)
        os << e+1 << " ";
    
    return os;
}

int main(int argc, char **argv)
{
    
    if (argc < 2)
    {
        std::cout << "Please specify mesh" << std::endl;
        return 1;
    }
    
    char *mesh_fn = argv[1];
    
    
    using RealType = double;
    
    typedef disk::generic_mesh<RealType,2> mesh_type;

    mesh_type msh = disk::load_fvca5_2d_mesh<RealType>(mesh_fn);
    
    dump_to_matlab(msh, "pre_split.m");
    
    auto old_storage = msh.backend_storage();
    
    size_t num_hexagonal_cells = 0;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        
        if ( fcs.size() == 6 )
            num_hexagonal_cells++;
    }
    
    std::cout << "Mesh has " << num_hexagonal_cells << " hex cells" << std::endl;
    
    std::vector<std::vector<size_t>> mfaces;
    for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
    {
        auto fc = *itor;
        
        auto ptids = fc.point_ids();
        
        mfaces.push_back( disk::convert_to<size_t>(ptids) );
    }
    
    std::vector<std::vector<size_t>> cells;
    size_t splitted_hexa_cells = 0;
    for (auto& cl : msh)
    {
        auto ptids = cl.point_ids();
        
        if (splitted_hexa_cells > num_hexagonal_cells/2 )
        {
            cells.push_back( disk::convert_to<size_t>(ptids) );
            continue;
        }
        
        auto fcs = faces(msh, cl);
        if ( fcs.size() != 6)
        {
            cells.push_back( disk::convert_to<size_t>(ptids) );
            continue;
        }
        
        mfaces.push_back({ptids[0], ptids[2]});
        mfaces.push_back({ptids[0], ptids[3]});
        mfaces.push_back({ptids[0], ptids[4]});
        cells.push_back({ptids[0], ptids[1], ptids[2]});
        cells.push_back({ptids[0], ptids[2], ptids[3]});
        cells.push_back({ptids[0], ptids[3], ptids[4]});
        cells.push_back({ptids[0], ptids[4], ptids[5]});
        splitted_hexa_cells++;
    }
    
    struct comp {
    
        bool operator()(const std::vector<size_t>& a, const std::vector<size_t>& b)
        {
            return (a.size() < b.size());
        }
        
    };
    
    std::sort(mfaces.begin(), mfaces.end());
    std::sort(cells.begin(), cells.end(), comp{});
    
    std::ofstream os("splitted.typ1");
    
    os << "vertices\n" << old_storage->points.size() << std::endl;
    for (auto& pt : old_storage->points)
        os << pt.x() << " " << pt.y() << std::endl;
    
    
    auto itor = cells.begin();
    
    os << "triangles" << std::endl;
    os << std::count_if(cells.begin(), cells.end(), [](const std::vector<size_t>& v) { return v.size() == 3; }) << std::endl;
    while ( (*itor).size() == 3 )
        os << (*itor++) << std::endl;
    
    os << "quadrangles" << std::endl;
    os << std::count_if(cells.begin(), cells.end(), [](const std::vector<size_t>& v) { return v.size() == 4; }) << std::endl;
    while ( (*itor).size() == 4 )
        os << (*itor++) << std::endl;
    
    os << "pentagons" << std::endl;
    os << std::count_if(cells.begin(), cells.end(), [](const std::vector<size_t>& v) { return v.size() == 5; }) << std::endl;
    while ( (*itor).size() == 5 )
        os << (*itor++) << std::endl;
    
    os << "hexagons" << std::endl;
    os << std::count_if(cells.begin(), cells.end(), [](const std::vector<size_t>& v) { return v.size() == 6; }) << std::endl;
    while ( (*itor).size() == 6 )
        os << (*itor++) << std::endl;
    
    os << "edges of the boundary" << std::endl;
    os << msh.boundary_faces_size() << std::endl;
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto fc = *itor;
        auto ptids = fc.point_ids();
        os << ptids[0]+1 << " " << ptids[1]+1 << std::endl;
    }
    
    
    std::vector<std::set<size_t>> face_to_cell;
    face_to_cell.resize( mfaces.size() );
    
    size_t cell_num = 1;
    for (auto& cl : cells)
    {
        for (size_t i = 0; i < cl.size(); i++)
        {
            auto sz = cl.size();
            auto p1 = cl.at(i);
            auto p2 = cl.at((i+1)%sz);
            
            if (p1 > p2)
                std::swap(p1, p2);
            
            std::vector<size_t> e({p1, p2});
            
            auto e_it = std::lower_bound(mfaces.begin(), mfaces.end(), e);
            auto e_num = std::distance(mfaces.begin(), e_it);
            face_to_cell.at(e_num).insert(cell_num);
        }
        cell_num++;
    }
    
    
    
    os << "all edges" << std::endl;
    os << mfaces.size() << std::endl;
    size_t i = 0;
    for (auto itor = mfaces.begin(); itor != mfaces.end(); itor++)
    {
        auto fc = *itor;
        os << fc << " ";
        for (auto& f : face_to_cell.at(i))
            os << f << " ";
        if (face_to_cell.at(i).size() == 1)
            os << "0";
        os << std::endl;
        i++;
    }
    
    
    os.close();
    return 0;

}




