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
#include <vector>
#include <array>
#include <map>
#include <list>
#include <random>

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"

namespace disk {

template<typename T>
class mesh_hierarchy
{
    typedef disk::simplicial_mesh<T,2>      mesh_type;
    typedef typename mesh_type::face        mesh_edge_type;
    typedef typename mesh_type::cell        mesh_triangle_type;
    typedef typename mesh_type::point_type  point_type;
    
    std::vector<mesh_type>                  meshes;
    std::vector< std::vector<std::array<size_t, 4>> >   coarse_to_fine;
    
    struct edge
    {
        edge(){}
        
        edge(const mesh_type& msh, const mesh_edge_type& e)
        {
            auto ptids  = e.point_ids();
            p0          = ptids[0];
            p1          = ptids[1];
            assert(p0 < p1);
            is_boundary = msh.is_boundary(e);
            is_broken   = false;
        }
        
        edge(size_t ap0, size_t ap1, bool bnd = false)
        {
            assert (ap0 != ap1);
            p0                  = ap0;
            p1                  = ap1;
            if (p0 > p1)
                std::swap(p0, p1);

            is_boundary         = bnd;
            is_broken           = false;
        }
        
        edge(size_t ap0, size_t ap1, size_t bid, bool bnd)
        {
            assert (ap0 != ap1);
            p0                  = ap0;
            p1                  = ap1;
            if (p0 > p1)
                std::swap(p0, p1);
            
            is_boundary         = bnd;
            is_broken           = false;
            boundary_id         = bid;
        }
        
        friend bool operator<(const edge& a, const edge& b) {
            assert(a.p0 < a.p1);
            assert(b.p0 < b.p1);
            return ( a.p0 < b.p0 or (a.p0 == b.p0 and a.p1 < b.p1) );
        }
        
        friend bool operator==(const edge& a, const edge& b) {
            return ( a.p0 == b.p0 and a.p1 == b.p1 );
        }
        
        friend std::ostream& operator<<(std::ostream& os, const edge& e) {
            os << "Edge: " << e.p0 << " " << e.p1;
            if (e.is_broken) os << ", broken";
            if (e.is_boundary) os << ", boundary";
            return os;
        }
        
        size_t  p0, p1, pb;
        bool    is_boundary, is_broken;
        size_t  boundary_id;
    };
    
    
    struct triangle
    {
        triangle() {}
        triangle(const mesh_type& msh, const mesh_triangle_type& tr)
        {
            auto ptids = tr.point_ids();
            std::copy(ptids.begin(), ptids.end(), p.begin());
            std::sort(p.begin(), p.end());
        }
        
        triangle(size_t ap0, size_t ap1, size_t ap2) : p{ap0, ap1, ap2}
        {
            std::sort(p.begin(), p.end());
        }
        
        std::array<size_t, 3> p;
        
        friend bool operator<(const triangle& a, const triangle& b) {
            return std::lexicographical_compare( a.p.begin(), a.p.end(),
                                                 b.p.begin(), b.p.end() );
        }
        
        friend std::ostream& operator<<(std::ostream& os, const triangle& t) {
            os << "Triangle: " << t.p[0] << " " << t.p[1] << " " << t.p[2];
            return os;
        }
    };
    
    /* Internal representation stuff */
    std::vector<point_type>                     m_points;
    std::vector<edge>                           m_edges;
    std::list<triangle>                         m_triangles;
    std::map<triangle, std::array<triangle,4>>  coarse_to_fine_int;
    
    void refine()
    {
        /* break all the edges */
        for (auto& e : m_edges)
        {
            assert(!e.is_broken);
            auto bar = (m_points[e.p0] + m_points[e.p1])*0.5;
            e.pb = m_points.size();
            m_points.push_back(bar);
            e.is_broken = true;
        }
        
        /* refine the triangles */
        size_t triangles_to_process = m_triangles.size();
        m_edges.reserve( m_edges.size() + triangles_to_process*12 );
        auto edges_end = m_edges.end();
        for (size_t i = 0; i < triangles_to_process; i++)
        {
            /* remove the original triangle */
            triangle t = m_triangles.front();
            m_triangles.pop_front();
            
            /* find the edges of the triangle */
            auto t_e0 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p[0], t.p[1]));
            assert(t_e0.is_broken);
            
            auto t_e1 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p[1], t.p[2]));
            assert(t_e1.is_broken);
            
            auto t_e2 = *std::lower_bound(m_edges.begin(), edges_end, edge(t.p[0], t.p[2]));
            assert(t_e2.is_broken);
            
            assert(t_e0.p1 == t_e1.p0);
            assert(t_e1.p1 == t_e2.p1);
            assert(t_e0.p0 == t_e2.p0);
            
            /* compute the edges of the new four triangles */
            /* first triangle */
            m_edges.push_back( edge(t_e0.p0, t_e0.pb, t_e0.boundary_id, t_e0.is_boundary) );
            m_edges.push_back( edge(t_e0.pb, t_e2.pb, false) );
            m_edges.push_back( edge(t_e0.p0, t_e2.pb, t_e2.boundary_id, t_e2.is_boundary) );
            triangle t0(t_e0.p0, t_e0.pb, t_e2.pb);
            m_triangles.push_back(t0);
            
            /* second triangle */
            m_edges.push_back( edge(t_e0.p1, t_e0.pb, t_e0.boundary_id, t_e0.is_boundary) );
            m_edges.push_back( edge(t_e0.pb, t_e1.pb, false) );
            m_edges.push_back( edge(t_e1.p0, t_e1.pb, t_e1.boundary_id, t_e1.is_boundary) );
            triangle t1(t_e0.p1, t_e0.pb, t_e1.pb);
            m_triangles.push_back(t1);
            
            /* third triangle */
            m_edges.push_back( edge(t_e1.p1, t_e1.pb, t_e1.boundary_id, t_e1.is_boundary) );
            m_edges.push_back( edge(t_e1.pb, t_e2.pb, false) );
            m_edges.push_back( edge(t_e2.p1, t_e2.pb, t_e2.boundary_id, t_e2.is_boundary) );
            triangle t2(t_e1.p1, t_e1.pb, t_e2.pb);
            m_triangles.push_back(t2);
            
            /* fourth triangle */
            triangle t3(t_e0.pb, t_e1.pb, t_e2.pb);
            m_triangles.push_back(t3);
            
            coarse_to_fine_int[t] = std::array<triangle,4>{t0,t1,t2,t3};
        }
        
        /* we don't need the broken edges anymore, discard them */
        auto edge_is_broken = [](const edge& e) -> bool { return e.is_broken; };
        std::remove_if(m_edges.begin(), m_edges.end(), edge_is_broken);
        
        /* sort the edges to allow fast lookups */
        disk::priv::sort_uniq(m_edges);
    }
    
    
    void convert_simplicial_to_internal(const mesh_type& msh)
    {
        m_points.clear();
        m_edges.clear();
        m_triangles.clear();
        
        m_points.resize( msh.points_size() );
        std::copy(msh.points_begin(), msh.points_end(), m_points.begin());
        
        auto bs = msh.backend_storage();
        m_edges.reserve( bs->edges.size() );
        
        for (auto itor = bs->edges.begin(); itor != bs->edges.end(); itor++)
            m_edges.push_back( edge(msh, *itor) );
        std::sort(m_edges.begin(), m_edges.end());
        
        for (auto itor = bs->surfaces.begin(); itor != bs->surfaces.end(); itor++)
            m_triangles.push_back( triangle(msh, *itor) );
    }
    
    mesh_type convert_internal_to_simplicial(void)
    {
        mesh_type refined_mesh;
        auto bs_refined = refined_mesh.backend_storage();
        
        bs_refined->points = std::move(m_points);
        bs_refined->nodes.reserve(bs_refined->points.size());
        for (size_t i = 0; i < bs_refined->points.size(); i++)
        {
            auto point_id = point_identifier<2>(i);
            bs_refined->nodes.push_back(typename mesh_type::node_type({point_id}));
        }
        bs_refined->edges.resize( m_edges.size() );
        bs_refined->boundary_info.resize( m_edges.size() );
        
        for (size_t i = 0; i < m_edges.size(); i++)
        {
            typedef point_identifier<2>    point_id_type;
            
            auto p0 = point_id_type(m_edges[i].p0);
            auto p1 = point_id_type(m_edges[i].p1);
            
            bs_refined->edges[i] = mesh_edge_type({p0, p1});
            bs_refined->boundary_info[i] = disk::bnd_info{1, m_edges[i].is_boundary};
        }
        
        bs_refined->surfaces.resize( m_triangles.size() );
        auto tri_tr = [](const triangle& t) -> mesh_triangle_type {
            typedef point_identifier<2>    point_id_type;
            auto p0 = point_id_type(t.p[0]);
            auto p1 = point_id_type(t.p[1]);
            auto p2 = point_id_type(t.p[2]);
            
            return mesh_triangle_type({p0, p1, p2});
        };
        
        std::transform(m_triangles.begin(), m_triangles.end(), bs_refined->surfaces.begin(), tri_tr);
        std::sort(bs_refined->surfaces.begin(), bs_refined->surfaces.end());
        
        return refined_mesh;
    }
    
    std::vector<std::array<size_t, 4>>
    make_coarse_to_fine_map(const mesh_type& coarse, const mesh_type& fine)
    {
        std::vector<std::array<size_t, 4>>  ret;
        ret.resize( coarse.cells_size() );
        
        for (auto& c2fi : coarse_to_fine_int)
        {
            point_identifier<2> p0(c2fi.first.p[0]);
            point_identifier<2> p1(c2fi.first.p[1]);
            point_identifier<2> p2(c2fi.first.p[2]);
            
            mesh_triangle_type tri({p0, p1, p2});
            
            auto coarse_triangle_id = coarse.lookup(tri);
            std::array<size_t, 4> fine_tri_ids;
            
            size_t i = 0;
            for (auto& ft : c2fi.second)
            {
                point_identifier<2> q0(ft.p[0]);
                point_identifier<2> q1(ft.p[1]);
                point_identifier<2> q2(ft.p[2]);
                
                mesh_triangle_type ftri({q0, q1, q2});
                auto fine_triangle_id = fine.lookup(ftri);
                fine_tri_ids[i++] = fine_triangle_id;
            }
            
            ret[coarse_triangle_id] = fine_tri_ids;
        }
        
        coarse_to_fine_int.clear();
        
        return ret;
    }
    
    mesh_type
    build_next(const mesh_type& msh)
    {
        convert_simplicial_to_internal(msh);
        refine();
        auto next_mesh = convert_internal_to_simplicial();
        
        return next_mesh;
    }
    
    void
    build_hierarchy(const mesh_type& initial_mesh, size_t levels)
    {
        meshes.clear();
        coarse_to_fine.clear();
        
        meshes.push_back(initial_mesh);
        
        for (size_t i = 0; i < levels; i++)
        {
            std::cout << "Building level " << i+1 << " mesh. H = ";
            auto next_mesh = build_next( meshes[i] );
            std::cout << average_diameter(next_mesh) << std::endl;
            //next_mesh.statistics();
            meshes.push_back( next_mesh );
            auto c2f = make_coarse_to_fine_map(meshes[i], meshes[i+1]);
            coarse_to_fine.push_back( std::move(c2f) );
        }
        
        m_edges.clear();
        m_triangles.clear();
        assert (meshes.size() == coarse_to_fine.size()+1);
    }
    
    bool
    is_inside(const mesh_type& msh, const mesh_triangle_type& t, const point_type& pt) const
    {
        auto pts = points(msh, t);
        
        auto v0 = (pts[1] - pts[0]).to_vector();
        auto v1 = (pts[2] - pts[0]).to_vector();
        auto v2 = (pt - pts[0]).to_vector();
        
        auto dot00 = v0.dot(v0);
        auto dot01 = v0.dot(v1);
        auto dot02 = v0.dot(v2);
        auto dot11 = v1.dot(v1);
        auto dot12 = v1.dot(v2);
        
        auto invDenom = 1. / (dot00 * dot11 - dot01 * dot01);
        auto u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        auto v = (dot00 * dot12 - dot01 * dot02) * invDenom;
        
        auto thresh = 0.0;

        return (u >= (0-thresh)) && (v >= (0-thresh)) && (u + v <= (1+thresh));
    }
    
    size_t
    locate_recursive(const point_type& pt, size_t cur_level, size_t max_level, size_t tri_num)
    {
        if (max_level > meshes.size())
            throw std::invalid_argument("inexistent mesh level");
        
        if (cur_level > max_level)
            return tri_num;
        
        auto current_mesh = meshes[cur_level];
        auto triangles = coarse_to_fine.at(cur_level-1).at(tri_num);
        
        for (auto& t : triangles)
        {
            auto cl = *std::next(current_mesh.cells_begin(), t);
            if (is_inside(current_mesh, cl, pt))
                return locate_recursive(pt, cur_level+1, max_level, t);
        }

        throw std::invalid_argument("point not found (s)");
    }
    
public:
    mesh_hierarchy(const mesh_type& initial_mesh, size_t num_levels)
    {
        build_hierarchy(initial_mesh, num_levels);
        //dump_to_matlab(meshes[2], "mesh2.m");
    }
    
    size_t
    locate_point(const point_type& pt, size_t level)
    {
        size_t tri_num;
        for (tri_num = 0; tri_num < meshes[0].cells_size(); tri_num++)
        {
            auto cl = *std::next(meshes[0].cells_begin(), tri_num);
            if ( is_inside(meshes[0], cl, pt) )
            {
                if (level > 0)
                    return locate_recursive(pt, 1, level, tri_num);
                else
                    return tri_num;
            }
        }
        
        
        throw std::invalid_argument("point not found");
    }

    auto meshes_begin() const { return meshes.begin(); }
    auto meshes_end()   const { return meshes.end(); }
    auto meshes_size()  const { return meshes.size(); }
};


template<typename T>
void
dump_netgen_format(const disk::simplicial_mesh<T,2>& msh, const std::string& filename)
{
    std::ofstream ofs(filename);
    
    ofs << msh.points_size() << std::endl;
    for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        ofs << (*itor).x() << " " << (*itor).y() << std::endl;
    
    ofs << msh.cells_size() << std::endl;
    for (auto itor = msh.cells_begin(); itor != msh.cells_end(); itor++)
    {
        auto cell = *itor;
        auto ptids = cell.point_ids();
        ofs << "1 ";
        for (auto& pi : ptids)
            ofs << pi+1 << " ";
        ofs << std::endl;
    }
    
    ofs << msh.boundary_faces_size() << std::endl;
    for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
    {
        auto face = *itor;
        auto ptids = face.point_ids();
        ofs << "1 ";
        for (auto& pi : ptids)
            ofs << pi+1 << " ";
        ofs << std::endl;
    }
    
    ofs.close();
}


} //namespace disk
