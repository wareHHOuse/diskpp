/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Karol Cascavita (C) 2018                     karol.cascavita@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */
#include <iostream>
#include <iomanip>
#include <regex>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"


using namespace disk;



enum vector_problem
{
    DRIVEN,
    COUETTE,
    VANE
};

std::ostream& operator<<(std::ostream & os, const vector_problem & p)
{
    os << "Problem : "<<std::endl;
    switch (p)
    {
        case DRIVEN : os << "* problem : DRIVEN"<< std::endl;         break;
        case COUETTE: os << "* problem : COUETTE"<< std::endl;        break;
        case VANE   : os << "* problem : VANE"<< std::endl;           break;
        default:
            os << "* problem : NOT SPECIFIED"<< std::endl;
            exit(1);
    }
    return os;
}

template<typename T, typename ProblemType>
struct bingham_data
{
     bingham_data(): Lref(1.), Vref(1.), Bn(0.), mu(1.), alpha(1.), yield(0.),
                         f(1.),  problem(DRIVEN)
     {}
     T f;                    //WK: Cuidado porque f deberia ser el valor externo de la fuente.
     T Lref;                 /* Charactetistic length */
     T Vref;                 /* Reference velocity */
     T Bn;                   /* Bingham number */
     T mu;                   /* viscosity */
     T alpha;                /* ALG parameter*/
     T yield;
     ProblemType problem;
     std::string info;

     friend std::ostream& operator<<(std::ostream& os, const bingham_data<T, ProblemType>& p)
     {
         os << p.problem <<std::endl;
         os << "Bingham data: "<<std::endl;
         os << "* f      : "<< p.f<< std::endl;
         os << "* Lref   : "<< p.Lref<< std::endl;
         os << "* Vref   : "<< p.Vref<< std::endl;
         os << "* mu     : "<< p.mu<< std::endl;
         os << "* Bi     : "<< p.Bn<< std::endl;
         os << "* yield  : "<< p.yield<< std::endl;
         os << "* alpha  : "<< p.alpha<< std::endl;
         os << "* info  : "<< p.info<< std::endl;

         return os;
     }

};

template<template<typename, size_t, typename> class Mesh,
      typename T, typename Storage>
void
renumber_boundaries(Mesh<T,2,Storage>& msh, T a = 1., T b = 1., T c = 0., T d = 0.)
{
    /*--------------------------------------------------------------------------
    *  Unification of boundary labels for a square domain
    *   Netgen     __1__          Medit     __1__              __1__
    *         4   |     | 2                |     |    =>    4 |     | 2
    *             |_____|                  |_____|            |____ |
    *                3                        2                  3
    *-------------------------------------------------------------------------*/

    auto storage = msh.backend_storage();

    for(size_t i = 0; i < msh.faces_size(); i++)
    {
        auto edge = *std::next(msh.faces_begin(), i);
        auto bar = barycenter(msh, edge);

        auto is_close_to = [](T val, T ref) -> bool {
            T eps = 1e-7;
            return std::abs(val - ref) < eps;
        };

        if (!storage->boundary_info.at(i).is_boundary())
            continue;

        size_t bid = 42;

        if ( is_close_to(bar.y(), a) ) bid = 1;
        if ( is_close_to(bar.x(), b) ) bid = 2;
        if ( is_close_to(bar.y(), c) ) bid = 3;
        if ( is_close_to(bar.x(), d) ) bid = 4;

        if (bid == 42)
            throw std::logic_error("Can not locate the edge");

        storage->boundary_info.at(i).id(bid);
    }
}


template<typename Mesh>
auto
make_scalar_solution_offset(const Mesh& msh, const hho_degree_info& hdi)
{
    auto cbs = scalar_basis_size( hdi.cell_degree(), Mesh::dimension);
    auto fbs = scalar_basis_size( hdi.face_degree(), Mesh::dimension-1);
    auto map = std::vector<size_t>(msh.cells_size());

    size_t sum = 0;
    size_t cl_id = 0;

    for(auto cl : msh)
    {
        map.at(cl_id++) = sum;
        auto cell_total_dofs = cbs + howmany_faces(msh, cl) * fbs;
        sum += cell_total_dofs;
    }
    return std::make_pair(sum, map);
}



template<typename Mesh>
eigen_compatible_stdvector<Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>>
tensor_initialize(const Mesh& msh, const size_t quad_degree, const size_t tsr_degree)
{
    using T = typename Mesh::coordinate_type;
    eigen_compatible_stdvector<Matrix<T, Dynamic, Dynamic>> ret(msh.cells_size());

    auto sbs = sym_matrix_basis_size(tsr_degree, Mesh::dimension, Mesh::dimension);
    auto cell_i = 0;
    for(auto cl : msh)
    {
        auto qps = integrate(msh, cl, quad_degree);
        ret.at(cell_i++) = Matrix<T, Dynamic, Dynamic>::Zero(sbs, qps.size());
    }
    return ret;
}

#if 0
template<typename Mesh>
size_t
tensor_quad_points_size(const Mesh& msh, const size_t quad_degree);
{}

 template<typename Mesh>
 class tensors_at_quad_pts_utils
 {
     std::vector<std::pair<size_t, size_t>> offsets_vec;
     size_t m_total_quads, m_quad_degree;

 public:
     tensors_at_quad_pts_utils(){};

     tensors_at_quad_pts_utils(const Mesh msh, const size_t quad_degree):
         m_quad_degree(quad_degree)
    {
         offsets_vec = std::vector<std::pair<size_t, size_t>>(msh.cells_size());

         size_t init_quads = 0;
         size_t cl_id = 0;

        for(auto cl : msh)
        {
            auto cell_quadpoints = integrate(msh, cl, quad_degree);
            auto number_quads = cell_quadpoints.size();
            offsets_vec.at(cl_id++) = std::make_pair(init_quads, number_quads);
            init_quads += cell_quadpoints.size();
        }
        m_total_quads = init_quads;
    }

     auto
     quad_degree()
     {
         return m_quad_degree;
     }
     auto
     num_total_quad_points()
     {
         return m_total_quads;
     }

     std::vector<std::pair<size_t, size_t>>
     offsets_vector()
     {
         return offsets_vec;
     }
 };
#endif

 template< typename T>
 std::string
 tostr(const T & var)
 {
     std::ostringstream  ostr;
     ostr << var;
     return ostr.str();
 }

template <typename T>
void
save_data(const Matrix< T, Dynamic, 1>  & vec,
          const std::string     & filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(size_t i = 0; i < vec.rows(); i++)
        ofs << vec(i) <<std::endl;
    ofs.close();
};
template <typename T>
void
save_data(const std::vector<T>  & vec,
          const std::string     & filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    ofs << vec.size()<<std::endl;

    for(auto& v :  vec)
        ofs << v <<std::endl;
    ofs.close();
};
template<typename Mesh>
void
save_coords(const Mesh& msh,
            const std::string & filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
        std::cout << "Error opening file"<<std::endl;

    for(auto& cell :  msh)
    {
        auto pts = points(msh, cell);
        for(auto p : pts)
            ofs << p.x() <<" ";
        for(auto p : pts)
            ofs << p.y() <<" ";
        ofs << std::endl;
    }
    ofs.close();
};


// Return polar angle of p with respect to origin o
template<typename T>
auto
to_angle(const point<T,2>& p, const point<T,2>& o)
{
  return atan2(p.y() - o.y(), p.x() - o.x());
}
template<typename Mesh, typename Points>
auto
sort_by_polar_angle(const Mesh & msh,
                    const typename Mesh::cell &  cl,
                    const Points& pts)
{
    typedef point<typename Mesh::coordinate_type,2> point_type;
    //Warningn this may work only on convex polygons
    auto h = barycenter(msh, cl);
    auto sorted_pts = pts;

    std::sort( sorted_pts.begin(), sorted_pts.end(),
                            [&](const point<typename Mesh::coordinate_type,2>& va,
                                 const point_type& vb )
        {
            auto theta_a = to_angle(va, h);
            auto theta_b = to_angle(vb, h);
            return (theta_a < theta_b);
        }
    );
    return sorted_pts;
}
/*
* Copyright 2000 softSurfer, 2012 Dan Sunday
* This code may be freely used and modified for any purpose
* providing that this copyright notice is included with it.
* SoftSurfer makes no warranty for this code, and cannot be held
* liable for any real or imagined damage resulting from its use.
* Users of this code must verify correctness for their application.
*===================================================================
* isLeft(): tests if a point is Left|On|Right of an infinite line.
*    Input:  three points P0, P1, and P2
*    Return: >0 for P2 left of the line through P0 and P1
*            =0 for P2  on the line
*            <0 for P2  right of the line
*    See: Algorithm 1 "Area of Triangles and Polygons"
*===================================================================
*/
template<typename T>
auto
isLeft( const point<T,2>& P0, const point<T,2>& P1, const point<T,2>& P2 )
{
    auto ret = ( (P1.x() - P0.x()) * (P2.y() - P0.y())
            - (P2.x() -  P0.x()) * (P1.y() - P0.y()) );
    return ret;
}
/*===================================================================
* wn_PnPoly(): winding number test for a point in a polygon
*      Input:   P = a point,
*               vts = vertex points of a polygon vts[n+1] with vts[n]=V[0]
*      Return:  the winding number (= 0 only when P is outside)
*===================================================================
*/
template<typename Mesh>
bool
wn_PnPoly(const Mesh& msh,
          const typename Mesh::cell& cl,
          const point<typename Mesh::coordinate_type,2>& P)
{
    auto vts = points(msh, cl);

    typedef typename Mesh::coordinate_type T;

    std::vector< point<T,2>> svts(vts.size()+1);
    auto svts_temp = sort_by_polar_angle(msh, cl, vts);
    std::copy(std::begin(svts_temp), std::end(svts_temp), std::begin(svts));
    svts[vts.size()] = svts_temp[0];

    int winding_number = 0;

    for (int i = 0 ; i < svts.size() - 1; i++)
    {
        if (svts.at(i).y() <= P.y() )
        {
            auto upward_crossing = svts[i+1].y()  > P.y();
            if (upward_crossing)
            {
                auto P_on_left = isLeft( svts[i], svts[i+1], P) > T(0);
                if (P_on_left)
                    ++winding_number; //valid up intersect
            }
        }
        else
        {
            // start y > P.y (no test needed)
            auto downward_crossing = svts[i+1].y()  <= P.y();
            if (downward_crossing)
            {
                auto P_on_right = isLeft( svts[i], svts[i+1], P) < T(0);
                if (P_on_right)
                    --winding_number;  //valid down intersect
            }
        }
    }
    if(winding_number != 0)
        return true;
    return false;
}

template<typename Mesh>
void
plot_over_line(const Mesh    & msh,
                const std::pair<point<typename Mesh::coordinate_type,2>,
                                point<typename Mesh::coordinate_type,2>> & e,
                const Matrix<typename Mesh::coordinate_type, Dynamic, 1> & vec,
                const size_t cell_degree,
                const std::string & filename)
{

    typedef Matrix<typename Mesh::coordinate_type, Dynamic, 1>  vector_type;
    typedef point<typename Mesh::coordinate_type,2>             point_type;

    std::ofstream pfs(filename);
    if(!pfs.is_open())
        std::cout << "Error opening file :"<< filename <<std::endl;

    auto N = 400;
    auto h = (e.second - e.first)/N;
    size_t i = 0 ;
    auto pts = std::vector<point_type>(N);
    auto dim = Mesh::dimension;
    for(auto& p: pts)
    {
        p = e.first + i * h;

        for(auto cl: msh)
        {
            if(wn_PnPoly( msh, cl, p))
            {
                auto cbs   = disk::vector_basis_size(cell_degree, dim, dim);
                auto cell_ofs = disk::offset(msh, cl);

                vector_type s = vec.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb  = disk::make_vector_monomial_basis(msh, cl, cell_degree);
                auto phi = cb.eval_functions(p);
                vector_type vel = phi.transpose() * s;
                pfs<< p.x() << " "<< p.y() << " "<< vel(0) << " "<< vel(1)<< std::endl;
            }
        }
        i++;
    }
    pfs.close();
    return;
}

template<typename Mesh>
void
find_values_at_points(const Mesh    & msh,
                const std::vector<typename Mesh::point_type>& pts,
                const Matrix<typename Mesh::coordinate_type, Dynamic, 1> & vec,
                const size_t cell_degree,
                const std::string & filename)
{

    typedef Matrix<typename Mesh::coordinate_type, Dynamic, 1>  vector_type;
    typedef point<typename Mesh::coordinate_type,2>             point_type;

    std::ofstream pfs(filename);
    if(!pfs.is_open())
        std::cout << "Error opening file :"<< filename <<std::endl;

    auto dim = Mesh::dimension;
    for(auto& p: pts)
    {
        for(auto cl: msh)
        {
            if(wn_PnPoly( msh, cl, p))
            {
                auto cbs   = vector_basis_size(cell_degree, dim, dim);
                auto cell_ofs = disk::offset(msh, cl);

                vector_type s = vec.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb  = make_vector_monomial_basis(msh, cl, cell_degree);
                auto phi = cb.eval_functions(p);
                vector_type vel = phi.transpose() * s;
                pfs<< p.x() << " "<< p.y() << " "<< vel(0) << " "<< vel(1)<< std::endl;
            }
        }
    }
    pfs.close();
    return;
}

template<typename Mesh>
void
compute_discontinuous_velocity(const Mesh& msh,
                        const dynamic_vector< typename Mesh::coordinate_type>& sol,
                        const typename disk::hho_degree_info& hdi,
                        const std::string& filename)
{
    typedef Mesh mesh_type;
    typedef typename Mesh::coordinate_type T;
    auto dim = Mesh::dimension;

    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = disk::vector_basis_size(cell_degree,
                                    dim, dim);
    // compute mesh for post-processing
    disk::PostMesh<mesh_type> post_mesh = disk::PostMesh<mesh_type>(msh);
    gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
    auto        storage = post_mesh.mesh().backend_storage();

    const static_vector<T, Mesh::dimension> vzero =
                            static_vector<T, Mesh::dimension>::Zero();

    const size_t nb_nodes(gmsh.getNumberofNodes());

    // first(number of data at this node), second(cumulated value)
    std::vector<std::pair<size_t, static_vector<T, Mesh::dimension>>>
                                    value(nb_nodes, std::make_pair(0, vzero));

    size_t cell_i = 0;

    for (auto& cl : msh)
    {
        auto cell_ofs = disk::offset(msh, cl);
        Matrix<T, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

        const auto cell_nodes = post_mesh.nodes_cell(cell_i);
        auto  cbas = disk::make_vector_monomial_basis(msh, cl, cell_degree);

       // Loop on the nodes of the cell
        for (auto& point_id : cell_nodes)
        {
            const auto pt = storage->points[point_id];

            const auto phi = cbas.eval_functions(pt);
            assert(phi.rows() == cbs);
            const auto depl = disk::eval(x, phi);

            // Add displacement at node
            value[point_id].first++;
            value[point_id].second += depl;
        }
        cell_i++;
    }

    std::vector<gmsh::Data>    data;    // create data
    std::vector<gmsh::SubData> subdata; // create subdata
    data.reserve(nb_nodes);             // data has a size of nb_node

    // Compute the average value and save it
    for (size_t i_node = 0; i_node < value.size(); i_node++)
    {
        const static_vector<T, Mesh::dimension> depl_avr =
                            value[i_node].second / double(value[i_node].first);

        const gmsh::Data tmp_data(i_node + 1, disk::convertToVectorGmsh(depl_avr));
        data.push_back(tmp_data);
    }

    // Create and init a nodedata view
    gmsh::NodeData nodedata(3, 0.0, "depl_node_cont", data, subdata);
    // Save the view
    nodedata.saveNodeData(filename, gmsh);

    return;
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename Mesh, typename T, typename Assembler>
void
quiver( const Mesh& msh, const disk::dynamic_vector<T>& sol, const Assembler& assembler,
        const typename disk::hho_degree_info & di, const std::string& filename)
{
    std::ofstream ofs(filename);

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    auto cbs = disk::vector_basis_size(di.cell_degree(),
                                        Mesh::dimension, Mesh::dimension);
    for (auto& cl: msh)
    {
        auto cell_ofs = disk::offset(msh, cl);
        Matrix<T, Dynamic, 1>  s = assembler.take_velocity(msh, cl, sol);

        auto cb  = disk::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto qps =  disk::integrate(msh, cl, 2 * di.face_degree());
        for (auto& qp : qps)
        {
            Matrix<T, Dynamic, 2>  phi = cb.eval_functions(qp.point());
            Matrix<T,  1, Dynamic>  st = (s.block(0,0,cbs,1)).transpose();
            Matrix<T, Dynamic, 1>  vt =  st * phi;

            ofs << qp.point().x() << "  "<< qp.point().y()<< "   ";
            ofs << vt(0) << "  "<< vt(1) << "    "<< std::endl;
        }
    }
    ofs.close();
}

template<typename T, size_t DIM, typename Storage>
class post_processing_base
{
    static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
};

template<typename T, typename Storage>
class post_processing_base<T,1, Storage>
{
    public:

    typedef disk::mesh<T,1,Storage>             mesh_type;
    typedef typename mesh_type::point_type      point_type;
    typedef std::vector<point_type>             point_vector_type;

    size_t m_degree;
    size_t m_num_sub_nodes,   m_num_sub_elems;
    size_t m_total_sub_nodes, m_total_sub_elems;

    point_vector_type   m_sub_nodes;
    point_vector_type   m_vertices;

//    post_processing_base(const mesh<T,1,Storage>& msh)
    post_processing_base()
    {}

    void
    plot_nodes(size_t n, const std::vector<point_type> & vers, std::vector<point_type> & nodes)
    {
        auto p1 = vers[0];
        auto p2 = vers[1];

        if(n == 0)
        {
            nodes.push_back(p1);
            nodes.push_back(p2);
        }
        else
        {
            T h = std::abs(p1.x()-p2.x())/n;
            for(size_t i = 0; i < n + 1;  i++)
            {
                point_type p = point<T,1>({h*i});
                nodes.push_back(p);
            }
        }
    }


    std::vector<point_type>
    make_vtk_points(const mesh_type& msh,const size_t degree)
    {
        m_degree = degree;
        //if(m_degree == 0)
        //    m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = m_degree+1;
        m_num_sub_elems   = 1;
        m_total_sub_nodes = 2 * m_num_sub_nodes * msh.cells_size();
        m_total_sub_elems = 2 * msh.cells_size();
        m_sub_nodes.reserve(m_num_sub_nodes);
        m_vertices.reserve(2); /*WK: do it without reserve*/
        m_vertices[0]     = point<T,1>({0.});
        m_vertices[1]     = point<T,1>({1.});

        plot_nodes(m_degree, m_vertices, m_sub_nodes);
        std::vector<point_type> test_points(m_total_sub_nodes);
        size_t cont = 0;

        for(auto& cl : msh)
        {
            auto cl_faces = faces(msh, cl);
            auto bar = barycenter(msh, cl);
            auto pts = points(msh, cl);
            auto h   = std::abs(pts[1].x() - pts[0].x())/2.;
            size_t ifc = 0;
            for(auto& fc:cl_faces)
            {
                for(size_t  i = 0; i < m_num_sub_nodes; i++)
                {
                    auto p     = m_sub_nodes[i];
                    auto c     = point<T,1> ({h*ifc});
                    int  idx   = cont * m_num_sub_nodes + i;

                    test_points[idx] = h * p + pts[0] + c;
                }
                ++cont;
                ++ifc;
            }
        }
        return test_points;
    }
};

template<typename T, typename Storage>
class post_processing_base<T,2,Storage>
{
public:
    typedef disk::mesh<T,2,Storage>                       mesh_type;
    typedef typename mesh_type::point_type          point_type;
    typedef std::vector<point_type>                 point_vector_type;

    size_t m_degree;
    size_t m_num_sub_nodes,   m_num_sub_elems;
    size_t m_total_sub_nodes, m_total_sub_elems;
    point_vector_type   m_sub_nodes;
    point_vector_type   m_vertices;


    //post_processing_base(const mesh<T,2,Storage>& msh)
    post_processing_base(const size_t degree): m_degree(degree)
    { }

    void
    plot_nodes(const size_t n, const std::vector<point_type>& vers, std::vector<point_type> & nodes)
    {
        /*WK: check for degree 0*/
        auto p1 = vers[0];
        auto p2 = vers[1];
        auto p3 = vers[2];

        if(n == 0)
        {
            nodes.push_back(p1);
            nodes.push_back(p2);
            nodes.push_back(p3);
        }
        else
        {
            auto p4 = (p1 + p2)/2.;
            auto p5 = (p2 + p3)/2.;
            auto p6 = (p1 + p3)/2.;

            point_vector_type vertices_T1 = point_vector_type({p1, p4, p6});
            point_vector_type vertices_T2 = point_vector_type({p4, p5, p6});
            point_vector_type vertices_T3 = point_vector_type({p4, p2, p5});
            point_vector_type vertices_T4 = point_vector_type({p6, p5, p3});

            plot_nodes(n-1, vertices_T1, nodes);
            plot_nodes(n-1, vertices_T2, nodes);
            plot_nodes(n-1, vertices_T3, nodes);
            plot_nodes(n-1, vertices_T4, nodes);
        }
    }

    std::vector<point_type>
    make_vtk_points(const mesh_type& msh, const size_t degree)
    {
        m_degree = degree;
        //if(m_degree == 0)
        //    m_degree = 1; /* Same treatment as linear functions*/

        m_num_sub_nodes   = 3 * std::pow(4,m_degree);
        m_num_sub_elems   = std::pow(4,m_degree);
        m_total_sub_nodes = m_num_sub_nodes*(msh.boundary_faces_size() +  2*msh.internal_faces_size() );
        m_total_sub_elems = m_num_sub_elems*(msh.boundary_faces_size() +  2*msh.internal_faces_size() );
        m_sub_nodes.reserve(m_num_sub_nodes);
        m_vertices.reserve(3);
        m_vertices[0]   =   point<T,2>({0.,0.});
        m_vertices[1]   =   point<T,2>({1.,0.});
        m_vertices[2]   =   point<T,2>({0.,1.});

        plot_nodes(m_degree, m_vertices, m_sub_nodes);

        std::vector<point_type> test_points(m_total_sub_nodes);
        size_t cont = 0;

        for(auto& cl : msh)
        {
            auto cl_faces = faces(msh, cl);
            auto bar = barycenter(msh, cl);
            for(auto& fc : cl_faces)
            {
                auto pts = points(msh, fc);
                for(size_t  i = 0; i < m_num_sub_nodes; i++)
                {
                    auto p    = m_sub_nodes[i];
                    int idx   = cont*m_num_sub_nodes;
                    test_points[idx+i] = ((bar + (pts[1] - bar) * p.x() ) + (pts[0] - bar) * p.y());
                }
            ++cont;
            }
        }
        return test_points;
    }
};


template<typename Mesh>
class paraview
{};

template<typename T, typename Storage>
class paraview<disk::mesh<T,2, Storage>>: public post_processing_base<T, 2,Storage>
{
public:

    typedef disk::mesh<T,2, Storage>                    mesh_type;
    typedef typename mesh_type::point_type              point_type;

    typedef post_processing_base<T,2,Storage>         pp_base;
    typedef disk::dynamic_vector<T>                           vector_type;
    typedef disk::dynamic_matrix<T>                           matrix_type;
    typedef std::vector<point_type>                     point_vector_type;

    point_vector_type sub_nodes;
    point_vector_type vertices;

    std::vector<point_type> test_points;
    size_t num_sub_elems;
    size_t num_sub_nodes;
    size_t total_sub_nodes;
    size_t total_sub_elems;
    size_t m_degree;
    size_t DIM;

    paraview(const size_t degree):
        pp_base(degree)
    {
        m_degree = degree;
        DIM = 2;
    }

    struct sizes
    {
        sizes(){};
        size_t num_nodes;
        size_t num_elems;
        size_t num_scalars;
        size_t num_vectors;
    };

    void vtk_vector(std::ofstream &ofs,
                    const mesh_type& msh,
                    const std::string& name,
                    const disk::dynamic_vector<T> & vec)
    {
        ofs << "VECTORS  Vector  double"<< std::endl;

        size_t cl_cont = 0;
        size_t fc_cont = 0;

        auto cbs_scalar = scalar_basis_size(m_degree, mesh_type::dimension);
        auto cbs = vector_basis_size(m_degree, mesh_type::dimension, mesh_type::dimension);

        for(auto& cl : msh)
        {
            auto cb = make_vector_monomial_basis(msh, cl, m_degree);

            auto cl_id = msh.lookup( cl);
            auto fcs   = faces(msh, cl);
            vector_type v = vec.block(cbs * cl_id, 0, cbs, 1);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto c_phi = cb.eval_functions(tp);

                    vector_type pot = vector_type::Zero(2);

                    for (size_t i = 0; i < cbs_scalar; i++)
                    {
                        pot(0)+=  v( 2 * i ) * c_phi(i);
                        pot(1)+=  v( 2 * i + 1) * c_phi(i);
                    }

                    for (size_t d = 0; d < DIM; d++)
                        ofs <<  pot(d) << " ";

                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        ofs<< 0. << " ";
                    ofs<< std::endl;
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        return;
    }

    void
    vtk_scalar( std::ofstream &ofs,
                const mesh_type& msh,
                const std::string& name,
                const disk::dynamic_vector<T> & vec)
    {}
        #if 0

        size_t cl_cont = 0;
        size_t fc_cont = 0;

        auto cbs = scalar_basis_size(m_degree, mesh_type::dimension);

        ofs << "SCALARS  Scalar  double 1"<< std::endl;
        ofs << "LOOKUP_TABLE default"<< std::endl;

        for(auto& cel : msh)
        {
            auto cb = make_scalar_monomial_basis(msh, cel, m_degree);

            vector_type dofs =  vec.at(cl_cont);

            auto fcs = faces(msh, cel);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto phi  = cb.eval_functions(tp);

                    assert(cbs == dofs.size());
                    for (size_t i = 0; i < cbs; i++)
                        pot  +=  phi[i] * dofs(i);
                    ofs << pot <<' ';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        return;
    }
    #endif

    void
    vtk_vector_magnitud( std::ofstream &ofs,
                const mesh_type& msh,
                const std::string& name,
                const disk::dynamic_vector<T> & vec)
    {
        size_t cl_cont = 0;
        size_t fc_cont = 0;

        auto cbs = scalar_basis_size(m_degree, mesh_type::dimension);

        ofs << "SCALARS  Scalar  double 1"<< std::endl;
        ofs << "LOOKUP_TABLE default"<< std::endl;

        for(auto& cel : msh)
        {
            auto cb = make_scalar_monomial_basis(msh, cel, m_degree);

            vector_type dofs =  vec.at(cl_cont);

            auto fcs = faces(msh, cel);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);

                    assert(cbs == dofs.size());
                    pot  =  std::sqrt( dofs(0) * dofs(0) + dofs(1) * dofs(1) );
                    ofs << pot <<' ';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        return;
    }

    void make_file(const mesh_type   & msh,
                const std::string & name,
                const disk::dynamic_vector<T> & vector_Th,
                const std::string & type)
    {
        test_points  = pp_base::make_vtk_points(msh, m_degree);

        total_sub_nodes = pp_base::m_total_sub_nodes;
        total_sub_elems = pp_base::m_total_sub_elems;
        num_sub_nodes   = pp_base::m_num_sub_nodes;
        sub_nodes       = pp_base::m_sub_nodes;
        vertices        = pp_base::m_vertices;

        size_t cl_cont = 0;
        size_t fc_cont = 0;
        //std::string file_name = name + "_i" + numstr;
        std::string file_name = name;

        std::vector<size_t> vtk_cell_type = {4,5,10}; /*WK: check this for others VTK types*/
        std::vector<size_t> vtk_nodes_inside_cell = {m_degree+1,3,4};/* WK: all cells are the same type */

        //vtk_writer(vtk_points, vtk_elems, vtk_scalar_data, vtk_vector_data, sz);

        /*VTK file*/
        std::string ext = ".vtk";
        std::ofstream  ofs(file_name + ext);
        ofs << "# vtk DataFile Version 3.0"<< std::endl;
        ofs << "#This file was generated by the DISK++ library"<< std::endl;
        ofs << "ASCII"<< std::endl;

        /* Data header */
        ofs << "DATASET UNSTRUCTURED_GRID\n" << std::endl;

        /* Points */
        ofs << "POINTS " << total_sub_nodes<<" double"<<std::endl;

        for(auto& cel : msh)
        {
            auto cl_faces = faces(msh, cel);

            for(auto& fc : cl_faces)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes*fc_cont;
                    auto tp   = test_points.at(idx);

                    for (size_t d = 0; d < DIM; d++)
                        ofs << tp.at(d) << " ";
                    for (size_t d = 0; d < 3 - DIM; d++) /* VTK only supports until 3D */
                        ofs<< 0. << " ";
                    ofs<< std::endl;
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        ofs<<std::endl;

        /* Cells */
        size_t el = 0;
        ofs << "CELLS " << total_sub_elems <<' ';
        ofs <<  total_sub_elems *(vtk_nodes_inside_cell[DIM-1] + 1)<< std::endl;
        for (size_t i = 0; i < total_sub_elems; i++)
        {
            ofs << vtk_nodes_inside_cell[DIM-1] << " ";
            for(size_t j=0; j < vtk_nodes_inside_cell[DIM-1]; j++, el++)
                ofs << el<< " ";
            ofs<< std::endl;
        }

        /* Types of cells*/
        ofs << "CELL_TYPES " << total_sub_elems << std::endl;
            for(size_t i = 0; i < total_sub_elems; i++)
                ofs << ' ' << vtk_cell_type[DIM-1];
        ofs << std::endl;

        /* Data */
        ofs << "POINT_DATA " << total_sub_nodes << std::endl;

        if(type == "scalar")
            vtk_scalar(ofs, msh, name, vector_Th);
        if(type == "vector")
            vtk_vector(ofs, msh, name, vector_Th);

        ofs.close();

        return;
    }
};
