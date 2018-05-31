/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"
#include "solvers/solver.hpp"

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
    typedef point<typename Mesh::scalar_type,2> point_type;
    //Warningn this may work only on convex polygons
    auto h = barycenter(msh, cl);
    auto sorted_pts = pts;

    std::sort( sorted_pts.begin(), sorted_pts.end(),
                            [&](const point<typename Mesh::scalar_type,2>& va,
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
          const point<typename Mesh::scalar_type,2>& P)
{
    auto vts = points(msh, cl);

    typedef typename Mesh::scalar_type T;

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
                const std::pair<point<typename Mesh::scalar_type,2>,
                                point<typename Mesh::scalar_type,2>> & e,
                const Matrix<typename Mesh::scalar_type, Dynamic, 1> & vec,
                const size_t cell_degree,
                const std::string & filename)
{

    typedef Matrix<typename Mesh::scalar_type, Dynamic, 1>  vector_type;
    typedef point<typename Mesh::scalar_type,2>             point_type;

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
                auto cbs   = revolution::vector_basis_size(cell_degree, dim, dim);
                auto cell_ofs = revolution::priv::offset(msh, cl);

                vector_type s = vec.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb  = revolution::make_vector_monomial_basis(msh, cl, cell_degree);
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
compute_discontinuous_velocity(const Mesh& msh,
                        const dynamic_vector< typename Mesh::scalar_type>& sol,
                        const typename revolution::hho_degree_info& hdi,
                        const std::string& filename)
{
    typedef Mesh mesh_type;
    typedef typename Mesh::scalar_type T;
    auto dim = Mesh::dimension;

    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = revolution::vector_basis_size(cell_degree,
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
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

        const auto cell_nodes = post_mesh.nodes_cell(cell_i);
        auto  cbas = revolution::make_vector_monomial_basis(msh, cl, cell_degree);

       // Loop on the nodes of the cell
        for (auto& point_id : cell_nodes)
        {
            const auto pt = storage->points[point_id];

            const auto phi = cbas.eval_functions(pt);
            assert(phi.rows() == cbs);
            const auto depl = revolution::eval(x, phi);

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
quiver( const Mesh& msh, const dynamic_vector<T>& sol, const Assembler& assembler,
        const typename revolution::hho_degree_info & di, const std::string& filename)
{
    std::ofstream ofs(filename);

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    auto cbs = revolution::vector_basis_size(di.cell_degree(),
                                        Mesh::dimension, Mesh::dimension);
    for (auto& cl: msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1>  s = assembler.take_velocity(msh, cl, sol);

        auto cb  = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto qps =  revolution::integrate(msh, cl, 2 * di.face_degree());
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


    #if 0

    template<typename T, size_t DIM, typename Storage>
    class post_processing_base
    {
        static_assert(DIM > 0 && DIM <= 3, "mesh: Allowed dimensions are 1, 2 and 3");
    };

    template<typename T, typename Storage>
    class post_processing_base<T,1, Storage>
    {
        public:

        typedef mesh<T,1,Storage>                   mesh_type;
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
        typedef mesh<T,2,Storage>                       mesh_type;
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

    void vtk_vector(std::ofstream &ofs,
                    const mesh_type& msh,
                    const std::string& name,
                    const std::vector<dynamic_vector<T>> & vec)
    {
        ofs << "VECTORS  Vector  double"<< std::endl;

        size_t cl_cont = 0;
        size_t fc_cont = 0;

        for(auto& cl : msh)
        {
            auto cl_id = msh.lookup( cl);
            auto fcs   = faces(msh, cl);
            vector_type v = vec.at(cl_id);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto c_phi = cb.eval_functions(msh, cl, tp);
                    matrix_type vec_phi = make_vectorial_matrix(c_phi, DIM);
                    vector_type dpot    = vec_phi * v;
                    for (size_t d = 0; d < DIM; d++)
                        ofs <<  dpot(d) << " ";

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
                const std::vector<dynamic_vector<T>> & vec)
    {
        size_t cl_cont = 0;
        size_t fc_cont = 0;

        ofs << "SCALARS  Scalar  double 1"<< std::endl;
        ofs << "LOOKUP_TABLE default"<< std::endl;

        for(auto& cel : msh)
        {
            vector_type dofs =  vec.at(cl_cont);

            auto fcs = faces(msh, cel);

            for(auto& fc : fcs)
            {
                for (size_t itp = 0; itp < num_sub_nodes; itp++)
                {
                    auto pot  = 0.;
                    auto idx  = itp + num_sub_nodes * fc_cont;
                    auto tp   = test_points.at(idx);
                    auto phi  = cb.eval_functions(msh, cel, tp);

                    assert(cb.size() == dofs.size());
                    for (size_t i = 0; i < cb.size(); i++)
                        pot  +=  phi[i] * dofs(i);
                    ofs << pot <<' ';
                }
                ++fc_cont;
            }
            ++cl_cont;
        }
        return;
    }
    void paraview(const mesh_type   & msh,
                  const std::string & name,
                  const std::vector<dynamic_vector<T>> & vector_Th,
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
#endif

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

enum problem_type
{
    DRIVEN,
    COUETTE,
    POISEUILLE
};

template<typename Mesh>
class augmented_lagrangian_viscoplasticity
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type T;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef dynamic_vector<T>       vector_type;
    typedef dynamic_matrix<T>       matrix_type;

    typedef Matrix<T, Mesh::dimension, Mesh::dimension>         tensor_type;
    typedef Matrix<T, Mesh::dimension, 1>                       vector2d_type;
    typedef std::function<vector2d_type (const point_type &)>   vector_funtion_type;
    typedef std::function<T   (const point_type &)>             scalar_funtion_type;

    vector_funtion_type     rhs_fun, velocity;
    vector_type             multiplier, auxiliar, auxiliar_old;

    typename revolution::hho_degree_info di;
    T             factor;
    T             viscosity;
    T             alpha;
    T             yield;
    size_t        cbs, fbs, pbs, sbs, dim;

public:
    vector_type             sol, sol_old;
    std::tuple<T, T, T>     convergence;
    bool                    use_sym_grad;

    augmented_lagrangian_viscoplasticity(const Mesh& msh,
                            const typename revolution::hho_degree_info & hdi,
                            const T& alpha_ext):
                            di(hdi), alpha(alpha_ext)
    {
        use_sym_grad = true;
        factor = (use_sym_grad)? 2. : 1.;
        viscosity = 2.;
        auto omegaExt = 2;
        T f = 1;
        T Lref = 1.;
        T Bn  =  2;
        yield =  Bn * viscosity;// * omegaExt; //* f * Lref;

        dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::sym_matrix_basis_size(di.face_degree(), dim, dim);
    };

    auto
    define_problem(const mesh_type& msh, const problem_type& problem )
    {
        boundary_type bnd(msh);

        auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>::Zero();
        };
        auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{1,0};
        };
        auto symmetryPlane = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
            return Matrix<T, Mesh::dimension, 1>{0,0};
        };

        switch (problem)
		{
            case DRIVEN:
                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    if( std::abs(p.y() - 1.) < 1.e-8 )
                        return Matrix<T, Mesh::dimension, 1>{1,0};
                    else
                        return Matrix<T, Mesh::dimension, 1>{0,0};
                };
                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                #if 0
                bnd.addDirichletBC(0, 1, movingWall);
                bnd.addDirichletBC(0, 2, wall);
                bnd.addDirichletBC(0, 3, wall);
                bnd.addDirichletBC(0, 4, wall);
                #endif
                bnd.addDirichletEverywhere(velocity);

               break;

            case POISEUILLE:
                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                        return Matrix<T, Mesh::dimension, 1>::Zero();
                    };
                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                        return Matrix<T, Mesh::dimension, 1>{1,0};
                    };

                    bnd.addDirichletBC(0, 1, wall);
                    bnd.addDirichletBC(0, 2, wall);

                    bnd.addNeumannBC(10, 3, symmetryPlane);
                    bnd.addNeumannBC(10, 4, symmetryPlane);
                break;


            case COUETTE:

                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());

                    T Rint(0.5), Rext(1.), omegaInt(2.),omegaExt(2.);
                    T r = std::sqrt(p.x()*p.x() + p.y()*p.y());

                    //1.All Solid
                    //auto u_theta_solid = (omegaExt/ omegaInt) * (r / Rint);
                    #if 0
                    //2.All Liquid
                    auto eta = Rint/Rext;
                    auto sigma_i = (2./(eta * eta - 1.)) * ( (1. - (omegaExt*eta / omegaInt))*((1 - eta)/eta)
                        + sgn((omegaInt*eta / omegaExt) - 1.) * Bn * std::log(eta) );
                    auto term_1 = u_theta_solid;
                    auto term_2 = 0.5 *(Rint * Rint) * sigma_i * r * ( 1./(Rext* Rext) - 1/(r *r)) ;
                    auto term_3 = Bn * r * std::log(Rext/r) * sgn(sigma_i);
                    u_theta = term_1 + term_2 + term_3;
                    //3. with plug
                    #endif

                    ret(0) =  -std::sin(theta) * 2.* r;
                    ret(1) =  std::cos(theta) * 2.* r;
                    return ret;
                };

                auto v1  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());
                    T r = std::sqrt(p.x()*p.x() + p.y()*p.y());

                    ret(0) =  -std::sin(theta) *  r;
                    ret(1) =  std::cos(theta) *  r;
                    return ret;
                };

                auto v2 = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());
                    T r = std::sqrt(p.x()*p.x() + p.y()*p.y());

                    //1.All Solid
                    ret(0) =  -std::sin(theta) *  2* r;
                    ret(1) =  std::cos(theta) *  2*r;
                    return ret;
                };

                bnd.addNeumannBC(10, 1, symmetryPlane);
                bnd.addNeumannBC(10, 2, symmetryPlane);
                bnd.addDirichletBC(0, 3, v2); //velocity);
                bnd.addDirichletBC(0, 4, v1);//velocity);

                break;

        }

        auto assembler = revolution::make_stokes_assembler_alg(msh, di, bnd);

        return assembler;
    }

    template<typename Assembler>
    auto
    initialize(const mesh_type& msh, const Assembler& assembler)
    {
        auto systsz = assembler.global_system_size();
        sol = vector_type::Zero(systsz);
        sol_old = vector_type::Zero(systsz);

        multiplier = vector_type::Zero(msh.cells_size() * sbs);
        auxiliar   = vector_type::Zero(msh.cells_size() * sbs);
        auxiliar_old = vector_type::Zero(msh.cells_size() * sbs);

        return;
    }

    template<typename Assembler>
    auto
    compute_errors( const mesh_type& msh,
                    const Assembler& assembler,
                    const bool& alg_finished)
    {
        T error_vel  = 0.;
        T error_temp = 0.;

        std::ofstream ofs;

        if(alg_finished)
        {
            ofs = std::ofstream("Gu_norm.data");
            if (!ofs.is_open())
                std::cout << "Error opening errors "<<std::endl;
        }

        for (auto& cl : msh)
        {
        	auto bar = barycenter(msh, cl);

            //energy error
            Matrix<T, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, sol);
            Matrix<T, Dynamic, 1> pvel = revolution::project_function(msh, cl, di, velocity);
            //Matrix<T, Dynamic, 1> svel_old =  assembler.take_velocity(msh, cl, sol_old);
            Matrix<T, Dynamic, 1> diff_vel = svel - pvel;
            auto gr = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<T, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

            Matrix<T, Dynamic, Dynamic> B = factor * (viscosity*G.second + viscosity*stab);

            error_vel += diff_vel.dot(B * diff_vel);
            auto Gu_norm = svel.dot(B * svel);
            error_temp += Gu_norm;

            if(alg_finished)
                ofs << bar.x()<< " " << bar.y() << " " << Gu_norm<< std::endl;
        }
        if(alg_finished)
            ofs.close();
        return std::make_pair(std::sqrt(error_vel), std::sqrt(error_temp));
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler,
                        const vector_type& velocity_dofs)
    {
        Matrix<T, Dynamic, 1> u_TF  = assembler.take_velocity(msh, cl, velocity_dofs);
        auto value = 1./(factor * (viscosity + alpha));
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

        auto cell_ofs    = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1> Gu = G.first * u_TF;
        Matrix<T, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);

        //Theta
        auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
        //barycenter only for k = 0; fix this for higher orders
        auto bar = barycenter(msh, cl);
        auto s_phi  = sb.eval_functions(bar);

        Matrix<T, Dynamic, 1>  theta  = stress  +  factor * alpha * Gu;
        Matrix<T, Mesh::dimension, Mesh::dimension> theta_eval = revolution::eval(theta, s_phi);
        T theta_norm  = std::sqrt((theta_eval.cwiseProduct(theta_eval)).sum());

        T theta_eigen = theta_eval.norm();
        assert(theta_norm == theta_eigen);

        //Gamma
        // A. Solid
        if(theta_norm <=  std::sqrt(2) * yield ||  std::abs(theta_norm - std::sqrt(2) * yield) < 1.e-8)
            return  Matrix<T, Dynamic, 1>::Zero(sbs);
        else  // B. liquid
            return  value * theta * (1. - std::sqrt(2) *(yield/theta_norm));
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler)
    {
        T conv_stress = 0.;
        T conv_gamma = 0.;

        auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            Matrix<T, Dynamic, 1> u_TF = assembler.take_velocity(msh, cl, sol);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = revolution::priv::offset(msh, cl);

            Matrix<T, Dynamic, 1> Gu = G.first * u_TF;
            Matrix<T, Dynamic, 1> gamma_old = auxiliar_old.block(cell_ofs *sbs, 0, sbs, 1);
            Matrix<T, Dynamic, 1> gamma = auxiliar.block(cell_ofs *sbs, 0, sbs, 1);

            Matrix<T, Dynamic, 1> diff_stress = factor * alpha * (Gu - gamma);
            Matrix<T, Dynamic, 1> diff_gamma  = alpha * (gamma - gamma_old);

            multiplier.block(cell_ofs * sbs, 0, sbs, 1) += diff_stress;

            auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, sb);

            conv_stress += diff_stress.dot(mass * diff_stress);
            conv_gamma  += diff_gamma.dot(mass * diff_gamma);
        }

        convergence = std::make_tuple(conv_stress + conv_gamma, conv_stress, conv_gamma);

        auxiliar_old = auxiliar;
        return;
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

        auto cell_ofs =  revolution::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        Matrix<T, Dynamic, 1> stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        Matrix<T, Dynamic, 1> gamma  = compute_auxiliar( msh,  cl, assembler, sol_old); //or sol at this point it's the same
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;

        vector_type str_agam = stress - factor * alpha * gamma;

        Matrix<T, Dynamic, Dynamic> mm = revolution::make_mass_matrix(msh, cl, sb);

        Matrix<T, Dynamic, 1> rhs =
                            Matrix<T, Dynamic, 1>::Zero(cbs + fbs * num_faces);

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun);

        //(stress - alpha * gamma, Gv)
        rhs -=  G.first.transpose() * mm * str_agam;

        return rhs;
    }

    template<typename Assembler>
    auto
    make_global_rhs(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_rhs();

        for (auto cl : msh)
        {
            Matrix<T, Dynamic, 1> local_rhs = make_rhs_alg(msh, cl, assembler);

            assembler.assemble_rhs(msh, cl, local_rhs);
        }
        assembler.finalize_rhs();

        return;
    }

    template<typename Assembler>
    auto
    make_global_matrix(const mesh_type& msh, Assembler& assembler)
    {
        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto G  = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto gr = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<T, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            Matrix<T, Dynamic, Dynamic> A = factor *(alpha * G.second + viscosity * stab);

            assembler.assemble_lhs(msh, cl, A, -dr.second);
        }

        assembler.finalize_lhs();

        return;
    }

    template<typename Assembler>
    auto
    run_stokes_like(const mesh_type& msh, Assembler& assembler, const size_t iter)
    {
        sol_old = sol;

        make_global_rhs(msh,assembler);

        if(iter == 0)
            make_global_matrix(msh, assembler);

        //dump_sparse_matrix(assembler.LHS, "stokes.txt");
        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        sol = vector_type::Zero(systsz);
        disk::solvers::pardiso_params<T> pparams;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

        return;
    }

    template<typename Assembler>
    void
    post_processing(const mesh_type& msh, Assembler& assembler,
                    const std::string & info,
                    const problem_type& problem)
    {
        auto dim = Mesh::dimension;
        auto rbs   = revolution::vector_basis_size(di.reconstruction_degree(), dim, dim);
        Matrix< T, Dynamic, 1> cell_sol(cbs * msh.cells_size());
        Matrix< T, Dynamic, 1> cell_rec_sol(rbs * msh.cells_size());
        Matrix< T, Dynamic, 1> press_vec(pbs * msh.cells_size());

        size_t cl_count = 0;
        for(auto cl : msh)
        {
            auto gr  = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = revolution::priv::offset(msh, cl);
            Matrix<T, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, sol);
            assert((gr.first * svel).rows() == rbs - dim);
            cell_rec_sol.block(cell_ofs * rbs + dim, 0, rbs - dim, 1) = gr.first * svel;
            cell_rec_sol.block(cell_ofs * rbs, 0, dim, 1) = svel.block(0,0, dim, 1);
            cell_sol.block(cell_ofs * cbs, 0, cbs, 1) = svel.block(0,0, cbs, 1);

            //this is only for k = 0, since there is only one velocity;
            auto bar = barycenter(msh, cl);
            Matrix<T, Dynamic, 1> spress =  assembler.take_pressure(msh, cl, sol);
            auto pb  = revolution::make_scalar_monomial_basis(msh, cl, di.face_degree());
            auto p_phi = pb.eval_functions(bar);
            press_vec(cl_count++) =  p_phi.dot(spress);
        }

        typedef point<T,2>             point_type;
        std::pair<point_type, point_type> p_x, p_y;
        auto eps = 1.e-4;
        switch(problem)
        {
            case COUETTE:
                p_x = std::make_pair(point_type({0., 0.}), point_type({0.866, 0.5}));
                p_y = std::make_pair(point_type({0.0, 0.0}), point_type({0.0, 0.0}));
                break;
            default:
                p_x = std::make_pair(point_type({0.0 + eps, 0.5 + eps}), point_type({1.0 + eps, 0.5 + eps}));
                p_y = std::make_pair(point_type({0.5 + eps, 0.0 + eps}), point_type({0.5 + eps, 1.0 + eps}));
                break;
        }
        plot_over_line(msh, p_x, cell_rec_sol, di.reconstruction_degree(), "plot_over_x_" + info + ".data");
        //plot_over_line(msh, p_y, cell_rec_sol, di.reconstruction_degree(), "plot_over_y_" + info + ".data");
        plot_over_line(msh, p_y, cell_sol, di.cell_degree(), "plot_over_y_" + info + ".data");
        compute_discontinuous_velocity( msh, cell_sol, di, "velocity_" + info +".msh");
        save_coords(msh, "Coords_"+ info + ".data");
        save_data(press_vec, "pressure_" + info + ".data");

        auto marks_sigma = std::vector<size_t>(msh.cells_size());
        auto marks_theta = std::vector<size_t>(msh.cells_size());

        size_t id = 0;
        for(auto & cl : msh)
        {
            Matrix<T, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, sol);
            auto value = 1./(factor * (viscosity + alpha));
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

            auto cell_ofs    = revolution::priv::offset(msh, cl);
            Matrix<T, Dynamic, 1> Gu = G.first * svel;
            Matrix<T, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);

            auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
            //barycenter only for k = 0; fix this for higher orders
            auto bar = barycenter(msh, cl);
            auto s_phi  = sb.eval_functions(bar);

            Matrix<T, Dynamic, 1>  theta  = stress  +  factor * alpha * Gu;
            Matrix<T, Mesh::dimension, Mesh::dimension> theta_eval = revolution::eval(theta, s_phi);
            T theta_norm  = std::sqrt((theta_eval.cwiseProduct(theta_eval)).sum());

            Matrix<T, Mesh::dimension, Mesh::dimension> sigma_eval = revolution::eval(stress, s_phi);
            T sigma_norm  = std::sqrt((sigma_eval.cwiseProduct(sigma_eval)).sum());

            if(theta_norm <=  std::sqrt(2) * yield ||  std::abs(theta_norm - std::sqrt(2) * yield) < 1.e-8)
                marks_theta.at(id++) = 1;
            if(sigma_norm <=  std::sqrt(2) * yield ||  std::abs(sigma_norm - std::sqrt(2) * yield) < 1.e-8)
                marks_sigma.at(id++) = 1;

        }

        save_data(marks_theta, "criterion_theta_"+ info + ".data");
        save_data(marks_sigma, "criterion_sigma_"+ info + ".data");

        quiver( msh, sol, assembler, di, "quiver_"+ info + ".data");
        return;
    }

};
