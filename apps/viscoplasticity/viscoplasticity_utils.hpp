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

 #include "bases/bases.hpp"
 #include "quadratures/quadratures.hpp"
 #include "methods/hho"

 using namespace revolution;

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
    typedef point<typename Mesh::coordinate_type, 2> point_type;
    //Warningn this may work only on convex polygons
    auto h = barycenter(msh, cl);
    auto sorted_pts = pts;

    std::sort(sorted_pts.begin(),
              sorted_pts.end(),
              [&](const point<typename Mesh::coordinate_type, 2>& va, const point_type& vb) {
                  auto theta_a = to_angle(va, h);
                  auto theta_b = to_angle(vb, h);
                  return (theta_a < theta_b);
              });
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
wn_PnPoly(const Mesh& msh, const typename Mesh::cell& cl, const point<typename Mesh::coordinate_type, 2>& P)
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
plot_over_line(const Mesh&                                                                                          msh,
               const std::pair<point<typename Mesh::coordinate_type, 2>,
               const point<typename Mesh::coordinate_type, 2>>& e,
               const Matrix<typename Mesh::coordinate_type, Dynamic, 1>&                                            vec,
               const size_t       cell_degree,
               const std::string& filename)
{

    typedef Matrix<typename Mesh::coordinate_type, Dynamic, 1>   vector_type;
    typedef point<typename Mesh::coordinate_type, 2>             point_type;

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
                auto cell_ofs = disk::priv::offset(msh, cl);

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
compute_discontinuous_velocity(const Mesh&                                           msh,
                               const dynamic_vector<typename Mesh::coordinate_type>& sol,
                               const typename disk::hho_degree_info&                 hdi,
                               const std::string&                                    filename)
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
        auto cell_ofs = disk::priv::offset(msh, cl);
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
quiver( const Mesh& msh, const dynamic_vector<T>& sol, const Assembler& assembler,
        const typename disk::hho_degree_info & di, const std::string& filename)
{
    std::ofstream ofs(filename);

    if (!ofs.is_open())
        std::cout << "Error opening errors "<<std::endl;

    auto cbs = disk::vector_basis_size(di.cell_degree(),
                                        Mesh::dimension, Mesh::dimension);
    for (auto& cl: msh)
    {
        auto cell_ofs = disk::priv::offset(msh, cl);
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
