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

#include "solvers/solver.hpp"

#include "output/gmshConvertMesh.hpp"
#include "output/gmshDisk.hpp"
#include "output/postMesh.hpp"

#include "output/silo.hpp"


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
                            [&](const point<typename Mesh::scalar_type,2>& va, const point_type& vb )
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

    typedef  typename Mesh::scalar_type scalar_type;

    std::vector< point<scalar_type,2>> svts(vts.size()+1);
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
                auto P_on_left = isLeft( svts[i], svts[i+1], P) > scalar_type(0);
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
                auto P_on_right = isLeft( svts[i], svts[i+1], P) < scalar_type(0);
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
                                point<typename Mesh::scalar_type,2>>   & e,
                const Matrix<typename Mesh::scalar_type, Dynamic, 1> & sol,
                const typename revolution::hho_degree_info& hdi,
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
                auto cbs   = revolution::vector_basis_size(hdi.cell_degree(), dim, dim);
                auto cell_ofs = revolution::priv::offset(msh, cl);

                vector_type s = sol.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb = revolution::make_vector_monomial_basis(msh,
                                                        cl, hdi.cell_degree());
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
    typedef typename Mesh::scalar_type scalar_type;
    const size_t cell_degree   = hdi.cell_degree();
    const size_t cbs = revolution::vector_basis_size(cell_degree,
                                    Mesh::dimension, Mesh::dimension);
    // compute mesh for post-processing
    disk::PostMesh<mesh_type> post_mesh = disk::PostMesh<mesh_type>(msh);
    gmsh::Gmesh gmsh    = disk::convertMesh(post_mesh);
    auto        storage = post_mesh.mesh().backend_storage();

    const static_vector<scalar_type, Mesh::dimension> vzero =
                            static_vector<scalar_type, Mesh::dimension>::Zero();

    const size_t nb_nodes(gmsh.getNumberofNodes());

    // first(number of data at this node), second(cumulated value)
    std::vector<std::pair<size_t, static_vector<scalar_type, Mesh::dimension>>>
                                    value(nb_nodes, std::make_pair(0, vzero));

    size_t cell_i = 0;

    for (auto& cl : msh)
    {
        auto cell_ofs = revolution::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> x = sol.block(cell_ofs * cbs, 0, cbs, 1);

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
        const static_vector<scalar_type, Mesh::dimension> depl_avr =
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

template<typename Mesh>
auto
run_stokes(const Mesh& msh, size_t degree, bool use_sym_grad = true)
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;
    typedef typename mesh_type::point_type  point_type;

    typedef dynamic_matrix<scalar_type>     matrix_type;
    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    auto rhs_fun = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };

    auto wall = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>::Zero();
    };
    auto flow_inlet = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>{1,0};
    };
    auto free_output = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
        return Matrix<scalar_type, 2, 1>{0,0};
    };


    typename revolution::hho_degree_info hdi(degree, degree);
    boundary_type bnd(msh);

    bnd.addNeumannBC(10, 1, free_output);
    //bnd.addDirichletBC(9, 1, wall);

    bnd.addDirichletBC(0, 3, flow_inlet);

    bnd.addDirichletBC(0, 2, wall);
    bnd.addDirichletBC(0, 4, wall);
    bnd.addDirichletBC(0, 5, wall);
    bnd.addDirichletBC(0, 6, wall);
    bnd.addDirichletBC(0, 7, wall);
    bnd.addDirichletBC(0, 8, wall);

    auto assembler = revolution::make_stokes_assembler(msh, hdi, bnd);

    auto viscosity = scalar_type(100);
    scalar_type factor = (use_sym_grad)? 2. * viscosity : viscosity;

    for (auto cl : msh)
    {
        auto gr = revolution::make_hho_stokes(msh, cl, hdi, use_sym_grad);
        Matrix<scalar_type, Dynamic, Dynamic> stab;
        stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, hdi);
        auto dr = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, hdi);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, hdi.cell_degree());
        auto rhs = make_rhs(msh, cl, cb, rhs_fun);
        assembler.assemble(msh, cl, factor * (gr.second + stab), -dr.second, rhs);
    }

    assembler.finalize();

    //dump_sparse_matrix(assembler.LHS, "stokes.txt");

    size_t systsz = assembler.LHS.rows();
    size_t nnz = assembler.LHS.nonZeros();

    dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);

    disk::solvers::pardiso_params<scalar_type> pparams;
    mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);

    std::ofstream ofs("velocity.dat");
    for (auto& cl : msh)
    {
        auto cbs   = revolution::vector_basis_size(degree, Mesh::dimension, Mesh::dimension);
        auto fbs   = revolution::vector_basis_size(degree, Mesh::dimension - 1, Mesh::dimension);
        auto cbs_B = revolution::scalar_basis_size(degree, Mesh::dimension);

        auto cell_ofs = revolution::priv::offset(msh, cl);
        auto num_other_faces = assembler.num_assembled_faces();
        auto offset_B = cbs * msh.cells_size() + fbs * num_other_faces + cbs_B * cell_ofs;

        Matrix<scalar_type, Dynamic, 1> s = sol.block(cell_ofs * cbs, 0, cbs, 1);
        Matrix<scalar_type, Dynamic, 1> press = sol.block( offset_B, 0, cbs_B, 1);
        auto bar = barycenter(msh, cl);
        ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " ";
        ofs << press(0) << std::endl;
        #if 0
        auto fcs = faces(msh, cl);
        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto fc = fcs[i];
            if (msh.is_boundary(fc))
            {
                auto face_offset = assembler.global_face_offset(msh, fc);
                Matrix<scalar_type, Dynamic, 1> s = sol.block(face_offset, 0, fbs, 1);

                auto bar = barycenter(msh, fc);
                ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " ";
                std::cout << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1)<< " " << std::endl;
                ofs << press(0) << std::endl;
            }
        }
        #endif
    }

    ofs.close();

    compute_discontinuous_velocity( msh, sol, hdi, "flow_over_cylinder2d.msh");

    return;

}

//#if 0
int main(int argc, char **argv)
{
    using RealType = double;
    bool use_sym_grad = false;

    char    *filename       = nullptr;
    int ch;
    size_t degree = 1;

    while ( (ch = getopt(argc, argv, "k:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1)
    {
        std::cout << "Please specify a 2D mesh" << std::endl;

        return 0;
    }

    filename = argv[0];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_stokes(msh, degree, use_sym_grad);
        std::cout << "fini" << std::endl;
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        run_stokes(msh, degree, use_sym_grad);
    }

    return 0;
}
//#endif
