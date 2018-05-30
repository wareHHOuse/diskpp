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

#include "solvers/solver.hpp"

#include "output/silo.hpp"

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

template<typename Mesh>
class augmented_lagrangian_stokes
{
    typedef Mesh mesh_type;
    typedef typename mesh_type::cell        cell_type;
    typedef typename mesh_type::face        face_type;
    typedef typename mesh_type::scalar_type scalar_type;

    typedef disk::mechanics::BoundaryConditions<mesh_type> boundary_type;

    using point_type = typename mesh_type::point_type;

    typedef dynamic_vector<scalar_type>       vector_type;
    typedef dynamic_matrix<scalar_type>       matrix_type;

    typedef Matrix<scalar_type, 2, 2>         tensor_type;
    typedef Matrix<scalar_type, 2, 1>         vector2d_type;
    typedef std::function<vector2d_type (const point_type &)> vector_funtion_type;
    typedef std::function<scalar_type   (const point_type &)> scalar_funtion_type;

    vector_funtion_type     velocity, rhs_fun;
    scalar_funtion_type     pressure;
    vector_type             multiplier, auxiliar;

    typename revolution::hho_degree_info di;
    bool                    use_sym_grad;
    scalar_type             factor;
    scalar_type             viscosity;
    scalar_type             alpha;
    size_t                  cbs, fbs, pbs, sbs;

public:
    vector_type             sol, sol_old;
    scalar_type             convergence;

    augmented_lagrangian_stokes(const Mesh& msh,
                            const typename revolution::hho_degree_info & hdi,
                            const bool& use_symmetric_gradient):
                            di(hdi), use_sym_grad(use_symmetric_gradient)
    {
        #if 0
        // 1.Stokes polynomial
        rhs_fun  = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            Matrix<scalar_type, 2, 1> ret;

            scalar_type x1 = p.x();
            scalar_type x2 = x1 * x1;
            scalar_type y1 = p.y();
            scalar_type y2 = y1 * y1;

            scalar_type ax =  x2 * (x2 - 2. * x1 + 1.);
            scalar_type ay =  y2 * (y2 - 2. * y1 + 1.);
            scalar_type bx =  x1 * (4. * x2 - 6. * x1 + 2.);
            scalar_type by =  y1 * (4. * y2 - 6. * y1 + 2.);
            scalar_type cx = 12. * x2 - 12.* x1 + 2.;
            scalar_type cy = 12. * y2 - 12.* y1 + 2.;
            scalar_type dx = 24. * x1 - 12.;
            scalar_type dy = 24. * y1 - 12.;

            ret(0) = - cx * by - ax * dy + 5.* x2 * x2;
            ret(1) = + cy * bx + ay * dx + 5.* y2 * y2;

            return ret;
        };
        velocity = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            Matrix<scalar_type, 2, 1> ret;

            scalar_type x1 = p.x();
            scalar_type x2 = x1 * x1;
            scalar_type y1 = p.y();
            scalar_type y2 = y1 * y1;

            ret(0) =  x2 * (x2 - 2. * x1 + 1.)  * y1 * (4. * y2 - 6. * y1 + 2.);
            ret(1) = -y2 * (y2 - 2. * y1 + 1. ) * x1 * (4. * x2 - 6. * x1 + 2.);

            return ret;
        };
        pressure = [](const point_type& p) -> scalar_type {
            return std::pow(p.x(), 5.)  +  std::pow(p.y(), 5.)  - 1./3.;
        };
        #endif 0

        //2. poiseuille
        velocity  = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
                return Matrix<scalar_type, 2, 1>::Zero();
        };

        rhs_fun  = [](const point_type& p) -> Matrix<scalar_type, 2, 1> {
            return Matrix<scalar_type, 2, 1>{1,0};
        };

        pressure = [](const point_type& p) -> scalar_type {
            return 0.;
        };

        factor = (use_sym_grad)? 2. : 1.;
        viscosity = 1.;
        alpha = 1.;
        convergence = 0.;
        auto dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::sym_matrix_basis_size(di.face_degree(), dim, dim);
    };

    auto
    define_assembler(const mesh_type& msh)
    {
        boundary_type bnd(msh);
        // 1.Stokes polynomial
        //bnd.addDirichletEverywhere(velocity);

        //2. poiseuille
        auto wall = [](const point_type& p) -> Matrix<scalar_type,2, 1> {
            return Matrix<scalar_type,2, 1>::Zero();
        };
        auto symmetryPlane = [](const point_type& p) -> Matrix<scalar_type,2, 1> {
            return Matrix<scalar_type,2, 1>{0,0};
        };


        bnd.addDirichletBC(0, 1, wall);
        bnd.addDirichletBC(0, 2, wall);

        bnd.addNeumannBC(10, 3, symmetryPlane);
        bnd.addNeumannBC(10, 4, symmetryPlane);

        auto assembler = revolution::make_stokes_assembler(msh, di, bnd);

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

        return;
    }

    template<typename Assembler>
    auto
    compute_errors( const mesh_type& msh,
                    const Assembler& assembler)
    {
        auto dim =  Mesh::dimension;

        scalar_type error(0), error_vel(0), error_pres(0);

        for (auto& cl : msh)
        {
        	auto bar = barycenter(msh, cl);
        	Matrix<scalar_type, Dynamic, 1> p = revolution::project_function(msh, cl, di, velocity);
        	auto cell_ofs = revolution::priv::offset(msh, cl);
        	Matrix<scalar_type, Dynamic, 1> s = sol.block(cell_ofs * cbs, 0, cbs, 1);
        	Matrix<scalar_type, Dynamic, 1> diff = s - p.head(cbs);
        	auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        	Matrix<scalar_type, Dynamic, Dynamic> mm = revolution::make_mass_matrix(msh, cl, cb);
        	error += diff.dot(mm*diff);
        	//ofs << bar.x() << " " << bar.y() << " " << s(0) << " " << s(1) << std::endl;

            //pressure error
            Matrix<scalar_type, Dynamic, 1> ppres = revolution::project_function(msh, cl, pressure, di.face_degree());
            auto pb  = revolution::make_scalar_monomial_basis(msh, cl, di.face_degree());

            Matrix<scalar_type, Dynamic, 1> spres = assembler.take_pressure(msh, cl, sol);
        	Matrix<scalar_type, Dynamic, 1> diff_pres = spres - ppres.head(pbs);
        	Matrix<scalar_type, Dynamic, Dynamic> scalar_mm = revolution::make_mass_matrix(msh, cl, pb);
        	error_pres += diff_pres.dot(scalar_mm*diff_pres);

            //energy error
            Matrix<scalar_type, Dynamic, 1> svel =  assembler.take_velocity(msh, cl, sol);
            Matrix<scalar_type, Dynamic, 1> diff_vel = svel - p;
            auto gr = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<scalar_type, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

            Matrix<scalar_type, Dynamic, Dynamic> B = factor * (viscosity*G.second + viscosity*stab);

            error_vel += diff_vel.dot(B * diff_vel);
        }

        //ofs.close();
        return std::make_pair(std::sqrt(error_vel), std::sqrt(error_pres));
    }

    template<typename Assembler>
    void
    save_auxiliar(  const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler)
    {
        auto gamma = compute_auxiliar(msh, cl, assembler);
        auto cell_ofs  = revolution::priv::offset(msh, cl);
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;
        return;
    }

    template<typename Assembler>
    Matrix<scalar_type, Dynamic, 1>
    compute_auxiliar(   const mesh_type& msh,
                        const cell_type& cl,
                        const Assembler& assembler)
    {
        auto u_TF  = assembler.take_velocity(msh, cl, sol_old);
        auto value = 1./(factor*(viscosity + alpha));
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);

        auto cell_ofs    = revolution::priv::offset(msh, cl);
        Matrix<scalar_type, Dynamic, 1> Gu = G.first * u_TF;
        Matrix<scalar_type, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);
        Matrix<scalar_type, Dynamic, 1> gamma = value * (stress +  factor * alpha * Gu);

        return gamma;
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler)
    {
        convergence = 0.;
        auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto u_TF = assembler.take_velocity(msh, cl, sol);
            auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto cell_ofs = revolution::priv::offset(msh, cl);

            Matrix<scalar_type, Dynamic, 1> Gu = G.first * u_TF;
            Matrix<scalar_type, Dynamic, 1> gamma = compute_auxiliar( msh,  cl, assembler);
            Matrix<scalar_type, Dynamic, 1> gamma_old = auxiliar.block(cell_ofs *sbs, 0, sbs, 1);

            Matrix<scalar_type, Dynamic, 1> diff_stress = factor * alpha * (Gu - gamma);
            Matrix<scalar_type, Dynamic, 1> diff_gamma  = alpha * (gamma - gamma_old);

            auto sb =  revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());

            Matrix<scalar_type, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, sb);

            convergence += diff_stress.dot(mass * diff_stress);//+ diff_gamma.dot(mass * diff_gamma);

            multiplier.block(cell_ofs * sbs, 0, sbs, 1) +=  diff_stress;
        }

        return;
    }

    template<typename Assembler>
    Matrix<scalar_type, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = revolution::make_sym_matrix_monomial_basis(msh, cl, di.face_degree());
        auto cell_ofs =  revolution::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        auto stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        auto gamma  = compute_auxiliar( msh,  cl, assembler);
        vector_type str_agam = stress - factor * alpha * gamma;

        Matrix<scalar_type, Dynamic, Dynamic> mm = revolution::make_mass_matrix(msh, cl, sb);

        Matrix<scalar_type, Dynamic, 1> rhs =
                    Matrix<scalar_type, Dynamic, 1>::Zero(cbs + fbs * num_faces);

        //(f, v_T)
        rhs.block( 0, 0, cbs, 1) = make_rhs(msh, cl, cb, rhs_fun);

        //(stress - alpha * gamma, Gv)
        rhs -=  G.first.transpose() * mm * str_agam;

        return rhs;
    }

    template<typename Assembler>
    auto
    run_stokes_like(const mesh_type& msh, Assembler& assembler)
    {
        auto celdeg = di.cell_degree();

        sol_old = sol;

        assembler.initialize();

        for (auto cl : msh)
        {
            auto cb  = revolution::make_vector_monomial_basis(msh, cl, celdeg);
            auto G   = revolution::make_hlow_stokes(msh, cl, di, use_sym_grad);
            auto gr  = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<scalar_type, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto dr  = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            Matrix<scalar_type, Dynamic, 1> local_rhs =
                Matrix<scalar_type, Dynamic, 1>::Zero(cbs + fbs * howmany_faces(msh,cl));
            local_rhs = make_rhs_alg(msh, cl, assembler);

            Matrix<scalar_type, Dynamic, Dynamic> A = factor * (alpha * G.second + viscosity * stab);

            assembler.assemble_alg(msh, cl, A, -dr.second, local_rhs);

            save_auxiliar(msh, cl, assembler);
        }

        assembler.finalize();

        //dump_sparse_matrix(assembler.LHS, "stokes.txt");

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        sol = vector_type::Zero(systsz);

        disk::solvers::pardiso_params<scalar_type> pparams;
        mkl_pardiso_ldlt(pparams, assembler.LHS, assembler.RHS, sol);
        //mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);

        //std::cout << "sol : "<< sol.transpose() << std::endl;
        //std::cout << sol.transpose()  << std::endl;
        //std::ofstream ofs("velocity.dat");

        //compute_convergence_velocity(msh, assembler);
        return;
    }

};


template<typename Mesh, typename T, typename Assembler>
void
quiver( const Mesh& msh, const dynamic_vector<T>& sol, const Assembler& assembler,
        const typename revolution::hho_degree_info & di)
{
    std::ofstream ofs("quiver_vel.data");

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

template<typename Mesh>
auto
run_alg_stokes(const Mesh& msh, size_t degree, std::ofstream & ofs,
                bool use_sym_grad )
{
    using T = typename Mesh::coordinate_type;
    T tolerance = 1.e-9, Ninf = 10.e+5;
    size_t max_iters = 50000;

    typename revolution::hho_degree_info hdi(degree, degree);
    augmented_lagrangian_stokes<Mesh> als(msh, hdi, use_sym_grad);

    auto assembler = als.define_assembler(msh);
    als.initialize(msh, assembler);

    T error_old = 0.;

    for(size_t i = 0; i < max_iters; i++)
    {

        als.run_stokes_like(msh, assembler);
        als.update_multiplier(msh, assembler);
        auto error = als.compute_errors(msh, assembler);
        auto convergence = std::sqrt(als.convergence);

        //ofs << "  i : "<< i<<"  " << convergence<< " ------  ";
        //ofs << std::scientific << std::setprecision(4) << error.first;
        //ofs << "     -- " << "          ";
        //ofs << std::scientific << std::setprecision(4) << error.second;
        //ofs << "     -- " << std::endl;

        if(i % 100 == 0)
            std::cout << "  i : "<< i<<"  - " << convergence<<std::endl;
        assert(convergence < Ninf);
        if( convergence < tolerance)
            break;
            //#if 0
        if(std::abs(error.first - error_old)/error.first < 10.e-8 && i > 1)
        {
            std::cout << "Break by convergence of velocity error : ";
            std::cout << std::abs(error.first - error_old)/error.first << std::endl;
            break;
        }
        error_old = error.first;

    }

    auto error_2 = als.compute_errors(msh, assembler);
    quiver( msh, als.sol, assembler, hdi);
    return error_2;
}

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

void convergence_test_typ1(void)
{
    using T = double;
    bool use_sym_grad = true;
    std::vector<std::string> meshfiles;

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_5.typ1");
    //meshfiles.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_6.typ1");
    */

    //meshfiles.push_back("../../../diskpp/meshes/2D_quads/medit/square_h00125.medit2d");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_quads/fvca5/mesh2_5.typ1");

    /*
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    meshfiles.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");
    */
    std::cout << "                   velocity H1-error";
    std::cout << "    -     pressure L2-error "<< std::endl;

    //for (size_t k = 0; k < 3; k++)
    size_t k= 0;
    {
        std::cout << "DEGREE " << k << std::endl;

        std::ofstream ofs("errors_k" + tostr(k) + ".data");
        if (!ofs.is_open())
            std::cout << "Error opening errors "<<std::endl;

        std::vector<T> mesh_hs;
        std::vector<std::pair<T,T>> errors;

        //for (size_t i = 0; i < 1; i++)
        //for (size_t i = 0; i < meshfiles.size(); i++)
        size_t i = 3;
        {
            typedef disk::generic_mesh<T, 2>  mesh_type;

            std::cout << " Mesh : "<< i << std::endl;
            mesh_type msh;
            disk::fvca5_mesh_loader<T, 2> loader;
            if (!loader.read_mesh(meshfiles.at(i)))
            {
                std::cout << "Problem loading mesh." << std::endl;
                //continue;
            }
            loader.populate_mesh(msh);

            auto error = run_alg_stokes(msh, k, ofs, use_sym_grad);

            mesh_hs.push_back( disk::mesh_h(msh) );
            errors.push_back(error);
            ofs << " " << std::endl;
        }

        ofs.close();

        for (size_t i = 0; i < mesh_hs.size(); i++)
        {
            if (i == 0)
            {
                std::cout << "    ";
                std::cout << std::scientific << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).first;
                std::cout << "     -- " << "          ";
                std::cout << std::scientific << std::setprecision(4) << errors.at(i).second;
                std::cout << "     -- " << std::endl;
            }
            else
            {
                auto rate = std::log( errors.at(i).first/errors.at(i-1).first ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << "    ";
                std::cout << std::scientific  << std::setprecision(4) << mesh_hs.at(i) << "    ";
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).first << "    ";
                std::cout << std::fixed<< std::setprecision(2) << rate << "          ";

                auto pres_rate = std::log( errors.at(i).second/errors.at(i-1).second ) /
                            std::log( mesh_hs.at(i)/mesh_hs.at(i-1) );
                std::cout << std::scientific  << std::setprecision(4) << errors.at(i).second << "    ";
                std::cout << std::fixed << std::setprecision(2) << pres_rate << std::endl;
            }
        }
    }
}

int main(void)
{
    convergence_test_typ1();
}

#if 0

int main(int argc, char **argv)
{
    using RealType = double;

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

        run_stokes(msh, degree);
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

        run_stokes(msh, degree);
    }

/*
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        test_bases(msh);
    }
    */

    return 0;
}
#endif
