
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

template<typename Mesh>
void
save_coords(const Mesh& msh,
            const std::string & filename)
{
    std::ofstream ofs("Coords.data" );
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
                                point<typename Mesh::scalar_type,2>>   & e,
                const Matrix<typename Mesh::scalar_type, Dynamic, 1> & sol,
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

                vector_type s = sol.block(cell_ofs * cbs, 0, cbs, 1);
                auto cb = revolution::make_vector_monomial_basis(msh, cl, cell_degree);
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

template< typename T>
std::string
tostr(const T & var)
{
    std::ostringstream  ostr;
    ostr << var;
    return ostr.str();
}

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
    size_t                  problem_num;
    bool                    use_sym_grad;

    augmented_lagrangian_viscoplasticity(const Mesh& msh,
                            const typename revolution::hho_degree_info & hdi,
                            const T& alpha_ext):
                            di(hdi), alpha(alpha_ext)
    {
        use_sym_grad = false;
        factor = 1.;//(use_sym_grad)? 2. : 1.;
        viscosity = 1;

        T f = 1;
        T Lref  = 1.;
        T Bn = 0.;
        yield =  Bn;// * f * Lref/ viscosity;

        dim =  Mesh::dimension;

        cbs = revolution::vector_basis_size(di.cell_degree(), dim, dim);
        fbs = revolution::vector_basis_size(di.face_degree(), dim - 1, dim);
        pbs = revolution::scalar_basis_size(di.face_degree(), dim);
        sbs = revolution::matrix_basis_size(di.face_degree(), dim, dim);
    };

    auto
    define_problem(const mesh_type& msh)
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

        problem_num = 1;

        //switch (problem_num)
		//{
        //    case 1: //driven
                auto wall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };
                auto movingWall = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>{1,0};
                };

                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    if( std::abs(p.y() - 1.) < 1.e-8 )
                        return Matrix<T, Mesh::dimension, 1>{1,0};
                    else
                        return Matrix<T, Mesh::dimension, 1>{0,0};
                };
                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                //#if 0
                bnd.addDirichletBC(0, 1, movingWall);
                bnd.addDirichletBC(0, 2, wall);
                bnd.addDirichletBC(0, 3, wall);
                bnd.addDirichletBC(0, 4, wall);
                //#endif
                //bnd.addDirichletEverywhere(velocity);

            //    break;

                #if 0
            case 2: //Taylor-Couette

                rhs_fun  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    return Matrix<T, Mesh::dimension, 1>::Zero();
                };

                velocity  = [](const point_type& p) -> Matrix<T, Mesh::dimension, 1> {
                    Matrix<T, Mesh::dimension, 1> ret = Matrix<T, Mesh::dimension, 1>::Zero();

                    auto theta  = std::atan2(p.y() , p.x());

                    T Rint(0.5), Rext(1.), omegaInt(2.),omegaExt(2.);
                    T r = std::sqrt(p.x()*p.x() + p.y()*p.y());

                    //1.All Solid
                    auto u_theta_solid = (omegaExt/ omegaInt) * (r / Rint);
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

                bnd.addNeumannBC(10, 1, symmetryPlane);
                bnd.addNeumannBC(10, 2, symmetryPlane);
                bnd.addDirichletBC(0, 3, velocity);
                bnd.addDirichletBC(0, 4, velocity);

                break;
                #endif
        //}

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
            auto G = make_hlow_vector_laplacian(msh, cl, di);

            Matrix<T, Dynamic, Dynamic> B = (viscosity*G.second + viscosity*stab);

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
        auto u_TF  = assembler.take_velocity(msh, cl, velocity_dofs);
        auto value = 1./(viscosity + alpha);
        auto G = make_hlow_vector_laplacian(msh, cl, di);

        auto cell_ofs    = revolution::priv::offset(msh, cl);
        Matrix<T, Dynamic, 1> Gu = G.first * u_TF;
        Matrix<T, Dynamic, 1> stress = multiplier.block(cell_ofs * sbs, 0, sbs, 1);

        //Theta
        auto sb =  revolution::make_matrix_monomial_basis(msh, cl, di.face_degree());
        //barycenter only for k = 0; fix this for higher orders
        auto bar = barycenter(msh, cl);
        auto s_phi  = sb.eval_functions(bar);

        Matrix<T, Dynamic, 1>  theta  = stress  +  alpha * Gu;
        Matrix<T, Mesh::dimension, Mesh::dimension> theta_eval = revolution::eval(theta, s_phi);
        T theta_norm  = theta_eval.norm();

        //Gamma
        // A. Solid
        Matrix<T, Dynamic, 1> gamma(sbs);
        // B. Liquid
        if(theta_norm > yield )
            gamma += value * (theta/theta_norm) * (theta_norm - yield);

        return gamma;
    }

    template<typename Assembler>
    auto
    update_multiplier(const mesh_type& msh, const Assembler& assembler)
    {
        T conv_total = 0.;
        T conv_stress = 0.;
        T conv_gamma = 0.;

        auto dim = Mesh::dimension;

        for(auto cl: msh)
        {
            auto u_TF = assembler.take_velocity(msh, cl, sol);
            auto G  = make_hlow_vector_laplacian(msh, cl, di);
            auto cell_ofs = revolution::priv::offset(msh, cl);

            Matrix<T, Dynamic, 1> Gu = G.first * u_TF;
            Matrix<T, Dynamic, 1> gamma_old = auxiliar_old.block(cell_ofs *sbs, 0, sbs, 1);
            Matrix<T, Dynamic, 1> gamma = auxiliar.block(cell_ofs *sbs, 0, sbs, 1);

            Matrix<T, Dynamic, 1> diff_stress = alpha * (Gu - gamma);
            Matrix<T, Dynamic, 1> diff_gamma  = alpha * (gamma - gamma_old);

            auto sb =  revolution::make_matrix_monomial_basis(msh, cl, di.face_degree());
            Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, sb);

            conv_stress += diff_stress.dot(mass * diff_stress);
            conv_gamma  += diff_gamma.dot(mass * diff_gamma);
            conv_total = conv_stress + conv_gamma;
            convergence = std::make_tuple(conv_total, conv_stress, conv_gamma);
            multiplier.block(cell_ofs * sbs, 0, sbs, 1) += alpha * diff_stress;
        }

        auxiliar_old = auxiliar;
        return;
    }

    template<typename Assembler>
    Matrix<T, Dynamic, 1>
    make_rhs_alg(   const mesh_type& msh,
                    const cell_type& cl,
                    const Assembler& assembler)
    {
        auto G = make_hlow_vector_laplacian(msh, cl, di);
        auto cb = revolution::make_vector_monomial_basis(msh, cl, di.cell_degree());
        auto sb = revolution::make_matrix_monomial_basis(msh, cl, di.face_degree());
        auto cell_ofs =  revolution::priv::offset(msh, cl);
        auto num_faces = howmany_faces(msh, cl);

        auto stress = multiplier.block( sbs * cell_ofs,  0, sbs, 1);
        auto gamma  = compute_auxiliar( msh,  cl, assembler, sol_old); //or sol at this point it's the same
        auxiliar.block(cell_ofs * sbs, 0, sbs, 1) = gamma;

        vector_type str_agam = stress - alpha * gamma;

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
        auto celdeg = di.cell_degree();

        assembler.initialize_lhs();

        for (auto cl : msh)
        {
            auto cb  = revolution::make_vector_monomial_basis(msh, cl, celdeg);
            auto G   = make_hlow_vector_laplacian(msh, cl, di);
            auto gr  = revolution::make_hho_stokes(msh, cl, di, use_sym_grad);
            Matrix<T, Dynamic, Dynamic> stab;
            stab = make_hho_fancy_stabilization_vector(msh, cl, gr.first, di);
            auto dr  = make_hho_divergence_reconstruction_stokes_rhs(msh, cl, di);

            Matrix<T, Dynamic, Dynamic> A = (alpha * G.second + viscosity * stab);

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

};
