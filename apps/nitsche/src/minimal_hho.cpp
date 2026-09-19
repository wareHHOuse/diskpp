/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2025
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho_assemblers.hpp"
#include "diskpp/solvers/direct_solvers.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/output/silo.hpp"
#include "asm.hpp"

/***************************************************************
 * Sources for the test problems
 */
template<typename Mesh>
struct source_functor;

template<disk::mesh_2D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        //return 2.0*M_PI*M_PI*sx*sy;
        return M_PI*M_PI*sx;
    }
};

template<disk::mesh_3D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return 3.0*M_PI*M_PI*sx*sy*sz;
    }
};


/***************************************************************
 * Analytical solutions for the test problems
 */
template<typename Mesh>
struct solution_functor;

template<disk::mesh_2D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        //return sx*sy;
        return sx;
    }
};

template<disk::mesh_3D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return sx*sy*sz;
    }
};

/***************************************************************
 * Helpers
 */
template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return source_functor<Mesh>();
}

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}



/***************************************************************
 * Nitsche-HHO reconstruction operator
 */
template<typename Mesh>
auto hho_minimal_reconstruction(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree,
    typename Mesh::coordinate_type eta, const std::vector<bc>& bcs)
{
    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
    auto rbs = rb.size();
    /* Cell space basis: same as reconstruction space */ 
    //auto cb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
    //auto cbs = cb.size();
    /* Face basis info */
    auto fcs = faces(msh, cl);
    auto fbs = disk::scalar_basis_size(degree, Mesh::dimension-1);
    auto n_allfacedofs = fcs.size() * fbs;

    disk::dynamic_matrix<scalar_type> LHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs-1, rbs-1);

    /* Stiffness */
    disk::dynamic_matrix<scalar_type> K =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);
    
    /* Robin */
    disk::dynamic_matrix<scalar_type> R =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);

    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs-1, rbs + n_allfacedofs);
    

    auto qps = disk::integrate(msh, cl, 2*degree);
    for (const auto& qp : qps) {
        /* (grad(v), grad(w))_T */ 
        auto dphi = rb.eval_gradients(qp.point());
        K += (qp.weight() * dphi) * dphi.transpose();
    }

    LHS.block(0,0,rbs-1,rbs-1) = K.block(1,1,rbs-1,rbs-1);
    RHS.block(0,0,rbs-1,rbs) = K.block(1,0,rbs-1,rbs);

    auto inv_hT = 1.0/diameter(msh, cl);
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        const auto& fc = fcs[fcnum];
        auto fb = disk::make_scalar_monomial_basis(msh, fc, degree);
        auto bi = msh.boundary_info(fc);
        auto ofs = rbs + fbs*fcnum;
        auto fqps = disk::integrate(msh, fc, 2*degree+2);
        auto n = normal(msh, cl, fc);
        auto fcid = offset(msh, fc);

        if (bi.is_boundary()) { /* Do "minimal hho" if on a domain boundary */

            if (bcs[fcid] == bc::dirichlet) {
                for (const auto& qp : fqps) {
                    auto cphi = rb.eval_functions(qp.point());
                    auto fphi = fb.eval_functions(qp.point());
                    auto dphi = rb.eval_gradients(qp.point());
                    disk::dynamic_vector<scalar_type> dphi_dot_n =
                        (dphi*n).tail(rbs-1);
                    RHS.block(0,   0, rbs-1, rbs) -= qp.weight() * dphi_dot_n * cphi.transpose();
                    RHS.block(0, ofs, rbs-1, fbs) += qp.weight() * dphi_dot_n * fphi.transpose();
                }
            }

            if (bcs[fcid] == bc::neumann) {

            }

            if (bcs[fcid] == bc::robin) {
                
                for (const auto& qp : fqps) {
                    auto cphi = rb.eval_functions(qp.point());
                    R += qp.weight() * cphi * cphi.transpose();
                }

            }

        } else { /* Do standard HHO if not on a domain boundary */
            for (const auto& qp : fqps) {
                auto cphi = rb.eval_functions(qp.point());
                auto fphi = fb.eval_functions(qp.point());
                auto dphi = rb.eval_gradients(qp.point());
                disk::dynamic_vector<scalar_type> dphi_dot_n =
                    (dphi*n).tail(rbs-1);
                RHS.block(0,   0, rbs-1, rbs) -= qp.weight() * dphi_dot_n * cphi.transpose();
                RHS.block(0, ofs, rbs-1, fbs) += qp.weight() * dphi_dot_n * fphi.transpose();
            }
        }
    }

    disk::dynamic_matrix<scalar_type> oper = LHS.ldlt().solve(RHS);
    disk::dynamic_matrix<scalar_type> data = oper.transpose() * RHS;

    data.block(0,0,rbs,rbs) += R;

    return std::pair{oper, data};
}

template<typename Mesh>
disk::dynamic_matrix<typename Mesh::coordinate_type>
hho_minimal_stabilization(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree, const std::vector<bc>& bcs)
{
    /* We use a standard Lehrenfeld-Schoeberl stabilization and we
     * need to stabilize only on the internal interfaces, not on
     * the domain boundary. */
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = degree+1;
    const auto cb = disk::make_scalar_monomial_basis(msh, cl, celdeg);
    const auto cbs = cb.size();

    const auto fcs = faces(msh, cl);
    const auto fbs = disk::scalar_basis_size(degree, Mesh::dimension-1);
    const auto num_faces_dofs = fbs*fcs.size();
    const auto total_dofs     = cbs + num_faces_dofs;

    matrix_type data = matrix_type::Zero(total_dofs, total_dofs);

    T hT = diameter(msh, cl);
    T stabparam = 1.0/hT;

    for (size_t i = 0; i < fcs.size(); i++) {
        size_t ofs = cbs+i*fbs;
        const auto fc = fcs[i];

        /* If the face is on the domain boundary, just skip to the next */
        auto bi = msh.boundary_info(fc);
        auto fcid = offset(msh, fc);
        if (bi.is_boundary() and (bcs[fcid] != bc::dirichlet)) {
            continue;
        }

        /* Compute standard L-S stabilization otherwise. */
        const auto facdeg = degree;
        const auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto fbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);

        const matrix_type If    = matrix_type::Identity(fbs, fbs);
        matrix_type       oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type       tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type       mass  = make_mass_matrix(msh, fc, fb);
        matrix_type       trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, ofs, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            trace += (qp.weight() * f_phi) * c_phi.transpose();
        }

        tr.block(0, ofs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs)      = trace;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * stabparam;
    }

    return data;
}

template<typename Mesh, typename SourceFun,
    typename NeumannFun, typename RobinFun>
disk::dynamic_vector<typename Mesh::coordinate_type>
hho_minimal_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
    SourceFun f, NeumannFun gN, RobinFun gR, size_t degree,
    typename Mesh::coordinate_type eta, const std::vector<bc>& bcs)
{
    using scalar_type = typename Mesh::coordinate_type;

    auto cb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
    auto cbs = cb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::scalar_basis_size(degree, Mesh::dimension-1);
    auto n_allfacedofs = fcs.size() * fbs;

    disk::dynamic_vector<scalar_type> ret =
        disk::dynamic_vector<scalar_type>::Zero(cbs + n_allfacedofs);
    
    auto qps = disk::integrate(msh, cl, 2*degree+2);
    for (auto& qp : qps) {
        auto phi = cb.eval_functions(qp.point());
        ret.head(cbs) += qp.weight() * f(qp.point()) * phi;
    }

    auto inv_hT = 1.0/diameter(msh, cl);
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        const auto& fc = fcs[fcnum];
        auto bi = msh.boundary_info(fc);
        if ( not bi.is_boundary() ) {
            continue;
        }

        auto fb = disk::make_scalar_monomial_basis(msh, fc, degree);
        auto ofs = cbs + fbs*fcnum;
        auto fqps = disk::integrate(msh, fc, 2*degree+2);
        auto n = normal(msh, cl, fc);

        auto fcid = offset(msh, fc);

        if (bcs[fcid] == bc::dirichlet) {
            /* Nothing */
        }

        if (bcs[fcid] == bc::neumann) {
            for (const auto& qp : fqps) {
                auto phi = cb.eval_functions(qp.point());
                auto gN_val = gN(qp.point());
                /* (gN, w)_F */
                ret.head(cbs) += gN_val * phi.transpose();
            }
        }

        if (bcs[fcid] == bc::robin) {
            for (const auto& qp : fqps) {
                auto phi = cb.eval_functions(qp.point());
                auto gR_val = gR(qp.point());
                /* (gN, w)_F */
                ret.head(cbs) += gR_val * phi.transpose();
            }
        }
    }

    return ret;
}

template<typename Mesh>
auto
minimal_hho_solver(const Mesh& msh, size_t degree, const std::vector<bc>& bcs)
{
    bool compute_cond = true;
    using scalar_type = typename Mesh::coordinate_type;

    scalar_type eta = 1.0;

    auto zerofun = [](const typename Mesh::point_type&) {
        return 0.0;
    };

    auto onefun = [](const typename Mesh::point_type&) {
        return 1.0;
    };

    auto u = [](const typename Mesh::point_type& pt) {
        auto x = pt.x();
        auto y = pt.y();
        //return std::exp(x) * std::sin(M_PI*y) + x*x + y;
        return x*(1-x)*y;
        //return std::cos(M_PI*x);
    };

    auto f = [](const typename Mesh::point_type& pt) {
        auto x = pt.x();
        auto y = pt.y();
        //return (M_PI*M_PI - 1) * std::exp(x)*std::sin(M_PI*y) - 2;
        return 2*y;
        //return M_PI*M_PI*std::cos(M_PI*x);
    };    

    auto gD = [&](const typename Mesh::point_type& pt) {
        return u(pt);
    };

    auto gN = [](const typename Mesh::point_type& pt) {
        auto x = pt.x();
        auto y = pt.y();
        //return M_E * std::sin(M_PI*y) + 2;
        return x*(1-x);
    };

    auto gR = [](const typename Mesh::point_type& pt) {
        auto x = pt.x();
        auto y = pt.y();
        auto alpha = 1.0;
        return 0;//(-M_PI*std::exp(x) - 1.0 + alpha*x*x);
    };


    disk::hho::slapl::degree_info di(degree+1, degree);

    const static size_t DIM = Mesh::dimension;
    auto cbs = disk::scalar_basis_size(degree+1, DIM);
    auto fbs = disk::scalar_basis_size(degree, DIM-1);

    std::vector<bool> dirfaces;
    dirfaces.resize( bcs.size() );

    auto df = [&](bc b) {
        return (b == bc::dirichlet) or (b == bc::neumann) or (b == bc::robin);
    };

    std::transform(bcs.begin(), bcs.end(), dirfaces.begin(), df);

    condensed_assembler assm(msh, fbs, dirfaces);


    timecounter tc;
    tc.tic();

    using MT = disk::dynamic_matrix<scalar_type>;
    using VT = disk::dynamic_vector<scalar_type>;
    std::vector<std::pair<MT, VT>> lcs;

    auto rhsfun = make_rhs_function(msh);

    for (auto& cl : msh) {
        auto [R, A] = hho_minimal_reconstruction(msh, cl, degree, eta, bcs);
        auto S = hho_minimal_stabilization(msh, cl, degree, bcs);
        disk::dynamic_matrix<scalar_type> lhs = A+S;

        disk::dynamic_vector<scalar_type> rhs =
            hho_minimal_rhs(msh, cl,
                f,     // source
                gN,    // neumann
                gR,    // robin
                degree, eta, bcs);


        disk::dynamic_vector<scalar_type> gD_rhs =
            disk::dynamic_vector<scalar_type>::Zero(A.rows());
        auto fcs = faces(msh, cl);
        auto ofs = cbs;
        for (auto& fc : fcs) {
            auto bi = msh.boundary_info(fc);
            if (bi.is_boundary()){
                auto boundary_id = bi.id();
                if ( bcs[offset(msh, fc)] == bc::dirichlet ) {
                    auto fb = disk::make_scalar_monomial_basis(msh, cl, degree);
                    auto fqps = disk::integrate(msh, fc, 2*degree);
                    disk::dynamic_matrix<scalar_type> M =
                        disk::dynamic_matrix<scalar_type>::Zero(fb.size(), fb.size());
                    disk::dynamic_vector<scalar_type> f_gD =
                        disk::dynamic_vector<scalar_type>::Zero(fb.size());
                    for (const auto& qp : fqps) {
                        auto phi = fb.eval_functions(qp.point());
                        M += qp.weight() * phi * phi.transpose();
                        f_gD += qp.weight() * gD(qp.point()) * phi;
                    }
                    gD_rhs.segment(ofs, fb.size()) += M.ldlt().solve(f_gD);
                }
            }
            ofs += fbs;
        }

        rhs += -lhs*gD_rhs;

        lcs.push_back({lhs, rhs});
        
        auto cbs = disk::scalar_basis_size(degree+1, Mesh::dimension);
        auto [Lc, Rc] = disk::static_condensation(lhs, rhs, cbs);
    
        assm.assemble(msh, cl, Lc, Rc);
    }
    assm.finalize();

    std::cout << "************" << std::endl;
    //std::cout << " Assembly time: " << tc.toc() << std::endl;
    auto bfsize = msh.faces_size() - msh.boundary_faces_size();
    std::cout << " Internal faces:    " << bfsize << ", fbs = " << fbs;
    std::cout << ", intfaces*fbs = " << bfsize * fbs << std::endl;
    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();
    std::cout << "  Solver: " << std::flush;
    disk::dynamic_vector<scalar_type> sol;
    disk::solvers::sparse_lu(assm.LHS, assm.RHS, sol);
    //std::cout << " Solver time: " << tc.toc() << std::endl;
    
    std::vector<scalar_type> u_data;
    std::vector<scalar_type> uex_data;
    std::vector<scalar_type> conditioning;
    auto solfun = make_solution_function(msh);

    scalar_type L2error = 0.0;
    auto u_sol = make_solution_function(msh);
    tc.tic();
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        const auto& [lhs, rhs] = lcs[cell_i++];
        auto locsolF = assm.take_local_solution(msh, cl, sol);
        auto cbs = disk::scalar_basis_size(degree+1, Mesh::dimension);
        disk::dynamic_vector<scalar_type> locsol =
            disk::static_decondensation(lhs, rhs, locsolF);
        u_data.push_back(locsol(0));
        uex_data.push_back( u(barycenter(msh,cl)) );

        disk::dynamic_vector<scalar_type> ana_sol =
            disk::project_function(msh, cl, degree+1, u);

        disk::dynamic_vector<scalar_type> diff = ana_sol - locsol.head(cbs);

        auto cb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
        disk::dynamic_matrix<scalar_type> mass = disk::make_mass_matrix(msh, cl, cb);

        if (compute_cond) {
            conditioning.push_back( cond(lhs) );
        }

        L2error += diff.dot(mass*diff);
    }
    //std::cout << " Postpro time: " << tc.toc() << std::endl;
    //std::cout << " L2-norm error: " << std::sqrt(L2error) << std::endl;

    disk::silo_database silo;
    silo.create("nitsche.silo");
    silo.add_mesh(msh, "mesh");
    silo.add_variable("mesh", "u", u_data, disk::zonal_variable_t);
    silo.add_variable("mesh", "u_ex", uex_data, disk::zonal_variable_t);
    if (compute_cond) {
        silo.add_variable("mesh", "cond", conditioning, disk::zonal_variable_t);
    }

    return std::sqrt(L2error);
}

int main(void)
{
    using T = double;
    using mesh_type = disk::cartesian_mesh<T,2>;


    for (size_t k = 1; k < 2; k++) {
        mesh_type msh;
        auto mesher = make_simple_mesher(msh);
        
        auto prev_err = 0.0;
        auto prev_h = 0.0;

        std::cout << "Minimal-HHO(k+1, k), k = " << k << std::endl;
        for (size_t i = 0; i < 4; i++) {
            mesher.refine();
            std::vector<bc> bcs;
            set_boundary(msh, bcs, bc::dirichlet, 0);
            set_boundary(msh, bcs, bc::dirichlet, 1);
            set_boundary(msh, bcs, bc::neumann, 2);
            set_boundary(msh, bcs, bc::dirichlet, 3);
            auto err = minimal_hho_solver(msh, k, bcs);
            auto h = disk::average_diameter(msh);

            if (i == 0) {
                std::cout << "  h = " << h << ", err = " << err << std::endl;
            }
            else {
                auto rate = std::log(prev_err/err)/std::log(prev_h/h);
                std::cout << "  h = " << h << ", err = " << err << ", rate = " << rate << std::endl;
            }
            prev_h = h;
            prev_err = err;
        }
    }

    return 0;
}