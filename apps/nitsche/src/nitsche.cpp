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
#include "mumps.hpp"
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
auto hho_nitsche_reconstruction(const Mesh& msh,
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

    /* Stiffness */
    disk::dynamic_matrix<scalar_type> K =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);

    /* Nitsche contributions (consistency, symmetry, penalization) */
    disk::dynamic_matrix<scalar_type> Nf =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);
    
    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs + n_allfacedofs);
    
    auto qps = disk::integrate(msh, cl, 2*degree);
    for (const auto& qp : qps) {
        /* (grad(v), grad(w))_T */ 
        auto dphi = rb.eval_gradients(qp.point());
        K += (qp.weight() * dphi) * dphi.transpose();
    }

    auto inv_hT = 1.0/diameter(msh, cl);
    for (size_t fcnum = 0; fcnum < fcs.size(); fcnum++) {
        const auto& fc = fcs[fcnum];
        auto fb = disk::make_scalar_monomial_basis(msh, fc, degree);
        auto bi = msh.boundary_info(fc);
        auto ofs = rbs + fbs*fcnum;
        auto fqps = disk::integrate(msh, fc, 2*degree+2);
        auto n = normal(msh, cl, fc);
        auto fcid = offset(msh, fc);

        if (bi.is_boundary()) { /* Do Nitsche if on a domain boundary */

            /* See chapter 4 of "Mathematical aspects of DG methods"
             * by Di Pietro & Ern, in particular (4.12) and (4.16). */
            if (bcs[fcid] == bc::dirichlet) {
                /* Dirichlet needs all the three "SIP-DG" boundary terms*/
                for (const auto& qp : fqps) {
                    auto phi = rb.eval_functions(qp.point());
                    auto dphi = rb.eval_gradients(qp.point());
                    auto dphi_dot_n = dphi*n;
                    /* (v, grad(w)).n)_F */
                    Nf -= (qp.weight() * dphi_dot_n) * phi.transpose();
                    /* (grad(v).n, w)_F */
                    Nf -= (qp.weight() * phi) * dphi_dot_n.transpose();
                    /* (eta/hT) * (v,w)_F */
                    Nf += (inv_hT*eta*qp.weight() * phi) * phi.transpose();
                }
            }

            if (bcs[fcid] == bc::neumann) {
                /* Only stiffness term */
            }

        } else { /* Do standard HHO if not on a domain boundary */
            for (const auto& qp : fqps) {
                auto cphi = rb.eval_functions(qp.point());
                auto fphi = fb.eval_functions(qp.point());
                auto dphi = rb.eval_gradients(qp.point());
                auto dphi_dot_n = dphi*n;
                RHS.block(0,   0, rbs, rbs) -= qp.weight() * dphi_dot_n * cphi.transpose();
                RHS.block(0, ofs, rbs, fbs) += qp.weight() * dphi_dot_n * fphi.transpose();
            }
        }
    }

    disk::dynamic_matrix<scalar_type> Nt = K + Nf;
    RHS.block(0, 0, rbs, rbs) += Nt;

    /* We didn't do any kind of average fixing: we rely on LDLT to
     * have one. This point must be clarified a bit to figure out
     * _which_ one. */
    disk::dynamic_matrix<scalar_type> oper = Nt.ldlt().solve(RHS);
    disk::dynamic_matrix<scalar_type> data = oper.transpose() * Nt * oper;

    return std::pair{oper, data};
}

template<typename Mesh>
disk::dynamic_matrix<typename Mesh::coordinate_type>
hho_nitsche_stabilization(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree)
{
    /* Nitsche-HHO as implemented here is mixed-order (k+1 on cells
     * and k on faces). We use a standard Lehrenfeld-Schoeberl
     * stabilization. We need to stabilize only on the internal
     * interfaces, not on the domain boundary. */
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
        size_t offset = cbs+i*fbs;
        const auto fc = fcs[i];

        /* If the face is on the domain boundary, just skip to the next */
        auto bi = msh.boundary_info(fc);
        if (bi.is_boundary()) {
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

        oper.block(0, offset, fbs, fbs) = -If;

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

        tr.block(0, offset, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs)      = trace;

        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);
        data += oper.transpose() * tr * stabparam;
    }

    return data;
}

template<typename Mesh, typename SourceFun, typename DirichetFun, typename NeumannFun>
disk::dynamic_vector<typename Mesh::coordinate_type>
hho_nitsche_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
    SourceFun f, DirichetFun gD, NeumannFun gN, size_t degree,
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

        /* Compute the Dirichlet contributions on the RHS.
         * Same stuff as SIP-DG */
        if (bcs[fcid] == bc::dirichlet) {
            for (const auto& qp : fqps) {
                auto phi = cb.eval_functions(qp.point());
                auto dphi = cb.eval_gradients(qp.point());
                auto dphi_dot_n = dphi*n;
                auto gD_val = gD(qp.point());
                /* (gD, grad(w).n)_F */
                ret.head(cbs) -= (qp.weight() * gD_val) * dphi_dot_n;
                /* (eta/hT) * (gD, w)_F */
                ret.head(cbs) += (inv_hT*eta*qp.weight() * gD_val) * phi.transpose();
            }
        }

        /* Compute the Neumann contributions on the RHS.
         * Same stuff as SIP-DG */
        if (bcs[fcid] == bc::neumann) {
            for (const auto& qp : fqps) {
                auto phi = cb.eval_functions(qp.point());
                auto gN_val = gN(qp.point());
                /* (gN, w)_F */
                ret.head(cbs) += gN_val * phi.transpose();
            }
        }
    }

    return ret;
}

template<typename Mesh>
auto
nitsche_hho_solver(const Mesh& msh, size_t degree, const std::vector<bc>& bcs)
{
    using scalar_type = typename Mesh::coordinate_type;

    scalar_type eta = 1.0;

    auto zerofun = [](const typename Mesh::point_type&) {
        return 0.0;
    };

    auto onefun = [](const typename Mesh::point_type&) {
        return 1.0;
    };

    disk::hho::slapl::degree_info di(degree+1, degree);
    auto assm = make_assembler(msh, di); // XXX

    timecounter tc;
    tc.tic();

    using MT = disk::dynamic_matrix<scalar_type>;
    using VT = disk::dynamic_vector<scalar_type>;
    std::vector<std::pair<MT, VT>> lcs;

    auto rhsfun = make_rhs_function(msh);

    for (auto& cl : msh) {
        auto [R, A] = hho_nitsche_reconstruction(msh, cl, degree, eta, bcs);
        auto S = hho_nitsche_stabilization(msh, cl, degree);
        disk::dynamic_matrix<scalar_type> lhs = A+S;

        disk::dynamic_vector<scalar_type> rhs =
            hho_nitsche_rhs(msh, cl,
                rhsfun,     // source
                zerofun,    // dirichlet
                zerofun,    // neumann
                degree, eta, bcs);

        lcs.push_back({lhs, rhs});
        
        auto cbs = disk::scalar_basis_size(degree+1, Mesh::dimension);
        auto [Lc, Rc] = disk::static_condensation(lhs, rhs, cbs);
    
        assm.assemble(msh, cl, Lc, Rc);
    }
    assm.finalize();

    //std::cout << " Assembly time: " << tc.toc() << std::endl;
    //std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    //std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();
    disk::dynamic_vector<scalar_type> sol = mumps_lu(assm.LHS, assm.RHS);
    //std::cout << " Solver time: " << tc.toc() << std::endl;
    
    std::vector<scalar_type> u_data;
    
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

        disk::dynamic_vector<scalar_type> ana_sol =
            disk::project_function(msh, cl, degree+1, solfun);

        disk::dynamic_vector<scalar_type> diff = ana_sol - locsol.head(cbs);

        auto cb = disk::make_scalar_monomial_basis(msh, cl, degree+1);
        disk::dynamic_matrix<scalar_type> mass = disk::make_mass_matrix(msh, cl, cb);

        L2error += diff.dot(mass*diff);
    }
    //std::cout << " Postpro time: " << tc.toc() << std::endl;
    //std::cout << " L2-norm error: " << std::sqrt(L2error) << std::endl;

    disk::silo_database silo;
    silo.create("nitsche.silo");
    silo.add_mesh(msh, "mesh");
    silo.add_variable("mesh", "u", u_data, disk::zonal_variable_t);

    return std::sqrt(L2error);
}

int main(void)
{
    using T = double;
    using mesh_type = disk::cartesian_mesh<T,2>;


    for (size_t k = 0; k < 5; k++) {
        mesh_type msh;
        auto mesher = make_simple_mesher(msh);
        
        auto prev_err = 0.0;
        auto prev_h = 0.0;

        std::cout << "Nitsche-HHO(k+1, k), k = " << k << std::endl;
        for (size_t i = 0; i < 4; i++) {
            mesher.refine();
            std::vector<bc> bcs;
            set_boundary(msh, bcs, bc::neumann, 0);
            set_boundary(msh, bcs, bc::dirichlet, 1);
            set_boundary(msh, bcs, bc::neumann, 2);
            set_boundary(msh, bcs, bc::dirichlet, 3);
            auto err = nitsche_hho_solver(msh, k, bcs);
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