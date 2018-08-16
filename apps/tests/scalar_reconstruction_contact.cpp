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

#include <xmmintrin.h>
//#define EIGEN_USE_MKL_ALL
#include "revolution/bases"
#include "revolution/quadratures"
#include "revolution/methods/hho"

#include "core/loaders/loader.hpp"

#include "common.hpp"


template<typename Mesh>
struct test_functor
{
    /* Expect k+1 convergence */
    typename Mesh::scalar_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::scalar_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = make_scalar_testing_data(msh);

        typename revolution::hho_degree_info hdi(degree);

        scalar_type error = 0.0;
        for (auto& cl : msh)
        {
            Matrix<scalar_type, Dynamic, 1> proj = revolution::project_function(msh, cl, hdi, f);
            auto gr = revolution::make_hho_scalar_laplacian(msh, cl, hdi);

            size_t rec_size = revolution::scalar_basis_size(hdi.reconstruction_degree(), Mesh::dimension);

            Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
            reconstr.tail(rec_size-1) = gr.first * proj;
            reconstr(0) = proj(0);

            auto cb = revolution::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);

            auto quad_degree = 2 * std::max(hdi.cell_degree(), hdi.reconstruction_degree());
            auto qps = revolution::integrate(msh, cl, quad_degree);

            for(auto& qp : qps)
            {
                auto c_phi = cb.eval_functions(qp.point());
                auto rec_f = reconstr.dot(c_phi);
                auto diff  = (rec_f - f(qp.point()));
                error += qp.weight() * diff * diff;
            }
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};


template<typename Mesh>
auto
make_is_dirichlet_vector(const Mesh& msh,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    //cells with contact faces
    auto num_cells = msh.cells_size();
    std::vector<size_t> ret = std::vector<size_t>(num_cells);
    size_t i =0;
    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        for (auto& fc:fcs)
        {
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

            if (bnd.is_dirichlet_face(face_id))
            {
                ret.at(i) = 1;
                continue;
            }
        }
        i++;
    }
    return ret;
}


template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_nitshce_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const revolution::hho_degree_info& di,
                          const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    const size_t DIM = Mesh::dimension;

    using vector_type   = Matrix<T, Dynamic, 1>;
    using matrix_type   = Matrix<T, Dynamic, Dynamic>;
    using gradient_type = Matrix<T, Dynamic, DIM>;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    auto cb = revolution::make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = revolution::scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = revolution::scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = revolution::scalar_basis_size(facdeg, Mesh::dimension-1);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = revolution::integrate(msh, cl, 2 * (recdeg-1));
    for (auto& qp : qps)
    {
        const auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    const auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];

        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_dirichlet_face(face_id))
            continue;

        const auto n  = normal(msh, cl, fc);
        auto fb = revolution::make_scalar_monomial_basis(msh, fc, facdeg);

        size_t quad_degree = std::max(recdeg - 1 + std::max(facdeg,celdeg), size_t(0));
        auto qps_f = revolution::integrate(msh, fc, quad_degree);
        for (auto& qp : qps_f)
        {
            vector_type     c_phi_tmp = cb.eval_functions(qp.point());
            vector_type     c_phi     = c_phi_tmp.head(cbs);
            gradient_type   c_dphi_tmp= cb.eval_gradients(qp.point());
            gradient_type   c_dphi    = c_dphi_tmp.block(1, 0, rbs-1, DIM);
            vector_type     f_phi     = fb.eval_functions(qp.point());

            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
        }
    }

    matrix_type oper = gr_lhs.llt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


template<typename Mesh>
struct test_functor_contact
{
    /* Expect k+1 convergence */
    typename Mesh::scalar_type
    operator()(const Mesh& msh, size_t degree) const
    {
        typedef Mesh mesh_type;
        typedef typename mesh_type::cell        cell_type;
        typedef typename mesh_type::face        face_type;
        typedef typename mesh_type::scalar_type scalar_type;
        typedef typename mesh_type::point_type  point_type;


        auto f = make_scalar_testing_data(msh);

        typename revolution::hho_degree_info hdi(degree);

        scalar_type error = 0.0;

        typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;
        boundary_type  bnd(msh);
        bnd.addDirichletEverywhere(f); //TOP
        auto has_boundary_vector = make_is_dirichlet_vector(msh, bnd);

        auto cl_count = 0;
        for (auto& cl : msh)
        {
            if(has_boundary_vector.at(cl_count) == 1)
            {
                Matrix<scalar_type, Dynamic, 1> proj = revolution::project_function(msh, cl, hdi, f);
                auto gr = make_hho_nitshce_scalar_laplacian(msh, cl, hdi, bnd);

                size_t rec_size = revolution::scalar_basis_size(hdi.reconstruction_degree(), Mesh::dimension);

                Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
                reconstr.tail(rec_size-1) = gr.first * proj;
                reconstr(0) = proj(0);

                auto cb = revolution::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
                Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);

                auto quad_degree = 2 * std::max(hdi.cell_degree(), hdi.reconstruction_degree());
                auto qps = revolution::integrate(msh, cl, quad_degree);

                for(auto& qp : qps)
                {
                    auto c_phi = cb.eval_functions(qp.point());
                    auto rec_f = reconstr.dot(c_phi);
                    auto diff  = (rec_f - f(qp.point()));
                    error += qp.weight() * diff * diff;
                }

            }
            else
            {
                Matrix<scalar_type, Dynamic, 1> proj = revolution::project_function(msh, cl, hdi, f);
                auto gr = revolution::make_hho_scalar_laplacian(msh, cl, hdi);

                size_t rec_size = revolution::scalar_basis_size(hdi.reconstruction_degree(), Mesh::dimension);

                Matrix<scalar_type, Dynamic, 1> reconstr = Matrix<scalar_type, Dynamic, 1>::Zero(rec_size);
                reconstr.tail(rec_size-1) = gr.first * proj;
                reconstr(0) = proj(0);

                auto cb = revolution::make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
                Matrix<scalar_type, Dynamic, Dynamic> mass = revolution::make_mass_matrix(msh, cl, cb);

                auto quad_degree = 2 * std::max(hdi.cell_degree(), hdi.reconstruction_degree());
                auto qps = revolution::integrate(msh, cl, quad_degree);

                for(auto& qp : qps)
                {
                    auto c_phi = cb.eval_functions(qp.point());
                    auto rec_f = reconstr.dot(c_phi);
                    auto diff  = (rec_f - f(qp.point()));
                    error += qp.weight() * diff * diff;
                }
            }
            cl_count++;
        }

        return std::sqrt(error);
    }

    size_t
    expected_rate(size_t k)
    {
        return k+1;
    }
};


int main(int argc, char **argv)
{
    tester<test_functor> tstr;
    int ch;

    while ( (ch = getopt(argc, argv, "sc")) != -1 )
    {
        switch(ch)
        {
            case 's':
                tester<test_functor> tstr_scalar;
                tstr_scalar.run();
            break;
            case 'c':
                tester<test_functor_contact> tstr_contact;
                tstr_contact.run();
        }
    }

    return 0;
}
