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
//#include "solvers/solver.hpp"

#include "cfem/cfem.hpp"

using namespace revolution;

enum method_type
{
    POSITIVE,
    NEGATIVE,
    ZERO
};

enum eval_solver_type
{
    EVAL_IN_CELLS,
    EVAL_ON_FACES,
    EVAL_IN_CELLS_AS_FACES,
    EVAL_WITH_PARAMETER
};

template<typename T>
void
dump_matrix(const dynamic_vector<T>& vec, const std::string& filename)
{
    std::ofstream ofs(filename);

    size_t i = 0;
    for (size_t j = 0;  j < vec.size(); j++)
        ofs << std::setprecision(16) << vec(j) << std::endl;
    ofs.close();
}

template<typename T>
struct algorithm_parameters
{
    algorithm_parameters():gamma_0(0.1), theta(-1.), solver(EVAL_ON_FACES), degree(1)
    {}

    T gamma_0;
    T theta;
    T degree;
    eval_solver_type solver;

    friend std::ostream& operator<<(std::ostream& os, const algorithm_parameters<T>& p)
    {
        os << "Algorithm parameters: "<<std::endl;
        os << "* gamma_0 : "<< p.gamma_0<< std::endl;
        os << "* theta   : "<< p.theta << std::endl;
        os << "* degree  : "<< p.degree<< std::endl;

        switch (p.solver)
        {
            case EVAL_IN_CELLS:
                os << "* Values : IN CELLS"<< std::endl;
                break;
            case EVAL_ON_FACES:
                os << "* Values : ON FACES"<< std::endl;
                break;
            case EVAL_IN_CELLS_AS_FACES:
                os << "* Values : IN CELLS AS FACES"<< std::endl;
                break;
            case EVAL_WITH_PARAMETER:
                os << "* Values : IN CELLS AND FACES (PARAMETER)"<< std::endl;
                break;
            default:
                os << "* Values : NOT SPECIFIED"<< std::endl;
                exit(1);
        }

        return os;
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

template<typename Mesh>
auto
full_offset(const Mesh& msh, const hho_degree_info& hdi)
{
    size_t  i = 1, dofs = 0;
    auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
    auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);

    std::vector<size_t> vec(msh.cells_size());

    for (auto itor = msh.cells_begin(); itor != msh.cells_end() - 1; itor++)
    {
        auto cl = *itor;
        dofs += howmany_faces(msh, cl) * fbs + cbs;
        vec.at(i++) = dofs;
    }
    return vec;
}

template<typename Mesh>
auto
make_is_contact_vector(const Mesh& msh,
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

            if (bnd.is_contact_face(face_id))
            {
                ret.at(i) = 1;
                continue;
            }
        }
        i++;
    }
    return ret;
}

template<typename T>
static_matrix<T, 3, 3>
make_fem_nitsche(const disk::simplicial_mesh<T, 2>& msh,
        const typename disk::simplicial_mesh<T, 2>::cell& cl,
        const disk::mechanics::BoundaryConditionsScalar<disk::simplicial_mesh<T, 2>>& bnd,
        const T & gamma_0,
        const T & theta)
{
    static_matrix<T, 3, 3> stiff = static_matrix<T, 3, 3>::Zero();

    auto gamma_N = gamma_0 / diameter(msh, cl);
    auto fcs = faces(msh, cl);

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id = eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto meas = measure(msh, fc);
            auto n    = normal(msh, cl, fc);
            auto dphi = disk::cfem::eval_basis_grad(msh, cl);
            static_vector<T, 3> dphi_n = dphi * n;
            stiff += meas * dphi_n * dphi_n.transpose();
        }
    }
    return stiff * (theta /gamma_N);
}

template <typename Mesh>
Matrix< typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_nitsche(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& rec,
            const typename Mesh::coordinate_type& gamma_0,
            const typename Mesh::coordinate_type& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();
    auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id = eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto n = normal(msh, cl, fc);
            auto qps = integrate(msh, fc, 2 * hdi.reconstruction_degree() - 2);

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                Matrix<T, Dynamic, 1>  sigma_n = rec.transpose() * (c_dphi * n);

                ret +=  qp.weight() * sigma_n * sigma_n.transpose();
            }
        }
    }

    return ret * (theta /gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_consist_diff(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T& gamma_0,
            const T& theta)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto celdeg =  hdi.cell_degree();
    auto facdeg =  hdi.face_degree();

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto cbs = scalar_basis_size(celdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, recdeg);
            auto n   = normal(msh, cl, fc);

            auto quad_degree = std::max(recdeg-1, celdeg);
            auto qps = integrate(msh, fc, 2 * quad_degree );
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                // grad * rec * u * n
                vector_type  sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_phi_temp.block(0, 0, cbs, 1);

                ret.block(0,0, cbs, num_total_dofs) += qp.weight() * c_phi * sigma_n.transpose();
            }
        }
    }
    return ret;
}


template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_consist_diff_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T& gamma_0,
            const T& theta)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto qps = integrate(msh, fc, 2* hdi.face_degree());
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  f_phi   = fb.eval_functions(qp.point());

                ret.block( face_ofs, 0, fbs, num_total_dofs) += qp.weight() * f_phi * sigma_n.transpose();
            }
        }
        fc_count++;
    }
    return ret;
}


template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_consist_diff_par(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T& gamma_0,
            const T& theta,
            const T& eta)
{

    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto quad_degree = std::max(hdi.face_degree(), hdi.cell_degree());
            auto qps = integrate(msh, fc, 2 * quad_degree + 2);
            for (auto& qp : qps)
            {

                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);

                vector_type  f_phi = fb.eval_functions(qp.point());
                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_phi_temp.block(0, 0, cbs, 1);

                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  vF =  (1. - eta) * f_phi;
                vector_type  vT = eta * c_phi_temp.block(0, 0, cbs, 1);

                //v = eta * vT + (1 -eta ) * vF
                vector_type v = vector_type::Zero(num_total_dofs);
                v.block(0, 0,cbs, 1)  = vT;
                v.block(face_ofs,0,fbs, 1) = vF;

                ret += qp.weight() * v * sigma_n.transpose();
                //ret += qp.weight() *  sigma_n * v.transpose();
            }
        }
        fc_count++;
    }
    return ret;
}


template <typename Mesh, typename T, typename Function>
std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, 1>>
make_hho_nitsche_diff(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const Function& dirichlet_bf)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto celdeg =  hdi.cell_degree();
    auto facdeg =  hdi.face_degree();

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto cbs = scalar_basis_size(celdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type Aret = matrix_type::Zero( num_total_dofs, num_total_dofs );
    vector_type Bret = vector_type::Zero( num_total_dofs);

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
            auto n   = normal(msh, cl, fc);

            auto quad_degree = std::max(recdeg-1, celdeg);
            auto qps = integrate(msh, fc, 2 * quad_degree );
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_phi_temp.block(0, 0, cbs, 1);

                // (theta * grad * rec * v * n - gamma_N* v_T)
                vector_type t_sigmav_n_gamma_v = theta * sigma_n;
                t_sigmav_n_gamma_v.block(0, 0,cbs, 1) -=  gamma_N * c_phi;

                Aret.block(0, 0, num_total_dofs, cbs) += qp.weight() * (t_sigmav_n_gamma_v) * c_phi.transpose();
                Bret += qp.weight() * (t_sigmav_n_gamma_v) * dirichlet_bf(qp.point());
            }
        }
    }

    return std::make_pair(Aret, Bret);
}

template <typename Mesh, typename T, typename Function>
std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, 1>>
make_hho_nitsche_diff_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const Function& dirichlet_bf)
{
    std::cout << " theta :"<< theta << std::endl;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type Aret = matrix_type::Zero( num_total_dofs, num_total_dofs );
    vector_type Bret = vector_type::Zero( num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto qps = integrate(msh, fc, 2* hdi.face_degree());
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  f_phi   = fb.eval_functions(qp.point());

                // (theta * grad * rec * v * n - gamma_N* v_F)
                vector_type t_sigmav_n_gamma_v = theta * sigma_n;
                t_sigmav_n_gamma_v.block(face_ofs, 0,fbs, 1) -=  gamma_N * f_phi;

                Aret.block(0, face_ofs, num_total_dofs, fbs) +=
                            qp.weight() * (t_sigmav_n_gamma_v) * f_phi.transpose();
                Bret += qp.weight() * (t_sigmav_n_gamma_v) * dirichlet_bf(qp.point());
            }
        }
        fc_count++;
    }
    return std::make_pair(Aret, Bret);
}


template <typename Mesh, typename T, typename Function>
std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, 1>>
make_hho_nitsche_diff_par(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T & theta,
            const Function& dirichlet_bf,
            const T& eta)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type Aret = matrix_type::Zero( num_total_dofs, num_total_dofs );
    vector_type Bret = vector_type::Zero( num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (msh.is_boundary(fc))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto quad_degree = std::max(hdi.face_degree(), hdi.cell_degree());
            auto qps = integrate(msh, fc, 2* quad_degree+2);
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  c_phi_temp   = cb.eval_functions(qp.point());

                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  vF  =  (1. - eta) * fb.eval_functions(qp.point());
                vector_type  vT  =  eta * c_phi_temp.block(0, 0, cbs, 1);

                //v = eta * vT + (1 -eta ) * vF
                vector_type u = vector_type::Zero(num_total_dofs);
                u.block(0, 0, cbs, 1)  = vT;
                u.block(face_ofs,0,fbs, 1) = vF;
                vector_type v = u;

                // (theta * grad * rec * v * n - gamma_N* v)
                vector_type t_sigmav_n_gamma_v = theta * sigma_n;
                t_sigmav_n_gamma_v -=  gamma_N * v;

                // u  ( theta * sigma_n - gamma_N *v)
                Aret += qp.weight() * (t_sigmav_n_gamma_v) * u.transpose();
                //Aret += qp.weight() * u * (t_sigmav_n_gamma_v).transpose();

                // u  ( theta * sigma_n - gamma_N *v)
                Bret += qp.weight() * (t_sigmav_n_gamma_v) * dirichlet_bf(qp.point());
            }
        }
        fc_count++;
    }
    return std::make_pair(Aret, Bret);
}


template<typename T>
static_matrix<T, 3, 3>
make_fem_heaviside(const disk::simplicial_mesh<T, 2>& msh,
        const typename disk::simplicial_mesh<T, 2>::cell& cl,
        const disk::mechanics::BoundaryConditionsScalar<disk::simplicial_mesh<T, 2>>& bnd,
        const T & gamma_0,
        const T & theta,
        const static_vector<T,3>& uloc)
{
    static_matrix<T, 3, 3> ret = static_matrix<T, 3, 3>::Zero();

    auto gamma_N = gamma_0 / diameter(msh, cl);
    auto fcs = faces(msh, cl);

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id = eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto n   = normal(msh, cl, fc);
            static_matrix<T, 3, 2> dphi = disk::cfem::eval_basis_grad(msh, cl);
            static_vector<T, 3> sigma_n = dphi * n;

            //(theta * grad * v * n - gamma_N * v)
            static_vector<T, 3> t_sigmav_n_gamma_v =  theta * sigma_n;
            static_vector<T, 3>   sigmau_n_gamma_u =  sigma_n;

            auto qps = integrate (msh, fc, 2);
            for (auto& qp : qps)
            {
                auto phi = disk::cfem::eval_basis(msh, cl, qp.point());

                t_sigmav_n_gamma_v -= gamma_N * phi;
                sigmau_n_gamma_u   -= gamma_N * phi;

                // Heaviside(-Pn(u) = grad * u * n - gamma_N* u)
                if(sigmau_n_gamma_u.dot(uloc) <= 0.)
                    ret += qp.weight() * t_sigmav_n_gamma_v * sigmau_n_gamma_u.transpose();
            }
        }
    }
    return ret * (1/gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_heaviside_other(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<T, Dynamic, 1>& uloc)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 /diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);
            auto qps = integrate(msh, fc, 2* hdi.face_degree());
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  f_phi   = fb.eval_functions(qp.point());
                vector_type  u_face  = uloc.block(face_ofs, 0, fbs, 1);

                // Heaviside(-Pn(u) = -D*recu*n + gN*u_F)
                T gamma_u  = gamma_N * f_phi.dot(u_face);
                T sigmau_n = sigma_n.dot(uloc);

                if (sigmau_n - gamma_u  <= 0.)
                {
                    // (theta * sigma * n * sigma * n)
                    matrix_type t_sigmavn_sigmaun = theta * sigma_n * sigma_n.transpose();
                    matrix_type g_vF_sigmaun = gamma_N * f_phi * sigma_n.transpose();
                    matrix_type t_sigmavn_g_uF = theta * gamma_N * sigma_n * f_phi.transpose();
                    matrix_type g2_vF_uF =  gamma_N * gamma_N * f_phi * f_phi.transpose();

                    ret += qp.weight() * t_sigmavn_sigmaun;
                    ret.block(face_ofs, 0, fbs, num_total_dofs) -= qp.weight() * g_vF_sigmaun;
                    ret.block(0, face_ofs, num_total_dofs, fbs) -= qp.weight() * t_sigmavn_g_uF;
                    ret.block(face_ofs, face_ofs, fbs, fbs) += qp.weight() * g2_vF_uF;
                }
            }
        }
        fc_count++;
    }
    return ret * (1./gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_heaviside(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<T, Dynamic, 1>& uloc)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto celdeg =  hdi.cell_degree();
    auto facdeg =  hdi.face_degree();

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto cbs = scalar_basis_size(celdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
            auto n   = normal(msh, cl, fc);

            auto quad_degree = std::max(recdeg-1, celdeg);
            auto qps = integrate(msh, fc, 3 * quad_degree );

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_phi_temp.block(0, 0, cbs, 1);

                // (grad * rec * u * n - gamma_N* u_T)
                vector_type sigmau_n_gamma_u  = sigma_n;
                sigmau_n_gamma_u.block(0, 0, cbs, 1) -= gamma_N * c_phi;

                // Heaviside(-Pn(u) = -D*recu*n + gN*u_T)
                if(sigmau_n_gamma_u.dot(uloc)  <= 0.)
                {
                    // (theta * grad * rec * v * n - gamma_N* v_T)

                    vector_type t_sigmav_n_gamma_v = vector_type::Zero(num_total_dofs);
                    t_sigmav_n_gamma_v = theta * sigma_n;
                    t_sigmav_n_gamma_v.block(0, 0,cbs, 1) -=  gamma_N * c_phi;

                    ret += qp.weight() * (t_sigmav_n_gamma_v) * (sigmau_n_gamma_u).transpose();
                }
            }
        }
    }
    return ret * (1./gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_heaviside_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<T, Dynamic, 1>& uloc)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto qps = integrate(msh, fc, 2* hdi.face_degree());
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  f_phi   = fb.eval_functions(qp.point());
                vector_type  u_face  = uloc.block(face_ofs, 0, fbs, 1);

                // Heaviside(-Pn(u) = -D*recu*n + gN*u_F)
                T gamma_u  = gamma_N * f_phi.dot(u_face);
                T sigmau_n = sigma_n.dot(uloc);

                if (sigmau_n - gamma_u  <= 0.)
                {
                    // (theta * grad * rec * v * n - gamma_N* v_F)
                    vector_type t_sigman_g_v = theta * sigma_n;
                    t_sigman_g_v.block(face_ofs, 0,fbs, 1) -=  gamma_N * f_phi;

                    // (grad * rec * u * n - gamma_N* u_F)
                    vector_type sigman_g_u = sigma_n;
                    sigman_g_u.block(face_ofs, 0,fbs, 1) -=  gamma_N * f_phi;

                    ret += qp.weight() * (t_sigman_g_v) * (sigman_g_u).transpose();
                }
            }
        }
        fc_count++;
    }
    return ret * (1./gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_heaviside_trace(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<T, Dynamic, 1>& uloc)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );
    matrix_type rec_full = matrix_type::Zero(rbs, num_total_dofs );

    rec_full(0,0) = 1.;
    rec_full.block(1,0, rbs -1, num_total_dofs)  =  rec;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto n = normal(msh, cl, fc);
            matrix_type stiff_n = matrix_type::Zero(rbs -1, rbs -1);
            matrix_type mm  = revolution::make_mass_matrix(msh, fc, cb);

            auto qps = integrate(msh, fc, 2 * recdeg);
            for (auto& qp : qps)
            {
                auto c_phi = cb.eval_functions(qp.point());

                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = (c_dphi * n).transpose() * rec;
                vector_type gamma_rec =  gamma_N * rec_full.transpose() * c_phi;

                // Heaviside(-Pn(u) = -D*recu*n + gN*u)
                T sigmau_n   = sigma_n.dot(uloc);
                T gamma_recu = gamma_rec.dot(uloc);

                if (sigmau_n - gamma_recu  <= 0.)
                {
                    // (theta * grad * rec * u * n - gamma_N* rec *u)
                    vector_type t_sigman_g_rec = theta * sigma_n - gamma_rec;
                    //t_sigman_g_rec.block(1, 0,rbs -1, 1) += theta * sigma_n;

                    // (grad * rec * v * n - gamma_N* rec * v)
                    vector_type sigman_g_rec = sigma_n - gamma_rec;
                    //sigman_g_rec.block(1, 0,rbs -1, 1) += sigma_n;

                    ret += qp.weight() * (t_sigman_g_rec) * (sigman_g_rec).transpose();
                }
            }
        }
    }
    return ret * (1./gamma_N);
}

template <typename T>
static_vector<T, 3>
make_fem_negative(const disk::simplicial_mesh<T, 2>& msh,
    const typename disk::simplicial_mesh<T, 2>::cell& cl,
    const disk::mechanics::BoundaryConditionsScalar<disk::simplicial_mesh<T, 2>>& bnd,
    const T & gamma_0,
    const T & theta,
    const static_vector<T, 3>& uloc)
{
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    static_vector<T, 3> rhs = static_vector<T, 3>::Zero();

    for (auto& fc : fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto n   = normal(msh, cl, fc);
            static_matrix<T, 3, 2> dphi = disk::cfem::eval_basis_grad(msh, cl);
            static_vector<T, 3> sigma_n = dphi * n;

            //(theta * grad * v * n - gamma_N * v)
            static_vector<T, 3> t_sigmav_n_gamma_v =  theta * sigma_n;
            static_vector<T, 3>   sigmau_n_gamma_u =  sigma_n;

            auto qps = integrate (msh, fc, 2);
            for (auto& qp : qps)
            {
                auto phi = disk::cfem::eval_basis(msh, cl, qp.point());

                t_sigmav_n_gamma_v -= gamma_N * phi;
                sigmau_n_gamma_u   -= gamma_N * phi;

                // [Pn(u)]_ = [grad * u * n - gamma_N* u]_
                T negative = std::min(sigmau_n_gamma_u.dot(uloc), 0.);

                rhs += qp.weight() * negative * t_sigmav_n_gamma_v;
            }
        }
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_hho_negative(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& rec,
            const typename Mesh::coordinate_type& gamma_0,
            const typename Mesh::coordinate_type& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<typename Mesh::coordinate_type, Dynamic, 1>& uloc)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto facdeg =  hdi.face_degree();
    auto celdeg =  hdi.cell_degree();

    auto fcs = faces(msh, cl);
    auto cbs = scalar_basis_size(celdeg, DIM);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fcs.size() * fbs;
    vector_type rhs = vector_type::Zero(num_total_dofs, 1);

    for (auto& fc : fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id = eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, recdeg);
            auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
            auto quad_degree = std::max( recdeg - 1, celdeg);
            auto qps = integrate (msh, fc, 3 * quad_degree);
            auto n = normal(msh, cl, fc);

            for (auto& qp : qps)
            {

                Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic,  DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_phi_temp.block(0, 0, cbs, 1);

                vector_type  sigma_n = rec.transpose() * (c_dphi * n);

                // (grad * rec * v * n - gamma_N * v_T)
                vector_type t_sigmav_n_gamma_v =  theta * sigma_n;
                t_sigmav_n_gamma_v.block(0, 0, cbs, 1) -= gamma_N * c_phi;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* u_T]_
                vector_type sigmau_n_gamma_u  = sigma_n;
                sigmau_n_gamma_u.block(0, 0, cbs, 1) -= gamma_N * c_phi;

                T negative = std::min(sigmau_n_gamma_u.dot(uloc), 0.);

                rhs += qp.weight() * negative * t_sigmav_n_gamma_v;
            }
        }
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_hho_negative_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& rec,
            const typename Mesh::coordinate_type& gamma_0,
            const typename Mesh::coordinate_type& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<typename Mesh::coordinate_type, Dynamic, 1>& uloc)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto facdeg =  hdi.face_degree();
    auto fcs = faces(msh, cl);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fcs.size() * fbs;
    vector_type rhs = vector_type::Zero(num_total_dofs, 1);

    auto fc_count = 0;
    for (auto& fc:fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, recdeg);
            auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
            auto qps = integrate (msh, fc, 2*facdeg);
            auto n = normal(msh, cl, fc);
            auto face_ofs  = cbs  +  fc_count * fbs;

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic,  DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  f_phi   = fb.eval_functions(qp.point());
                vector_type  sigma_n = rec.transpose() * (c_dphi * n);

                // (grad * rec * v * n - gamma_N * v_F)
                vector_type t_sigmav_n_gamma_v =  theta * sigma_n;
                t_sigmav_n_gamma_v.block(face_ofs, 0, fbs, 1) -= gamma_N * f_phi;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* u_F]_
                vector_type sigmau_n_gamma_u  = sigma_n;
                sigmau_n_gamma_u.block(face_ofs, 0, fbs, 1) -= gamma_N * f_phi;

                T negative = std::min(sigmau_n_gamma_u.dot(uloc), 0.);

                rhs += qp.weight() * negative * t_sigmav_n_gamma_v;
            }
        }
        fc_count++;
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_hho_negative_trace(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& rec,
            const typename Mesh::coordinate_type& gamma_0,
            const typename Mesh::coordinate_type& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<typename Mesh::coordinate_type, Dynamic, 1>& uloc)
{
    using T = typename Mesh::coordinate_type;
    using vector_type = Matrix<T, Dynamic, 1>;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;

    const size_t DIM = Mesh::dimension;

    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto facdeg =  hdi.face_degree();
    auto fcs = faces(msh, cl);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fcs.size()* fbs;

    vector_type rhs = vector_type::Zero(num_total_dofs, 1);
    matrix_type rec_full = matrix_type::Zero(rbs, num_total_dofs );
    rec_full(0,0) = 1.;
    rec_full.block(1,0, rbs -1, num_total_dofs)  =  rec;

    auto fc_count = 0;
    for (auto& fc:fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
            auto cb = make_scalar_monomial_basis(msh, cl, recdeg);
            auto qps = integrate (msh, fc, 2 * recdeg);
            auto n = normal(msh, cl, fc);
            auto face_ofs  = cbs  +  fc_count * fbs;

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic,  DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type     c_phi = cb.eval_functions(qp.point());

                vector_type sigma_n   =  rec.transpose() * (c_dphi * n) ;
                vector_type gamma_rec = gamma_N * rec_full.transpose() * c_phi;

                assert((rec_full * uloc).size() == rbs);

                // (grad * rec * v * n - gamma_N* rec * v)
                vector_type t_sigman_g_recv = theta * sigma_n - gamma_rec;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* rec * u]_
                T sigmau_n   = sigma_n.dot(uloc);
                T gamma_recu = gamma_rec.dot(uloc);
                T negative   = std::min(sigmau_n - gamma_recu, 0.);

                rhs += qp.weight() * negative * t_sigman_g_recv;
            }
        }
        fc_count++;
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_hho_heaviside_par(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<T, Dynamic, Dynamic>& rec,
            const T & gamma_0,
            const T& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<T, Dynamic, 1>& uloc,
            const typename Mesh::coordinate_type& eta)
{
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    T gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();

    matrix_type ret = matrix_type::Zero( num_total_dofs, num_total_dofs );

    auto fc_count = 0;

    for (auto& fc: fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());
            auto fb  = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

            auto face_ofs  = cbs  +  fc_count * fbs;
            auto n = normal(msh, cl, fc);

            auto qps = integrate(msh, fc, 2* std::max(hdi.face_degree(), hdi.cell_degree()));
            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi   = c_phi_temp.block(0, 0, cbs, 1);
                vector_type  u_cell  = uloc.block(0, 0, cbs, 1);

                vector_type  f_phi   = fb.eval_functions(qp.point());
                vector_type  u_face  = uloc.block(face_ofs, 0, fbs, 1);

                // Heaviside(-Pn(u) = -D*recu*n +  gamma_N* (eta* u_T + (1-eta)*uF)
                T uF_1_eta  = (1. - eta) * f_phi.dot(u_face);
                T uT_eta    = eta * c_phi.dot(u_cell);
                T sigmau_n = sigma_n.dot(uloc);

                if (sigmau_n - gamma_N * (uT_eta +  uF_1_eta)  <= 0.)
                {
                    // (theta * grad * rec * v * n - gamma_N*(eta* v_T + (1-eta)*vF))
                    vector_type t_sigman_g_v = theta * sigma_n;
                    t_sigman_g_v.block(face_ofs, 0,fbs, 1) -= (1. - eta) *  gamma_N * f_phi;
                    t_sigman_g_v.block(0, 0,cbs, 1) -=  eta * gamma_N * c_phi;

                    // (grad * rec * u * n - gamma_N* (eta* u_T + (1-eta)*uF)
                    vector_type sigman_g_u = sigma_n;
                    sigman_g_u.block(face_ofs, 0,fbs, 1) -= (1. - eta) * gamma_N * f_phi;
                    sigman_g_u.block(0, 0,cbs, 1) -=  eta * gamma_N * c_phi;

                    ret += qp.weight() * (t_sigman_g_v) * (sigman_g_u).transpose();
                }
            }
        }
        fc_count++;
    }
    return ret * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_hho_negative_par(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& rec,
            const typename Mesh::coordinate_type& gamma_0,
            const typename Mesh::coordinate_type& theta,
            const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
            const Matrix<typename Mesh::coordinate_type, Dynamic, 1>& uloc,
            const typename Mesh::coordinate_type& eta)
{
    using T = typename Mesh::coordinate_type;
    using matrix_type = Matrix<T, Dynamic, Dynamic>;
    using vector_type = Matrix<T, Dynamic, 1>;

    const size_t DIM = Mesh::dimension;
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto recdeg =  hdi.reconstruction_degree();
    auto facdeg =  hdi.face_degree();
    auto celdeg =  hdi.cell_degree();

    auto fcs = faces(msh, cl);
    auto cbs = scalar_basis_size(celdeg, DIM);
    auto rbs = scalar_basis_size(recdeg, DIM);
    auto fbs = scalar_basis_size(facdeg, DIM-1);
    auto num_total_dofs = cbs + fcs.size() * fbs;
    vector_type rhs = vector_type::Zero(num_total_dofs, 1);

    auto fc_count = 0;
    for (auto& fc:fcs)
    {
        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
        {
            auto cb  = make_scalar_monomial_basis(msh, cl, recdeg);
            auto fb  = make_scalar_monomial_basis(msh, fc, facdeg);
            auto qps = integrate (msh, fc, 2*std::max(facdeg, celdeg));
            auto n = normal(msh, cl, fc);
            auto face_ofs  = cbs  +  fc_count * fbs;

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic,  DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi   = c_phi_temp.block(0, 0, cbs, 1);
                vector_type  f_phi   = fb.eval_functions(qp.point());
                vector_type  sigma_n = rec.transpose() * (c_dphi * n);

                // ( theta * grad * rec * v * n - gamma_N* (eta* v_T + (1-eta)*vF))
                vector_type t_sigmav_n_gamma_v =  theta * sigma_n;
                t_sigmav_n_gamma_v.block(face_ofs, 0, fbs, 1) -= (1. - eta) * gamma_N * f_phi;
                t_sigmav_n_gamma_v.block(0, 0, cbs, 1) -= eta * gamma_N * c_phi;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* (eta* u_T + (1-eta)*uF)]_
                vector_type sigmau_n_gamma_u  = sigma_n;
                sigmau_n_gamma_u.block(face_ofs, 0, fbs, 1) -= (1. - eta) * gamma_N * f_phi;
                sigmau_n_gamma_u.block(0, 0, cbs, 1) -= eta * gamma_N * c_phi;

                T negative = std::min(sigmau_n_gamma_u.dot(uloc), 0.);

                rhs += qp.weight() * negative * t_sigmav_n_gamma_v;
            }
        }
        fc_count++;
    }
    return rhs * (1./gamma_N);
}


template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_contact_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di,
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

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2 * (recdeg-1));
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

        if (bnd.is_contact_face(face_id))
            continue;

        const auto n  = normal(msh, cl, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        size_t quad_degree = std::max(recdeg - 1 + std::max(facdeg,celdeg), size_t(0));
        auto qps_f = integrate(msh, fc, quad_degree);
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
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_nitshce_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di,
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

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    const auto num_faces = howmany_faces(msh, cl);

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type gr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type gr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2 * (recdeg-1));
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
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        size_t quad_degree = std::max(recdeg - 1 + std::max(facdeg,celdeg), size_t(0));
        auto qps_f = integrate(msh, fc, quad_degree);
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
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hdg_nitsche_stabilization(const Mesh& msh,
        const typename Mesh::cell_type& cl, const hho_degree_info& di,
        const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto num_faces = howmany_faces(msh, cl);
    const auto total_dofs = cbs + num_faces * fbs;

    matrix_type         data = matrix_type::Zero(total_dofs, total_dofs);
    const matrix_type   If   = matrix_type::Identity(fbs, fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto fcs = faces(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        const auto fc = fcs[i];

        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_dirichlet_face(face_id))
            continue;

        const auto h  = measure(msh, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type    oper  = matrix_type::Zero(fbs, total_dofs);
        matrix_type    tr    = matrix_type::Zero(fbs, total_dofs);
        matrix_type    mass  = matrix_type::Zero(fbs, fbs);
        matrix_type    trace = matrix_type::Zero(fbs, cbs);

        oper.block(0, cbs + i  * fbs, fbs, fbs) = -If;

        const auto qps = integrate(msh, fc, facdeg + celdeg);
        for (auto& qp : qps)
        {
            const auto c_phi = cb.eval_functions(qp.point());
            const auto f_phi = fb.eval_functions(qp.point());

            assert(c_phi.rows() == cbs);
            assert(f_phi.rows() == fbs);
            assert(c_phi.cols() == f_phi.cols());

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        tr.block(0, cbs + i * fbs, fbs, fbs) = -mass;
        tr.block(0, 0, fbs, cbs) = trace;

        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);
        data += oper.transpose() * tr * (1./h);
    }

    return data;
}


template<typename Mesh>
class diffusion_condensed_assembler_nitsche_cells
{
    using T = typename Mesh::coordinate_type;

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;

    boundary_type                       m_bnd;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;
    std::vector< Triplet<T> >           triplets;

    size_t       num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t       system_size;


    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

    public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    diffusion_condensed_assembler_nitsche_cells(const Mesh& msh, hho_degree_info hdi,
                                                const boundary_type& bnd)
        : di(hdi),  m_bnd(bnd)
    {
        auto is_dirichlet = [&](const typename Mesh::face& fc) -> bool {
            auto fc_id = msh.lookup(fc);
            return bnd.is_dirichlet_face(fc_id);
        };

        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ( !is_dirichlet(fc) )
            {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }


        auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
        auto system_size = fbs * num_other_faces;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );

    }

    #if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
    #endif

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            auto face_id = msh.lookup(fc);
            bool dirichlet = m_bnd.is_dirichlet_face(face_id);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
            }

            RHS(asm_map[i]) += rhs(i);
        }
    }


    Matrix<T, Dynamic,1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
            const Matrix<T, Dynamic, 1>& solution)

    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

            if (m_bnd.is_dirichlet_face(face_id))
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                auto dirichlet_bf = m_bnd.dirichlet_boundary_func(face_id);
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = priv::offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};

template<typename Mesh>
class diffusion_condensed_assembler_nitsche_faces
{
    using T = typename Mesh::coordinate_type;

    typedef disk::mechanics::BoundaryConditionsScalar<Mesh> boundary_type;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    hho_degree_info                     di;
    std::vector< Triplet<T> >           triplets;

    size_t       num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t       system_size;


    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

    public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    diffusion_condensed_assembler_nitsche_faces(const Mesh& msh, hho_degree_info hdi)
        : di(hdi)
    {
        num_all_faces = msh.faces_size();
        num_other_faces = num_all_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            compress_table.at(i) = compressed_offset;
            expand_table.at(compressed_offset) = i;
            compressed_offset++;
        }

        auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);
        auto system_size = fbs * num_other_faces;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );

    }

    #if 0
    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }
    #endif

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs,
             const Matrix<T, Dynamic, 1>& rhs)
    {
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = priv::offset(msh, fc);
            auto face_LHS_offset = compress_table.at(face_offset)*fbs;

            auto face_id = msh.lookup(fc);
            bool dirichlet = false;

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    throw std::invalid_argument("It shouldn't come here!!");
            }

            RHS(asm_map[i]) += rhs(i);
        }
    }


    Matrix<T, Dynamic,1>
    take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl,
            const Matrix<T, Dynamic, 1>& solution)

    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
            const auto face_id=eid.second;

            auto face_offset = priv::offset(msh, fc);
            auto face_SOL_offset = compress_table.at(face_offset)*fbs;
            ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};

template<typename Mesh>
auto make_diffusion_assembler_nitsche_cells(const Mesh& msh, hho_degree_info hdi,
                const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd,
                const algorithm_parameters<typename Mesh::scalar_type>& ap)
{
    return diffusion_condensed_assembler_nitsche_cells<Mesh>(msh, hdi, bnd);
}

template<typename Mesh>
auto make_diffusion_assembler_nitsche_faces(const Mesh& msh, hho_degree_info hdi,
                const algorithm_parameters<typename Mesh::scalar_type>& ap)
{
    return diffusion_condensed_assembler_nitsche_faces<Mesh>(msh, hdi);
}
