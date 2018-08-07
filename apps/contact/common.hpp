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
    const size_t DIM = Mesh::dimension;
    auto gamma_N = gamma_0 / diameter(msh, cl);

    auto fcs = faces(msh, cl);
    auto rbs = scalar_basis_size(hdi.reconstruction_degree(), DIM);
    auto cbs = scalar_basis_size(hdi.cell_degree(), DIM);
    auto fbs = scalar_basis_size(hdi.face_degree(), DIM-1);
    auto num_total_dofs = cbs + fbs * fcs.size();
    auto cb  = make_scalar_monomial_basis(msh, cl, hdi.reconstruction_degree());

    Matrix<T, Dynamic, Dynamic> ret =
            Matrix<T, Dynamic, Dynamic>::Zero( num_total_dofs, num_total_dofs );

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

            auto qps = integrate(msh, fc, 2* std::max(hdi.face_degree(), hdi.face_degree()));
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
