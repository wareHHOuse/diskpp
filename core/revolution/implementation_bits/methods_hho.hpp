/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#pragma once

#include "revolution/bases"
#include "revolution/quadratures"

namespace revolution {

class hho_degree_info
{
    size_t  cell_deg, face_deg, reconstruction_deg;

public:
    hho_degree_info()
        : cell_deg(1), face_deg(1), reconstruction_deg(2)
    {}

    explicit hho_degree_info(size_t degree)
        : cell_deg(degree), face_deg(degree), reconstruction_deg(degree+1)
    {}

    hho_degree_info(size_t cd, size_t fd)
    {
        bool c1 = fd > 0  && (cd == fd-1 || cd == fd || cd == fd+1);
        bool c2 = fd == 0 && (cd == fd || cd == fd+1);
        if ( c1 || c2 )
        {
            cell_deg            = cd;
            face_deg            = fd;
            reconstruction_deg  = fd+1;

        }
        else
        {
            std::cout << "Invalid cell degree. Reverting to equal-order" << std::endl;
            cell_deg            = fd;
            face_deg            = fd;
            reconstruction_deg  = fd+1;
        }

        std::cout << cell_deg << " " << face_deg << " " << reconstruction_deg << std::endl;
    }

    size_t cell_degree() const
    {
        return cell_deg;
    }

    size_t face_degree() const
    {
        return face_deg;
    }

    size_t reconstruction_degree() const
    {
        return reconstruction_deg;
    }
};

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff += qp.second * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    auto fcs = faces(msh, cl);
    auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = ns[i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> c_phi_tmp = cb.eval_basis(qp.first);
            Matrix<T, Dynamic, 1> c_phi = c_phi_tmp.head(cbs);
            Matrix<T, Dynamic, 2> c_dphi_tmp = cb.eval_gradients(qp.first);
            Matrix<T, Dynamic, 2> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, 2);
            Matrix<T, Dynamic, 1> f_phi = fb.eval_basis(qp.first);
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.second * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.second * (c_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_naive_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    auto fcs = faces(msh, cl);
    auto h = measure(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+num_faces*fbs);
        Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> trace = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

        oper.block(0, cbs+i*fbs, fbs, fbs) = -If;

        auto qps = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);

            mass += qp.second * f_phi * f_phi.transpose();
            trace += qp.second * f_phi * c_phi.transpose();
        }

        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);

        data += oper.transpose() * mass * oper * (1./h);
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_fancy_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl,
                             const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> reconstruction,
                             const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension);

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    Matrix<T, Dynamic, Dynamic> mass_mat = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    auto cell_quadpoints = integrate(msh, cl, 2*recdeg);
    for (auto& qp : cell_quadpoints)
    {
        auto c_phi = cb.eval_basis(qp.first);
        mass_mat += qp.second * c_phi * c_phi.transpose();
    }

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    //Step 1: compute \pi_T^k p_T^k v (third term).
    Matrix<T, Dynamic, Dynamic> M1 = mass_mat.block(0, 0, cbs, cbs);
    Matrix<T, Dynamic, Dynamic> M2 = mass_mat.block(0, 1, cbs, rbs-1);
    Matrix<T, Dynamic, Dynamic> proj1 = -M1.llt().solve(M2*reconstruction);

    //Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    Matrix<T, Dynamic, Dynamic> I_T = Matrix<T, Dynamic, Dynamic>::Identity(cbs, cbs);
    proj1.block(0, 0, cbs, cbs) += I_T;

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto h = diameter(msh, /*fcs[face_i]*/cl);
        auto fc = fcs[face_i];
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> face_mass_matrix    = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> face_trace_matrix   = Matrix<T, Dynamic, Dynamic>::Zero(fbs, rbs);

        auto face_quadpoints = integrate(msh, fc, 2*facdeg);
        for (auto& qp : face_quadpoints)
        {
            auto f_phi = fb.eval_basis(qp.first);
            auto c_phi = cb.eval_basis(qp.first);
            auto q_f_phi = qp.second * f_phi;
            face_mass_matrix += q_f_phi * f_phi.transpose();
            face_trace_matrix += q_f_phi * c_phi.transpose();
        }

        LLT<Matrix<T, Dynamic, Dynamic>> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR1 = face_trace_matrix.block(0, 1, fbs, rbs-1);

        Matrix<T, Dynamic, Dynamic> proj2 = piKF.solve(MR1*reconstruction);
        Matrix<T, Dynamic, Dynamic> I_F = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);
        proj2.block(0, cbs+face_i*fbs, fbs, fbs) -= I_F;

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        Matrix<T, Dynamic, Dynamic> MR2 = face_trace_matrix.block(0, 0, fbs, cbs);
        Matrix<T, Dynamic, Dynamic> proj3 = piKF.solve(MR2*proj1);
        Matrix<T, Dynamic, Dynamic> BRF = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / h;
    }

    return data;
}


} // revolution





