/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */



namespace disk {

template<typename T>
Matrix<T, Dynamic, 3>
vcross(const static_vector<T, 3>& normal, const Matrix<T, Dynamic, 3>& field)
{
    Matrix<T, Dynamic, 3>   ret;
    ret = Matrix<T, Dynamic, 3>::Zero(field.rows(), field.cols());

    for (size_t i = 0; i < field.rows(); i++)
    {
        static_vector<T, 3> f;
        f(0) = field(i,0);
        f(1) = field(i,1);
        f(2) = field(i,2);

        auto r = normal.cross( f );
        ret(i,0) = r(0);
        ret(i,1) = r(1);
        ret(i,2) = r(2);
    }

    return ret;
}

template<typename T>
Matrix<T, Dynamic, 3>
vcross(const Matrix<T, Dynamic, 3>& field, const static_vector<T, 3>& normal)
{
    Matrix<T, Dynamic, 3>   ret;
    ret = Matrix<T, Dynamic, 3>::Zero(field.rows(), field.cols());

    for (size_t i = 0; i < field.rows(); i++)
    {
        static_vector<T, 3> f;
        f(0) = field(i,0);
        f(1) = field(i,1);
        f(2) = field(i,2);

        auto r = f.cross( normal );
        ret(i,0) = r(0);
        ret(i,1) = r(1);
        ret(i,2) = r(2);
    }

    return ret;
}

template<typename T>
Matrix<T, Dynamic, 1>
vcross(const static_vector<T, 2>& normal, const Matrix<T, Dynamic, 2>& field)
{
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(field.rows());

    for (size_t i = 0; i < field.rows(); i++)
        ret(i) = normal(0)*field(i,1) - normal(1)*field(i,0);

    return ret;
}

template<typename T>
Matrix<T, Dynamic, 1>
vcross(const Matrix<T, Dynamic, 2>& field, const static_vector<T, 2>& normal)
{
    return -vcross(normal, field);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_mass_oper(const Mesh&                       msh,
                      const typename Mesh::cell_type&   cl,
                      const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension-1);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type mass = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*celdeg);
    for (auto& qp : qps)
    {
        const auto phi = cb.eval_functions(qp.point());
        mass.block(0,0,cbs,cbs) += qp.weight() * phi * phi.transpose();
    }

    return mass;
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_hho_grad_impl(const Mesh&                       msh,
                          const typename Mesh::cell_type&   cl,
                          const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);

    auto fcs = faces(msh, cl);

    const auto num_face_dofs = fbs * fcs.size();

    const matrix_type gr_lhs = make_mass_matrix(msh, cl, rb);
    matrix_type       gr_rhs = matrix_type::Zero(rbs, cbs + num_face_dofs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for(auto& qp : qps)
    {
        const auto c_phi = cb.eval_functions( qp.point() );
        const auto r_phi = rb.eval_divergences( qp.point() );
        gr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * r_phi * c_phi.transpose();
    }

    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc     = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fbs    = scalar_basis_size(facdeg, Mesh::dimension - 1);
        const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto offset = cbs + i*fbs;

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            const vector_type f_phi = fb.eval_functions( qp.point() );
            const auto        r_phi = rb.eval_functions( qp.point() );
            gr_rhs.block(0, offset, rbs, fbs) += qp.weight() * (r_phi*n) * f_phi.transpose();
        }
    }

    return gr_lhs.llt().solve(gr_rhs);
}

template<typename T>
using matrixpair = std::pair<dynamic_matrix<T>, dynamic_matrix<T>>;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
matrixpair<typename Mesh<T,2,Storage>::coordinate_type>
curl_reconstruction(const Mesh<T,2,Storage>&                     msh,
                    const typename Mesh<T,2,Storage>::cell_type& cl,
                    const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<T,2,Storage>;
    using scalar_type       = typename mesh_type::coordinate_type;
    using matrix_type       = Matrix<T, Dynamic, Dynamic>;
    using vector_type       = Matrix<T, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;

    static const size_t DIM = 2;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto rb = make_scalar_monomial_basis(msh, cl, recdeg);
    const auto fbs = scalar_basis_size(facdeg, DIM - 1);
    const auto cbs = scalar_basis_size(celdeg, DIM);
    const auto rbs = scalar_basis_size(recdeg, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type cr_lhs = matrix_type::Zero(rbs-1, rbs-1);
    matrix_type cr_rhs = matrix_type::Zero(rbs-1, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        const auto cphi = rb.eval_curls2(qp.point());
        stiff += qp.weight() * cphi * cphi.transpose();
    }

    cr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    cr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_scalar_monomial_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 2> cphi      = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 1> cphi_n    = vcross(cphi, n).block(1,0,rbs-1,1);
            Matrix<T, Dynamic, 1> f_phi     = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> c_phi     = cb.eval_functions(qp.point());

            cr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * cphi_n * c_phi.transpose();
            cr_rhs.block(0, offset, rbs-1, fbs) += qp.weight() * cphi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    matrix_type oper = cr_lhs.ldlt().solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
matrixpair<typename Mesh<T,3,Storage>::coordinate_type>
curl_reconstruction(const Mesh<T,3,Storage>&                     msh,
                    const typename Mesh<T,3,Storage>::cell_type& cl,
                    const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<T,3,Storage>;
    using scalar_type       = typename mesh_type::coordinate_type;
    using matrix_type       = Matrix<T, Dynamic, Dynamic>;
    using vector_type       = Matrix<T, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;

    static const size_t DIM = 3;
    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);
    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs = vector_basis_size(facdeg, DIM - 1, DIM - 1);
    const auto cbs = vector_basis_size(celdeg, DIM, DIM);
    const auto rbs = vector_basis_size(recdeg, DIM, DIM);

    const auto      fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type cr_lhs = matrix_type::Zero(rbs-3, rbs-3);
    matrix_type cr_rhs = matrix_type::Zero(rbs-3, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        const auto cphi = rb.eval_curls2(qp.point());
        stiff += qp.weight() * cphi * cphi.transpose();
    }

    cr_lhs = stiff.block(3, 3, rbs-3, rbs-3);
    cr_rhs.block(0, 0, rbs-3, cbs) = stiff.block(3, 0, rbs-3, cbs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_vector_monomial_tangential_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> cphi      = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 3> cphi_n    = vcross(cphi, n).block(3,0,rbs-3,3);
            Matrix<T, Dynamic, 3> f_phi     = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

            cr_rhs.block(0, 0, rbs-3, cbs) -= qp.weight() * cphi_n * c_phi.transpose();
            cr_rhs.block(0, offset, rbs-3, fbs) += qp.weight() * cphi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    matrix_type oper = cr_lhs.ldlt().solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_curl_impl_pk(const Mesh&                       msh,
                          const typename Mesh::cell_type&   cl,
                          const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);
    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension-1);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);

    const auto      fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type mass = matrix_type::Zero(rbs, rbs);
    matrix_type cr_rhs = matrix_type::Zero(rbs, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        const auto phi = rb.eval_functions(qp.point());
        mass += qp.weight() * phi * phi.transpose();

        const auto cphi = cb.eval_curls2(qp.point());
        cr_rhs.block(0, 0, rbs, cbs) += qp.weight() * phi * cphi.transpose();
    }

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_vector_monomial_tangential_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> r_phi     = rb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> r_phi_n   = vcross(r_phi, n);
            Matrix<T, Dynamic, 3> f_phi     = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

            cr_rhs.block(0, 0, rbs, cbs) -= qp.weight() * r_phi_n * c_phi.transpose();
            cr_rhs.block(0, offset, rbs, fbs) += qp.weight() * r_phi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    matrix_type oper = mass.ldlt().solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_hho_curl_impl_nedelec(const Mesh&                       msh,
                          const typename Mesh::cell_type&   cl,
                          const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_vector_monomial_basis(msh, cl, celdeg);
    const auto rb = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs = nedelec_tangential_basis_size(facdeg);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto rbs = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);

    const auto      fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stiff = matrix_type::Zero(rbs, rbs);
    matrix_type cr_lhs = matrix_type::Zero(rbs-3, rbs-3);
    matrix_type cr_rhs = matrix_type::Zero(rbs-3, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        const auto cphi = rb.eval_curls2(qp.point());
        stiff += qp.weight() * cphi * cphi.transpose();
    }

    cr_lhs = stiff.block(3, 3, rbs-3, rbs-3);
    cr_rhs.block(0, 0, rbs-3, cbs) = stiff.block(3, 0, rbs-3, cbs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_vector_monomial_nedelec_tangential_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> cphi      = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 3> cphi_n    = vcross(cphi, n).block(3,0,rbs-3,3);
            Matrix<T, Dynamic, 3> f_phi     = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

            cr_rhs.block(0, 0, rbs-3, cbs) -= qp.weight() * cphi_n * c_phi.transpose();
            cr_rhs.block(0, offset, rbs-3, fbs) += qp.weight() * cphi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    matrix_type oper = cr_lhs.ldlt().solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
make_vector_lapl_H0t_oper(const Mesh&                       msh,
                          const typename Mesh::cell_type&   cl,
                          const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb    = make_vector_monomial_basis(msh, cl, celdeg);
    const auto rb    = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs_b = scalar_basis_size(facdeg, Mesh::dimension-1);
    const auto fbs_i = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);
    const auto cbs   = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto rbs   = vector_basis_size(recdeg, Mesh::dimension, Mesh::dimension);

    const auto fcs = faces(msh, cl);
    size_t num_faces_dofs = 0;
    for (auto& fc : fcs)
    {
        if (msh.is_boundary(fc))
            num_faces_dofs += fbs_b;
        else
            num_faces_dofs += fbs_i;
    }

    matrix_type stiff = make_stiffness_matrix(msh, cl, rb);
    matrix_type cr_lhs = matrix_type::Zero(rbs-3, rbs-3);
    matrix_type cr_rhs = matrix_type::Zero(rbs-3, cbs + num_faces_dofs);

    cr_lhs = stiff.block(3, 3, rbs-3, rbs-3);
    cr_rhs.block(0, 0, rbs-3, cbs) = stiff.block(3, 0, rbs-3, cbs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);

        if (msh.is_boundary(fc))
        {
            const auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

            const auto qps_f = integrate(msh, fc, 2*recdeg);
            for (auto& qp : qps_f)
            {
                auto dphi = rb.eval_gradients(qp.point());
                
                Matrix<T, Dynamic, 3> dphi_n   = Matrix<T, Dynamic, 3>::Zero(rbs-3, 3);
                for (size_t di = 3; di < rbs; di++)
                    dphi_n.row(di-3) = (dphi.at(di)*n).transpose();

                Matrix<T, Dynamic, 1> f_phi_tmp = fb.eval_functions(qp.point());
                Matrix<T, Dynamic, 3> f_phi     = f_phi_tmp * n.transpose();
                Matrix<T, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

                cr_rhs.block(0, 0, rbs-3, cbs) -= qp.weight() * dphi_n * c_phi.transpose();
                cr_rhs.block(0, offset, rbs-3, fbs_b) += qp.weight() * dphi_n * f_phi.transpose();
            }
            offset += fbs_b;
        }
        else
        {
            const auto fb = make_vector_monomial_basis(msh, fc, facdeg);

            const auto qps_f = integrate(msh, fc, 2*recdeg);
            for (auto& qp : qps_f)
            {
                auto dphi = rb.eval_gradients(qp.point());
                
                Matrix<T, Dynamic, 3> dphi_n   = Matrix<T, Dynamic, 3>::Zero(rbs-3, 3);
                for (size_t di = 3; di < rbs; di++)
                    dphi_n.row(di-3) = (dphi.at(di)*n).transpose();

                Matrix<T, Dynamic, 3> f_phi     = fb.eval_functions(qp.point());
                Matrix<T, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

                cr_rhs.block(0, 0, rbs-3, cbs) -= qp.weight() * dphi_n * c_phi.transpose();
                cr_rhs.block(0, offset, rbs-3, fbs_i) += qp.weight() * dphi_n * f_phi.transpose();
            }
            offset += fbs_i;
        }
    }

    matrix_type oper = cr_lhs.ldlt().solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh>
dynamic_matrix<typename Mesh::coordinate_type>
make_vector_lapl_H0t_stab(const Mesh&                       msh,
                          const typename Mesh::cell_type&   cl,
                          const hho_degree_info&       cell_infos)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb    = make_vector_monomial_basis(msh, cl, celdeg);
    const auto fbs_b = scalar_basis_size(facdeg, Mesh::dimension-1);
    const auto fbs_i = vector_basis_size(facdeg, Mesh::dimension-1, Mesh::dimension);
    const auto cbs   = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto fcs = faces(msh, cl);
    size_t num_faces_dofs = 0;
    for (auto& fc : fcs)
    {
        if (msh.is_boundary(fc))
            num_faces_dofs += fbs_b;
        else
            num_faces_dofs += fbs_i;
    }

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);

        if (msh.is_boundary(fc))
        {
            const auto fb = make_scalar_monomial_basis(msh, fc, facdeg);
            const auto hf = diameter(msh, fc);

            matrix_type mass = matrix_type::Zero(fbs_b, fbs_b);
            matrix_type trace = matrix_type::Zero(fbs_b, cbs);
            matrix_type rhs = matrix_type::Zero(fbs_b, cbs + num_faces_dofs);

            rhs.block(0, offset, fbs_b, fbs_b) = matrix_type::Identity(fbs_b, fbs_b);

            const auto qps_f = integrate(msh, fc, celdeg+facdeg);
            for (auto& qp : qps_f)
            {
                Matrix<T, Dynamic, 1> f_phi_tmp = fb.eval_functions(qp.point());
                Matrix<T, Dynamic, 3> f_phi     = f_phi_tmp * n.transpose();
                Matrix<T, Dynamic, 3> c_phi = cb.eval_functions(qp.point());

                mass += qp.weight() * f_phi * f_phi.transpose();
                trace += qp.weight() * f_phi * c_phi.transpose();
            }

            rhs.block(0,0,fbs_b,cbs) = -mass.ldlt().solve(trace);
            stab += rhs.transpose() * mass * rhs / hf;

            offset += fbs_b;
        }
        else
        {
            const auto fb = make_vector_monomial_basis(msh, fc, facdeg);
            const auto hf = diameter(msh, fc);

            matrix_type mass = matrix_type::Zero(fbs_i, fbs_i);
            matrix_type trace = matrix_type::Zero(fbs_i, cbs);
            matrix_type rhs = matrix_type::Zero(fbs_i, cbs + num_faces_dofs);

            rhs.block(0, offset, fbs_i, fbs_i) = matrix_type::Identity(fbs_i, fbs_i);

            const auto qps_f = integrate(msh, fc, celdeg+facdeg);
            for (auto& qp : qps_f)
            {
                Matrix<T, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
                Matrix<T, Dynamic, 3> c_phi = cb.eval_functions(qp.point());

                mass += qp.weight() * f_phi * f_phi.transpose();
                trace += qp.weight() * f_phi * c_phi.transpose();
            }

            rhs.block(0,0,fbs_i,cbs) = -mass.ldlt().solve(trace);
            stab += rhs.transpose() * mass * rhs / hf;

            offset += fbs_i;
        }
    }

    return stab;
}


template<typename T>
std::pair<dynamic_matrix<T>, dynamic_vector<T>>
static_condensation(const dynamic_matrix<T>& lhs, const dynamic_vector<T>& rhs,
                    const size_t tf_bnd)
{
    auto t_size = tf_bnd;
    auto f_size = lhs.rows() - tf_bnd;

    dynamic_matrix<T> LTT = lhs.block(     0,      0, t_size, t_size);
    dynamic_matrix<T> LTF = lhs.block(     0, tf_bnd, t_size, f_size);
    dynamic_matrix<T> LFT = lhs.block(tf_bnd,      0, f_size, t_size);
    dynamic_matrix<T> LFF = lhs.block(tf_bnd, tf_bnd, f_size, f_size);

    dynamic_vector<T> bT = rhs.segment(     0, t_size);
    dynamic_vector<T> bF = rhs.segment(tf_bnd, f_size);

    LDLT<dynamic_matrix<T>> ldlt_LTT(LTT);
    if (ldlt_LTT.info() != Eigen::Success)
        throw std::invalid_argument("Can't factorize matrix for static condensation");

    dynamic_matrix<T> LC = LFF - LFT*ldlt_LTT.solve(LTF);
    dynamic_vector<T> bC = bF - LFT*ldlt_LTT.solve(bT);

    return std::make_pair(LC, bC);
}

template<typename T>
dynamic_vector<T>
static_decondensation(const dynamic_matrix<T>& lhs, const dynamic_vector<T>& rhs,
                      const dynamic_vector<T>& uF)
{
    auto t_size = lhs.rows() - uF.rows();
    auto tf_bnd = t_size;
    auto f_size = uF.rows();

    dynamic_matrix<T> LTT = lhs.block(0,      0, t_size, t_size);
    dynamic_matrix<T> LTF = lhs.block(0, tf_bnd, t_size, f_size);

    dynamic_vector<T> bT = rhs.segment(0, t_size);

    dynamic_vector<T> ret = dynamic_vector<T>::Zero(lhs.rows());
    LDLT<dynamic_matrix<T>> ldlt_LTT(LTT);
    if (ldlt_LTT.info() != Eigen::Success)
        throw std::invalid_argument("Can't factorize matrix for static condensation");

    ret.segment(     0, t_size) = ldlt_LTT.solve(bT - LTF*uF);
    ret.segment(tf_bnd, f_size) = uF;

    return ret;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
dynamic_matrix<typename Mesh<T,2,Storage>::coordinate_type>
curl_hdg_stabilization(const Mesh<T,2,Storage>&                     msh,
                       const typename Mesh<T,2,Storage>::cell_type& cl,
                       const hho_degree_info&                       cell_infos)
{
    /* This is actually the plain HDG stabilization */
    using mesh_type         = Mesh<T,2,Storage>;
    using scalar_type       = typename mesh_type::coordinate_type;
    using matrix_type       = Matrix<T, Dynamic, Dynamic>;
    using vector_type       = Matrix<T, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;

    static const size_t DIM = 2;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb  = make_scalar_monomial_basis(msh, cl, celdeg);
    const auto fbs = scalar_basis_size(facdeg, DIM - 1);
    const auto cbs = scalar_basis_size(celdeg, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_scalar_monomial_basis(msh, fc, facdeg);
        const auto hf = diameter(msh, fc);

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);
        matrix_type rhs = matrix_type::Zero(fbs, cbs + num_faces_dofs);

        rhs.block(0, offset, fbs, fbs) = matrix_type::Identity(fbs, fbs);

        const auto qps_f = integrate(msh, fc, 2*std::max(celdeg,facdeg));
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> f_phi = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> c_phi = cb.eval_functions(qp.point());

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        rhs.block(0,0,fbs,cbs) = -mass.ldlt().solve(trace);
        stab += rhs.transpose() * mass * rhs / hf;

        offset += fbs;
    }

    return stab;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
dynamic_matrix<typename Mesh<T,3,Storage>::coordinate_type>
curl_hdg_stabilization(const Mesh<T,3,Storage>&                     msh,
                       const typename Mesh<T,3,Storage>::cell_type& cl,
                       const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<T,3,Storage>;
    using scalar_type       = typename mesh_type::coordinate_type;
    using matrix_type       = Matrix<T, Dynamic, Dynamic>;
    using vector_type       = Matrix<T, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;

    static const size_t DIM = 3;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb  = make_vector_monomial_basis(msh, cl, celdeg);
    const auto fbs = vector_basis_size(facdeg, DIM - 1, DIM - 1);
    const auto cbs = vector_basis_size(celdeg, DIM, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_tangential_basis(msh, fc, facdeg);
        const auto hf = diameter(msh, fc);

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);
        matrix_type rhs = matrix_type::Zero(fbs, cbs + num_faces_dofs);

        rhs.block(0, offset, fbs, fbs) = matrix_type::Identity(fbs, fbs);

        const auto qps_f = integrate(msh, fc, 2*std::max(celdeg,facdeg));
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi_tmp = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi = vcross(n, vcross(c_phi_tmp, n));

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        rhs.block(0,0,fbs,cbs) = -mass.ldlt().solve(trace);
        stab += rhs.transpose() * mass * rhs / hf;

        offset += fbs;
    }

    return stab;
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
dynamic_matrix<typename Mesh<T,3,Storage>::coordinate_type>
curl_hho_stabilization(const Mesh<T,3,Storage>&                     msh,
                       const typename Mesh<T,3,Storage>::cell_type& cl,
                       const Matrix<T,Dynamic,Dynamic>              recop,
                       const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<T,3,Storage>;
    using scalar_type       = typename mesh_type::coordinate_type;
    using matrix_type       = Matrix<T, Dynamic, Dynamic>;
    using vector_type       = Matrix<T, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;

    static const size_t DIM = 3;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb  = make_vector_monomial_basis(msh, cl, celdeg);
    const auto rb  = make_vector_monomial_basis(msh, cl, recdeg);
    const auto fbs = vector_basis_size(facdeg, DIM - 1, DIM - 1);
    const auto cbs = vector_basis_size(celdeg, DIM, DIM);
    const auto rbs = vector_basis_size(recdeg, DIM, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);
    
    matrix_type Rtrace = matrix_type::Zero(cbs, rbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto c_phi = cb.eval_functions(qp.point());
        auto r_phi = rb.eval_functions(qp.point());

        Rtrace += qp.weight() * c_phi * r_phi.transpose();
    }

    matrix_type Rtrace1 = Rtrace.block(0,3,cbs,rbs-3);
    matrix_type Cmass = Rtrace.block(0,0,cbs,cbs);

    /* project reconstruction on cell */
    matrix_type RprojC = Cmass.ldlt().solve(Rtrace1 * recop);
    /* subtract vT */
    RprojC.block(0,0,cbs,cbs) -= matrix_type::Identity(cbs,cbs);
    /* RprojC is now π_T(R(v)) - vT*/

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_tangential_basis(msh, fc, facdeg);
        const auto hf = diameter(msh, fc);

        matrix_type Fmass = matrix_type::Zero(fbs, fbs);
        matrix_type CtraceF = matrix_type::Zero(fbs, cbs);
        matrix_type RtraceF = matrix_type::Zero(fbs, rbs-3);
        matrix_type rhs = matrix_type::Zero(fbs, cbs + num_faces_dofs);

        rhs.block(0, offset, fbs, fbs) = matrix_type::Identity(fbs, fbs);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi_tmp = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi = vcross(n, vcross(c_phi_tmp, n));
            Matrix<T, Dynamic, 3> r_phi = rb.eval_functions(qp.point()).block(3,0,rbs-3,3);

            Fmass += qp.weight() * f_phi * f_phi.transpose();
            CtraceF += qp.weight() * f_phi * c_phi.transpose();
            RtraceF += qp.weight() * f_phi * r_phi.transpose();
        }

        LDLT<matrix_type> Fmass_ldlt(Fmass);

        /* project reconstruction on face */
        matrix_type RprojF = -Fmass_ldlt.solve(RtraceF * recop);
        /* add vF */
        RprojF.block(0,offset,fbs,fbs) += matrix_type::Identity(fbs, fbs);
        /* RprojF is now vF - π_F(R(v)) */

        matrix_type rf = Fmass_ldlt.solve(CtraceF * RprojC);
        /* rf is now π_F(π_T(R(v)) - vT)*/

        rhs = RprojF + rf;
        stab += rhs.transpose() * Fmass * rhs / hf;

        offset += fbs;
    }

    return stab;
}

template<typename Mesh, typename Element, typename Basis>
Matrix<typename Basis::scalar_type, Dynamic, Dynamic>
make_curl_curl_matrix(const Mesh& msh, const Element& elem, const Basis& basis, size_t di = 0)
{
    const auto degree     = basis.degree();
    const auto basis_size = basis.size();

    using T = typename Basis::scalar_type;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(basis_size, basis_size);

    const auto qps = integrate(msh, elem, 2*(degree+di));
    for (auto& qp : qps)
    {
        const auto cphi    = basis.eval_curls2(qp.point());
        ret += qp.weight() * cphi * cphi.transpose();
    }

    return ret;
}


}



