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

template<typename T, typename U>
Matrix<T, Dynamic, 3>
vcross(const static_vector<U, 3>& normal, const Matrix<T, Dynamic, 3>& field)
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

template<typename T, typename U>
Matrix<T, Dynamic, 3>
vcross(const Matrix<T, Dynamic, 3>& field, const static_vector<U, 3>& normal)
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

template<typename ScalT = double, typename Mesh>
dynamic_matrix<ScalT>
make_vector_mass_oper(const Mesh&                       msh,
                      const typename Mesh::cell_type&   cl,
                      const hho_degree_info&       cell_infos)
{
    using mesh_type = Mesh;
    using coord_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    using cell_type = typename Mesh::cell_type;

    typedef Matrix<scalar_type, Dynamic, Dynamic> matrix_type;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, celdeg);
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

template<typename ScalT = double, typename Mesh>
dynamic_matrix<ScalT>
make_impedance_term(const Mesh& msh, const typename Mesh::face_type& fc,
    size_t facdeg)
{
    using mesh_type = Mesh;
    using coord_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;

    using matrix_type = Matrix<scalar_type, Dynamic, Dynamic>;
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension-1);

    matrix_type Y = matrix_type::Zero(fbs, fbs);

    const auto fb = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);
    const auto qps = integrate(msh, fc, 2*facdeg);
    for (auto& qp : qps)
    {
        Matrix<scalar_type, Dynamic, 3> fphi = fb.eval_functions(qp.point());
        Y += qp.weight() * fphi * fphi.transpose();
    }

    return Y;
}

template<typename ScalT = double, typename Mesh, typename Function>
std::pair<dynamic_matrix<ScalT>, dynamic_vector<ScalT>>
make_plane_wave_term(const Mesh& msh, const typename Mesh::face_type& fc,
    size_t facdeg, Function f)
{
    using mesh_type = Mesh;
    using coord_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;

    using matrix_type = Matrix<scalar_type, Dynamic, Dynamic>;
    using vector_type = Matrix<scalar_type, Dynamic, 1>;
    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension-1);

    matrix_type Y = matrix_type::Zero(fbs, fbs);
    vector_type y = vector_type::Zero(fbs);

    const auto fb = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);
    const auto qps = integrate(msh, fc, 2*facdeg);
    for (auto& qp : qps)
    {
        Matrix<scalar_type, Dynamic, 3> fphi = fb.eval_functions(qp.point());
        Y += qp.weight() * fphi * fphi.transpose();
        y += qp.weight() * fphi * f(qp.point());
    }

    return std::make_pair(Y, y);
}

template<typename ScalT = double, typename Mesh>
std::pair<dynamic_matrix<ScalT>, dynamic_vector<ScalT>>
make_impedance_term(const Mesh& msh, const typename Mesh::cell_type& cl,
    const hho_degree_info& hdi, size_t bndid, ScalT z)
{
    using mesh_type = Mesh;
    using coord_type = typename Mesh::coordinate_type;
    using scalar_type = ScalT;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;

    typedef Matrix<scalar_type, Dynamic, Dynamic> matrix_type;
    typedef Matrix<scalar_type, Dynamic, 1> vector_type;

    const auto facdeg = hdi.face_degree();
    const auto celdeg = hdi.cell_degree();
    const auto recdeg = hdi.reconstruction_degree();

    const auto fbs = vector_basis_size(facdeg, Mesh::dimension - 1, Mesh::dimension-1);
    const auto cbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);
    const auto rbs = vector_basis_size(celdeg, Mesh::dimension, Mesh::dimension);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    const auto rb = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, recdeg);

    matrix_type Y = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);
    vector_type yr = vector_type::Zero(cbs + num_faces_dofs);

    auto offset = cbs;

    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];

        if ( msh.is_boundary(fc) and msh.boundary_id(fc) == bndid )
        {
            const auto fb = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);
            const auto qps_f = integrate(msh, fc, 2*facdeg);
            for (auto& qp : qps_f)
            {
                Matrix<scalar_type, Dynamic, 3> fphi = fb.eval_functions(qp.point());
                Y.block(offset, offset, fbs, fbs) += qp.weight() * z * fphi * fphi.transpose();
                yr.segment(offset, fbs) += qp.weight() * z * fphi.col(2);
            }
        }

        offset += fbs;
    }

    return std::make_pair(Y, yr);
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

template<typename ScalT = double, template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
matrixpair<ScalT>
curl_reconstruction(const Mesh<CoordT,3,Storage>&                     msh,
                    const typename Mesh<CoordT,3,Storage>::cell_type& cl,
                    const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<CoordT,3,Storage>;
    using scalar_type       = ScalT;
    using matrix_type       = Matrix<scalar_type, Dynamic, Dynamic>;
    using vector_type       = Matrix<scalar_type, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;
    using cell_type         = typename mesh_type::cell_type;
    using face_type         = typename mesh_type::face_type;

    static const size_t DIM = 3;
    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();
    const auto recdeg = cell_infos.reconstruction_degree();

    const auto cb = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, celdeg);
    const auto rb = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, recdeg);
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
        const auto fb     = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<scalar_type, Dynamic, 3> cphi      = rb.eval_curls2(qp.point());
            Matrix<scalar_type, Dynamic, 3> cphi_n    = vcross(cphi, n).block(3,0,rbs-3,3);
            Matrix<scalar_type, Dynamic, 3> f_phi     = fb.eval_functions(qp.point());
            Matrix<scalar_type, Dynamic, 3> c_phi     = cb.eval_functions(qp.point());

            cr_rhs.block(0, 0, rbs-3, cbs) -= qp.weight() * cphi_n * c_phi.transpose();
            cr_rhs.block(0, offset, rbs-3, fbs) += qp.weight() * cphi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    FullPivLU<matrix_type> ldlt_lhs(cr_lhs);
    //if (ldlt_lhs.info() != Eigen::Success)
    //    throw std::invalid_argument("Can't factorize matrix for curl reconstruction");

    matrix_type oper = ldlt_lhs.solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename Mesh, typename Function>
dynamic_vector<std::complex<double>>
maxwell_prj(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi, const Function& f)
{
    using T = std::complex<double>;

    const auto fbs = vector_basis_size(hdi.face_degree(), Mesh::dimension - 1, Mesh::dimension-1);
    const auto cbs = vector_basis_size(hdi.cell_degree(), Mesh::dimension, Mesh::dimension);
    
    const auto fcs = faces(msh, cl);

    dynamic_vector<T> ret = dynamic_vector<T>::Zero(cbs+fcs.size() * fbs);

    dynamic_matrix<T> cmat = dynamic_matrix<T>::Zero(cbs, cbs);
    dynamic_matrix<T> crhs = dynamic_vector<T>::Zero(cbs);

    const auto cb = make_vector_monomial_basis(msh, cl, hdi.cell_degree());
    const auto qps = integrate(msh, cl, 2*hdi.cell_degree());
    for (const auto& qp : qps)
    {
        auto phi = cb.eval_functions(qp.point());
        cmat += qp.weight() * phi * phi.transpose();
        crhs += qp.weight() * phi * f(qp.point());
    }

    ret.segment(0, cbs) = cmat.ldlt().solve(crhs);

    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto& fc = fcs[i];
        const auto fb = make_vector_monomial_tangential_basis(msh, fc, hdi.face_degree());
        const auto qps_f = integrate(msh, cl, 2*hdi.face_degree());
        dynamic_matrix<T> fmat = dynamic_matrix<T>::Zero(fbs, fbs);
        dynamic_matrix<T> frhs = dynamic_vector<T>::Zero(fbs);
        for (const auto& qp : qps_f)
        {
            auto phi = fb.eval_functions(qp.point());
            fmat += qp.weight() * phi * phi.transpose();
            frhs += qp.weight() * phi * f(qp.point());
        }

        ret.segment(cbs+i*fbs, fbs) = fmat.ldlt().solve(frhs);
    }

    return ret;
}

template<typename Mesh>
std::pair<dynamic_matrix<typename Mesh::coordinate_type>, dynamic_matrix<typename Mesh::coordinate_type>>
curl_reconstruction_pk(const Mesh&                       msh,
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
        const auto r_phi = rb.eval_functions(qp.point());
        mass += qp.weight() * r_phi * r_phi.transpose();

        const auto r_cphi = rb.eval_curls2(qp.point());
        const auto c_phi = cb.eval_functions(qp.point());
        cr_rhs.block(0, 0, rbs, cbs) += qp.weight() * r_cphi * c_phi.transpose();
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

            cr_rhs.block(0, offset, rbs, fbs) += qp.weight() * r_phi_n * f_phi.transpose();
        }

        offset += fbs;
    }

    LDLT<matrix_type> ldlt_lhs(mass);
    if (ldlt_lhs.info() != Eigen::Success)
        throw std::invalid_argument("Can't factorize matrix for curl reconstruction");

    matrix_type oper = ldlt_lhs.solve(cr_rhs);
    matrix_type data = cr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
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

template<typename ScalT=double, template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
dynamic_matrix<ScalT>
curl_hdg_stabilization(const Mesh<CoordT,3,Storage>&                     msh,
                       const typename Mesh<CoordT,3,Storage>::cell_type& cl,
                       const hho_degree_info&                       cell_infos)
{
    using mesh_type         = Mesh<CoordT,3,Storage>;
    using scalar_type       = ScalT;
    using matrix_type       = Matrix<scalar_type, Dynamic, Dynamic>;
    using vector_type       = Matrix<scalar_type, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;
    using cell_type         = typename mesh_type::cell_type;
    using face_type         = typename mesh_type::face_type;

    static const size_t DIM = 3;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb  = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, celdeg);
    const auto fbs = vector_basis_size(facdeg, DIM - 1, DIM - 1);
    const auto cbs = vector_basis_size(celdeg, DIM, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    auto ht = diameter(msh, cl);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);
        const auto hf = diameter(msh, fc);

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);
        matrix_type rhs = matrix_type::Zero(fbs, cbs + num_faces_dofs);

        rhs.block(0, offset, fbs, fbs) = matrix_type::Identity(fbs, fbs);

        const auto qps_f = integrate(msh, fc, 2*std::max(celdeg,facdeg));
        for (auto& qp : qps_f)
        {
            Matrix<scalar_type, Dynamic, 3> f_phi = fb.eval_functions(qp.point());
            Matrix<scalar_type, Dynamic, 3> c_phi_tmp = cb.eval_functions(qp.point());
            Matrix<scalar_type, Dynamic, 3> c_phi = vcross(n, vcross(c_phi_tmp, n));

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * c_phi.transpose();
        }

        rhs.block(0,0,fbs,cbs) = -mass.ldlt().solve(trace);
        stab += rhs.transpose() * mass * rhs / ht;

        offset += fbs;
    }

    return stab;
}

/* This essentialy enforces an impedance condition between the function on
 * the cell and the function on the face:
 *
 *                      (∇×Eᵗ)×n - iωμY(n×Eᶠ)×n = 0.
 *
 * Looks like it does not work however. Still need to determine if it is
 * actually supposed to not work or if it is just PEBKAC.
 */
template<typename ScalT=double, template<typename, size_t, typename> class Mesh, typename CoordT, typename Storage>
dynamic_matrix<ScalT>
curl_Z_stabilization(const Mesh<CoordT,3,Storage>&                     msh,
                       const typename Mesh<CoordT,3,Storage>::cell_type& cl,
                       const hho_degree_info&                       cell_infos,
                       ScalT mur, ScalT jwmu0, ScalT Z)
{
    using mesh_type         = Mesh<CoordT,3,Storage>;
    using scalar_type       = ScalT;
    using matrix_type       = Matrix<scalar_type, Dynamic, Dynamic>;
    using vector_type       = Matrix<scalar_type, Dynamic, 1>;
    using point_type        = typename mesh_type::point_type;
    using cell_type         = typename mesh_type::cell_type;
    using face_type         = typename mesh_type::face_type;

    static const size_t DIM = 3;

    const auto facdeg = cell_infos.face_degree();
    const auto celdeg = cell_infos.cell_degree();

    const auto cb  = make_vector_monomial_basis<mesh_type, cell_type, scalar_type>(msh, cl, celdeg);
    const auto fbs = vector_basis_size(facdeg, DIM - 1, DIM - 1);
    const auto cbs = vector_basis_size(celdeg, DIM, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    matrix_type stab = matrix_type::Zero(cbs + num_faces_dofs, cbs + num_faces_dofs);

    auto ht = diameter(msh, cl);

    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = normal(msh, cl, fc);
        const auto fb = make_vector_monomial_tangential_basis<mesh_type, face_type, scalar_type>(msh, fc, facdeg);
        const auto hf = diameter(msh, fc);

        matrix_type mass = matrix_type::Zero(fbs, fbs);
        matrix_type trace = matrix_type::Zero(fbs, cbs);
        matrix_type rhs = matrix_type::Zero(fbs, cbs + num_faces_dofs);

        /* Evaluate -2iωμ0Y(n×Eᶠ)×n */
        rhs.block(0, offset, fbs, fbs) = (-/*2.**/jwmu0/Z)*matrix_type::Identity(fbs, fbs);

        const auto qps_f = integrate(msh, fc, 2*std::max(celdeg,facdeg));
        for (auto& qp : qps_f)
        {
            /* Evaluate iωμ0Y(n×Eᵗ)×n */
            Matrix<scalar_type, Dynamic, 3> phi_tmp = cb.eval_functions(qp.point());
            Matrix<scalar_type, Dynamic, 3> n_x_phi_x_n = disk::vcross(n, disk::vcross(phi_tmp, n));

            /* Evaluate (1/mur)*(∇×Eᵗ)×n */
            //Matrix<scalar_type, Dynamic, 3> cphi_tmp = cb.eval_curls2(qp.point());
            //Matrix<scalar_type, Dynamic, 3> cphi_x_n = (1./mur)*disk::vcross(cphi_tmp, n);

            Matrix<scalar_type, Dynamic, 3> f_phi = fb.eval_functions(qp.point());

            mass += qp.weight() * f_phi * f_phi.transpose();
            trace += qp.weight() * f_phi * n_x_phi_x_n.transpose();
            //trace += qp.weight() * f_phi * cphi_x_n.transpose();
        }

        rhs.block(0,0,fbs,cbs) = (jwmu0/Z)*mass.ldlt().solve(trace);
        stab += rhs.transpose() * mass * rhs;

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








template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
matrixpair<typename Mesh<T,2,Storage>::coordinate_type>
curl_reconstruction_div(const Mesh<T,2,Storage>&                     msh,
                    const typename Mesh<T,2,Storage>::cell_type& cl,
                    const hho_degree_info&                       cell_infos)
{
}

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
matrixpair<typename Mesh<T,3,Storage>::coordinate_type>
curl_reconstruction_div(const Mesh<T,3,Storage>&                     msh,
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

    const auto db = make_scalar_monomial_basis(msh, cl, recdeg);
    const auto dbs = scalar_basis_size(recdeg, DIM);

    const auto fcs = faces(msh, cl);
    const auto num_faces_dofs = fcs.size() * fbs;

    const auto nunk = rbs+dbs+3;

    matrix_type lhs = matrix_type::Zero(nunk, nunk);
    matrix_type rhs = matrix_type::Zero(nunk, cbs + num_faces_dofs);

    const auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        const auto      Rphi = rb.eval_functions(qp.point());
        const auto curl_Rphi = rb.eval_curls2(qp.point());
        const auto  div_Rphi = rb.eval_divergences(qp.point());
        
        const auto      Dphi = db.eval_functions(qp.point());
        const auto grad_Dphi = db.eval_gradients(qp.point());
        
        const Matrix<T, Dynamic, 3>      Cphi = Rphi.block(0, 0, cbs, 3);
        const Matrix<T, Dynamic, 3> curl_Cphi = curl_Rphi.block(0, 0, cbs, 3);
        const Matrix<T, Dynamic, 1>  div_Cphi = div_Rphi.segment(0, cb.size());

        lhs.block(0,   0, rbs, rbs) += qp.weight() * curl_Rphi * curl_Rphi.transpose();
        rhs.block(0,   0, rbs, cbs) += qp.weight() * curl_Rphi * curl_Cphi.transpose();

        lhs.block(0, rbs, rbs, dbs) += qp.weight() * Rphi * grad_Dphi.transpose();
        lhs.block(rbs, 0, dbs, rbs) += qp.weight() * grad_Dphi * Rphi.transpose();
        rhs.block(rbs, 0, dbs, cbs) -= qp.weight() * Dphi * div_Cphi.transpose();

        lhs.block(rbs+dbs, 0, 3, rbs) += qp.weight() * Rphi.transpose();
        rhs.block(rbs+dbs, 0, 3, cbs) += qp.weight() * Cphi.transpose();
    }

    lhs.block(0, rbs+dbs, rbs, 3) = lhs.block(rbs+dbs, 0, 3, rbs).transpose();


    size_t offset = cbs;
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n      = normal(msh, cl, fc);
        const auto fb     = make_vector_monomial_tangential_basis(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, 2*recdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 3> phi           = rb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> phi_n         = vcross(phi, n);
            Matrix<T, Dynamic, 3> curl_phi      = rb.eval_curls2(qp.point());
            Matrix<T, Dynamic, 3> curl_phi_n    = vcross(curl_phi, n);
            Matrix<T, Dynamic, 3> f_phi         = fb.eval_functions(qp.point());
            Matrix<T, Dynamic, 3> c_phi         = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> d_phi         = db.eval_functions(qp.point());

            rhs.block(0,        0, rbs, cbs) -= qp.weight() * curl_phi_n * c_phi.transpose();
            rhs.block(0,   offset, rbs, fbs) += qp.weight() * curl_phi_n * f_phi.transpose();
            rhs.block(rbs,      0, dbs, cbs) += qp.weight() * d_phi * (c_phi * n).transpose();
        }

        offset += fbs;
    }

    FullPivLU<matrix_type> ldlt_lhs(lhs);
    //if (ldlt_lhs.info() != Eigen::Success)
    //    throw std::invalid_argument("Can't factorize matrix for curl reconstruction");

    matrix_type oper = ldlt_lhs.solve(rhs);
    matrix_type roper = oper.block(0,0,rbs,oper.cols());
    matrix_type data = rhs.block(0,0,rbs,rhs.cols()).transpose() * roper;

    return std::make_pair(roper, data);
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

    matrix_type real_recop = recop/*.block(3, 0, rbs-3, recop.cols())*/;

    /* project reconstruction on cell */
    matrix_type RprojC = Cmass.ldlt().solve(Rtrace1 * real_recop);
//RprojC -= RprojC; /*DEBUG,leave commented*/
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
            Matrix<T, Dynamic, 3> r_phi_tmp = rb.eval_functions(qp.point()).block(3,0,rbs-3,3);
            Matrix<T, Dynamic, 3> r_phi = vcross(n, vcross(r_phi_tmp, n));

            Fmass += qp.weight() * f_phi * f_phi.transpose();
            CtraceF += qp.weight() * f_phi * c_phi.transpose();
            RtraceF += qp.weight() * f_phi * r_phi.transpose();
        }

        LDLT<matrix_type> Fmass_ldlt(Fmass);

        /* project reconstruction on face */
        matrix_type RprojF = -Fmass_ldlt.solve(RtraceF * real_recop);
//RprojF -= RprojF; /*DEBUG,leave commented*/
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













}



