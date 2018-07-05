template <typename Mesh>
Matrix< typename Mesh::coordinate_type, Dynamic, Dynamic>
make_contact_nitsche(const Mesh& msh, const typename Mesh::cell_type& cl,
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
            auto qps = integrate(msh, fc, 2* hdi.face_degree());

            for (auto& qp : qps)
            {
                Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                Matrix<T, Dynamic, 1>  sigma_n = rec.transpose() * (c_dphi * n);

                ret +=  sigma_n * sigma_n.transpose();
            }
        }
    }

    return ret * (theta /gamma_N);
}

template <typename Mesh, typename T>
Matrix< T, Dynamic, Dynamic>
make_contact_heaviside_other(const Mesh& msh, const typename Mesh::cell_type& cl,
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
make_contact_heaviside(const Mesh& msh, const typename Mesh::cell_type& cl,
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
                vector_type  sigma_n   = rec.transpose() * (c_dphi * n);

                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_dphi_tmp.block(0, 0, cbs, 1);
                vector_type  u_cell  = uloc.block(0, 0, cbs, 1);

                // Heaviside(-Pn(u) = -D*recu*n + gN*u_T)
                T gamma_u  = gamma_N * c_phi.dot(u_cell);
                T sigmau_n = sigma_n.dot(uloc);

                if (sigmau_n - gamma_u  <= 0.)
                {
                    // (theta * grad * rec * v * n - gamma_N* v_T)
                    vector_type t_sigman_g_v = theta * sigma_n;
                    t_sigman_g_v.block(0, 0,cbs, 1) -=  gamma_N * c_phi;

                    // (grad * rec * u * n - gamma_N* u_T)
                    vector_type sigman_g_u = sigma_n;
                    sigman_g_u.block(0, 0,cbs, 1) -=  gamma_N * c_phi;

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
make_contact_heaviside_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
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
make_contact_heaviside_trace(const Mesh& msh, const typename Mesh::cell_type& cl,
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

            auto qps = integrate(msh, fc, 2* hdi.face_degree());
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

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_contact_negative(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> rec,
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

            for (auto& qp : qps)
            {

                Matrix<T, Dynamic,  DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
                Matrix<T, Dynamic,  DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
                vector_type  c_phi_temp   = cb.eval_functions(qp.point());
                vector_type  c_phi = c_dphi_tmp.block(0, 0, cbs, 1);

                vector_type  sigma_n = rec.transpose() * (c_dphi * n);

                // (grad * rec * v * n - gamma_N * v_T)
                vector_type t_sigman_gamma_v =  theta * sigma_n;
                t_sigman_gamma_v.block(0, 0, cbs, 1) -= gamma_N * c_phi;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* u_T]_
                vector_type  u_cell  = uloc.block(0, 0, cbs, 1);
                T gamma_u  = gamma_N * c_phi.dot(u_cell);
                T sigmau_n = sigma_n.dot(uloc);
                T negative = std::min(sigmau_n - gamma_u, 0.);

                rhs += qp.weight() * negative * t_sigman_gamma_v;
            }
        }
        fc_count++;
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_contact_negative_faces(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> rec,
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
                vector_type  u_face  = uloc.block(face_ofs, 0, fbs, 1);
                vector_type  sigma_n = rec.transpose() * (c_dphi * n);

                // (grad * rec * v * n - gamma_N * v_F)
                vector_type t_sigman_gamma_v =  theta * sigma_n;
                t_sigman_gamma_v.block(face_ofs, 0, fbs, 1) -= gamma_N * f_phi;

                // [Pn(u)]_ = [grad * rec * u * n - gamma_N* u_F]_
                T gamma_u  = gamma_N * f_phi.dot(u_face);
                T sigmau_n = sigma_n.dot(uloc);
                T negative = std::min(sigmau_n - gamma_u, 0.);

                rhs += qp.weight() * negative * t_sigman_gamma_v;
            }
        }
        fc_count++;
    }
    return rhs * (1./gamma_N);
}

template <typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_contact_negative_trace(const Mesh& msh, const typename Mesh::cell_type& cl,
            const hho_degree_info& hdi,
            const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic> rec,
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
            auto qps = integrate (msh, fc, 2*facdeg);
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

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_contact_hho_scalar_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl,
                          const hho_degree_info& di,
                          const disk::mechanics::BoundaryConditionsScalar<Mesh>& bnd)
{
    using T = typename Mesh::coordinate_type;
    const size_t DIM = Mesh::dimension;

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    auto fbs = scalar_basis_size(facdeg, Mesh::dimension-1);

    auto num_faces = howmany_faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, cbs + num_faces*fbs);

    auto qps = integrate(msh, cl, 2*recdeg);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.point());
        stiff += qp.weight() * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];

        auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
        if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
        const auto face_id=eid.second;

        if (bnd.is_contact_face(face_id))
            continue;

        auto n = normal(msh, cl, fc);
        auto fb = make_scalar_monomial_basis(msh, fc, facdeg);

        auto qps_f = integrate(msh, fc, 2*facdeg);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> c_phi_tmp = cb.eval_functions(qp.point());
            Matrix<T, Dynamic, 1> c_phi = c_phi_tmp.head(cbs);
            Matrix<T, Dynamic, DIM> c_dphi_tmp = cb.eval_gradients(qp.point());
            Matrix<T, Dynamic, DIM> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, DIM);
            Matrix<T, Dynamic, 1> f_phi = fb.eval_functions(qp.point());
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.weight() * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.weight() * (c_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}
