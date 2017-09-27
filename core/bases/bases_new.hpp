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

template<typename Mesh, typename T>
class abstract_basis
{
public:
    abstract_basis()
    {}
};

template<typename Basis>
class gradient_op
{
    const Basis&    m_target_basis;

public:
    typedef typename Basis::gradient_type       gradient_type;

    gradient_op(const Basis& basis)
        m_target_basis(basis)
    {}

    gradient_type operator()(const point_type& pt)
    {
        return m_target_basis->evaluate_gradients(pt);
    }

};

template<typename Basis>
auto grad(const Basis& basis)
{
    return gradient_op<Basis>(basis);
}



template<typename Mesh, typename T>
class scaled_monomial_scalar_basis
{
public:
    typedef Mesh                    mesh_type;
    typedef mesh_type::point_type   point_type;
    typedef T                       value_type;
    typedef T                       function_type;

    scaled_monomial_scalar_basis()
    {}

    function_type operator()(const point_type& pt)
    {
        return evaluate_functions(pt);
    }
};


template<typename Basis>
class scaled_monomial_scalar_basis_view
{

};


template<typename Mesh, typename T>
class hho_global_space
{

public:
    hho_global_space() {}
    hho_global_space(size_t degree) {}
    hho_global_space(size_t cell_degree, size_t face_degree) {}

    get_local_space(T) {}
};


template<typename Mesh, typename T>
class hho_basis
{
    typedef T                               value_type;
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::cell_type   cell_type;
    typedef typename mesh_type::face_type   face_type;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>  cell_basis_type;
    typedef scaled_monomial_scalar_basis<mesh_type, face_type>  face_basis_type;

    mesh_type       m_msh;
    cell_type       m_cl;

public:
    hho_basis(const mesh_type& msh, const cell_type& cl)
        : m_msh(msh), m_cl(cl)
    {}

    hho_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
        : m_msh(msh), m_cl(cl)
    {}

    hho_basis(const mesh_type& msh, const cell_type& cl, const std::vector<size_t>& degrees)
        : m_msh(msh), m_cl(cl)
    {}

    cell_basis_type T() const
    {
        return cell_basis;
    }

    face_basis_type F(size_t which) const
    {
        return face_basis;
    }

};
