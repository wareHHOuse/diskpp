/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2018, 2019                 nicolas.pignet@enpc.fr
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

#include <tuple>
#include <vector>

#include "diskpp/bases/bases.hpp"
#include "diskpp/common/eigen.hpp"
#include "diskpp/mesh/point.hpp"

namespace disk
{

enum DirichletType : size_t
{
    DIRICHLET = 0,
    CLAMPED   = 1,
    DX        = 2,
    DY        = 3,
    DZ        = 4,
    DXDY      = 5,
    DXDZ      = 6,
    DYDZ      = 7,
    NOTHING   = 8
};

enum NeumannType : size_t
{
    NEUMANN = 9,
    FREE    = 10,
};
enum RobinType : size_t
{
    ROBIN    = 11,
    WHATEVER = 12,
};

enum ContactType : size_t
{
    SIGNORINI      = 13,
    SIGNORINI_FACE = 14,
    SIGNORINI_CELL = 15,
    NO_CONTACT     = 16,
};

namespace priv
{
template<typename T>
T
bnd_product(const T& fact, const T& func)
{
    return fact * func;
}

template<typename T, int N>
Matrix<T, N, 1>
bnd_product(const T& fact, const Matrix<T, N, 1>& func)
{
    return fact * func;
}

template<typename T, int N, int M>
Matrix<T, N, M>
bnd_product(const T& fact, const Matrix<T, N, M>& func)
{
    return fact * func;
}

template<typename scalar_type, size_t DIM, bool FunctionScalarType>
struct FunctionType
{
    typedef scalar_type function_type;
};

template<typename scalar_type, size_t DIM>
struct FunctionType<scalar_type, DIM, false>
{
    typedef static_vector<scalar_type, DIM> function_type;
};

template<bool ScalarBoundary>
struct num_face_dofs
{
    static size_t
    size(const size_t& face_degree, const size_t dimension)
    {
        return scalar_basis_size(face_degree, dimension - 1);
    }
};

template<>
struct num_face_dofs<false>
{
    static size_t
    size(const size_t& face_degree, const size_t dimension)
    {
        return vector_basis_size(face_degree, dimension - 1, dimension);
    }
};

template<bool ScalarBoundary>
struct imposed_dofs
{
    static size_t
    dirichlet_imposed_dofs(const size_t& btype, const size_t& num_face_dofs, const size_t dimension)
    {
        switch (btype)
        {
            case DIRICHLET: return num_face_dofs; break;
            case NOTHING: return 0; break;
            default: throw std::logic_error("Unknown Boundary condition"); break;
        }
    }
};

template<>
struct imposed_dofs<false>
{
    static size_t
    dirichlet_imposed_dofs(const size_t& btype, const size_t& num_face_dofs, const size_t dimension)
    {
        const size_t num_dim_dofs = num_face_dofs / dimension;
        switch (btype)
        {
            case DIRICHLET: return num_face_dofs; break;
            case CLAMPED: return num_face_dofs; break;
            case DX: return num_dim_dofs; break;
            case DY: return num_dim_dofs; break;
            case DZ:
                if (dimension != 3)
                    throw std::invalid_argument("You are not in 3D");
                return num_dim_dofs;
                break;
            case DXDY: return 2 * num_dim_dofs; break;
            case DXDZ:
                if (dimension != 3)
                    throw std::invalid_argument("You are not in 3D");
                return 2 * num_dim_dofs;
                break;
            case DYDZ:
                if (dimension != 3)
                    throw std::invalid_argument("You are not in 3D");
                return 2 * num_dim_dofs;
                break;
            case NOTHING: return 0; break;
            default: throw std::logic_error("Unknown Boundary condition"); break;
        }
    }
};
}

// class to create and impose boudary conditions
// ScalarBoundary = true for scalar problem like diffusion
// ScalarBoundary = false for vectorial problem like linear_elasticity

template<typename MeshType, bool ScalarBoundary>
class BoundaryConditions
{
  public:
    typedef MeshType                                                                                      mesh_type;
    typedef typename mesh_type::cell_type                                                                 cell_type;
    typedef typename mesh_type::face_type                                                                 face_type;
    typedef typename mesh_type::coordinate_type                                                           scalar_type;
    typedef point<scalar_type, mesh_type::dimension>                                                      point_type;
    typedef typename priv::FunctionType<scalar_type, mesh_type::dimension, ScalarBoundary>::function_type function_type;

  private:
    const mesh_type& m_msh;

    std::vector<std::function<function_type(point_type)>> m_dirichlet_func;
    std::vector<std::function<function_type(point_type)>> m_neumann_func;
    std::vector<std::function<function_type(point_type)>> m_robin_func;
    std::vector<std::function<scalar_type(point_type)>>   m_contact_func;

    // 1) bool to know if a boundary condition is associated to the face
    // 2) type of boundary conditions
    // 3) boundary id of the face
    // 4) boundary function id of the face
    typedef std::vector<std::tuple<bool, size_t, size_t, size_t>> bnd_storage_type;
    bnd_storage_type                                              m_faces_is_dirichlet;
    bnd_storage_type                                              m_faces_is_neumann;
    bnd_storage_type                                              m_faces_is_robin;
    std::vector<std::tuple<bool, size_t, size_t, int>>            m_faces_is_contact;

    size_t m_dirichlet_faces, m_neumann_faces, m_robin_faces, m_contact_faces;

    scalar_type m_factor;

    // search faces that have the boundary id "b_id"
    std::vector<size_t>
    search_faces(const size_t b_id) const
    {
        std::vector<size_t> list_faces;

        for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++)
        {
            const auto bfc = *itor;

            const auto face_id = m_msh.lookup(bfc);

            if (m_msh.boundary_id(face_id) == b_id)
            {
                list_faces.push_back(face_id);
            }
        }

        list_faces.shrink_to_fit();

        return list_faces;
    }

  public:
    BoundaryConditions() = delete;

    BoundaryConditions(const mesh_type& msh) :
      m_msh(msh), m_dirichlet_faces(0), m_neumann_faces(0), m_robin_faces(0), m_contact_faces(0), m_factor(1)
    {
        m_faces_is_dirichlet.assign(m_msh.faces_size(), std::make_tuple(false, NOTHING, 0, 0));
        m_faces_is_neumann.assign(m_msh.faces_size(), std::make_tuple(false, FREE, 0, 0));
        m_faces_is_robin.assign(m_msh.faces_size(), std::make_tuple(false, WHATEVER, 0, 0));
        m_faces_is_contact.assign(m_msh.faces_size(), std::make_tuple(false, NO_CONTACT, 0, -1));
    }

    void
    addContactBC(const size_t& btype, const size_t& b_id)
    {
        const auto list_faces = search_faces(b_id);

        for (size_t face_id : list_faces)
        {
            m_faces_is_contact.at(face_id) = std::make_tuple(true, btype, b_id, -1);
            m_contact_faces++;
        }
    }

    template<typename Function>
    void
    addContactBC(const size_t& btype, const size_t& b_id, const Function& bcf)
    {

        const size_t bcf_id = m_contact_func.size();
        m_contact_func.push_back(bcf);

        const auto list_faces = search_faces(b_id);

        for (size_t face_id : list_faces)
        {
            m_faces_is_contact.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_contact_faces++;
        }
    }

    template<typename Function>
    void
    addDirichletEverywhere(const Function& bcf)
    {
        const size_t bcf_id = m_dirichlet_func.size();
        m_dirichlet_func.push_back(bcf);

        for (auto itor = m_msh.boundary_faces_begin(); itor != m_msh.boundary_faces_end(); itor++)
        {
            const auto bfc = *itor;

            const auto face_id = m_msh.lookup(bfc);

            m_faces_is_dirichlet.at(face_id) = std::make_tuple(true, DIRICHLET, 0, bcf_id);
            m_dirichlet_faces++;
        }
    }

    template<typename Function>
    void
    addRobinBC(const size_t& btype, const size_t& b_id, const Function& bcf)
    {

        const size_t bcf_id = m_robin_func.size();
        m_robin_func.push_back(bcf);

        const auto list_faces = search_faces(b_id);

        for (size_t face_id : list_faces)
        {
            m_faces_is_robin.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_robin_faces++;
        }
    }

    template<typename Function>
    void
    addDirichletBC(const size_t& btype, const size_t& b_id, const Function& bcf)
    {
        const size_t bcf_id = m_dirichlet_func.size();
        m_dirichlet_func.push_back(bcf);

        const auto list_faces = search_faces(b_id);

        for (size_t face_id : list_faces)
        {
            m_faces_is_dirichlet.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_dirichlet_faces++;
        }
    }

    template<typename Function>
    void
    addNeumannBC(const size_t& btype, const size_t& b_id, const Function& bcf)
    {
        const size_t bcf_id = m_neumann_func.size();
        m_neumann_func.push_back(bcf);

        const auto list_faces = search_faces(b_id);

        for (size_t face_id : list_faces)
        {
            m_faces_is_neumann.at(face_id) = std::make_tuple(true, btype, b_id, bcf_id);
            m_neumann_faces++;
        }
    }

    template<typename Function>
    void
    updateDirichletFunction(const Function& bcf, size_t bcf_index)
    {
        assert(0 <= bcf_index && bcf_index < m_dirichlet_func.size());
        m_dirichlet_func[bcf_index] = bcf;
    }

    void
    multiplyAllFunctionsByAFactor(const scalar_type& factor)
    {
        m_factor = factor;
    }

    size_t
    nb_faces_boundary() const
    {
        return m_dirichlet_faces + m_neumann_faces + m_robin_faces + m_contact_faces;
    }

    size_t
    nb_faces_dirichlet() const
    {
        return m_dirichlet_faces;
    }

    size_t
    nb_faces_neumann() const
    {
        return m_neumann_faces;
    }

    size_t
    nb_faces_robin() const
    {
        return m_robin_faces;
    }
    size_t
    nb_faces_contact() const
    {
        return m_contact_faces;
    }

    bool
    is_dirichlet_face(const size_t face_i) const
    {
        return std::get<0>(m_faces_is_dirichlet.at(face_i));
    }

    bool
    is_dirichlet_face(const face_type& fc) const
    {
        return is_dirichlet_face(m_msh.lookup(fc));
    }

    bool
    is_neumann_face(const size_t face_i) const
    {
        return std::get<0>(m_faces_is_neumann.at(face_i));
    }

    bool
    is_neumann_face(const face_type& fc) const
    {
        return is_neumann_face(m_msh.lookup(fc));
    }

    bool
    is_robin_face(const size_t face_i) const
    {
        return std::get<0>(m_faces_is_robin.at(face_i));
    }

    bool
    is_robin_face(const face_type& fc) const
    {
        return is_robin_face(m_msh.lookup(fc));
    }

    bool
    is_contact_face(const size_t face_i) const
    {
        return std::get<0>(m_faces_is_contact.at(face_i));
    }

    bool
    is_contact_face(const face_type& fc) const
    {
        return is_contact_face(m_msh.lookup(fc));
    }

    size_t
    dirichlet_boundary_type(const size_t face_i) const
    {
        return std::get<1>(m_faces_is_dirichlet.at(face_i));
    }

    size_t
    dirichlet_boundary_type(const face_type& fc) const
    {
        return dirichlet_boundary_type(m_msh.lookup(fc));
    }

    size_t
    neumann_boundary_type(const size_t face_i) const
    {
        return std::get<1>(m_faces_is_neumann.at(face_i));
    }

    size_t
    neumann_boundary_type(const face_type& fc) const
    {
        return neumann_boundary_type(m_msh.lookup(fc));
    }

    size_t
    robin_boundary_type(const size_t face_i) const
    {
        return std::get<1>(m_faces_is_robin.at(face_i));
    }

    size_t
    robin_boundary_type(const face_type& fc) const
    {
        return robin_boundary_type(m_msh.lookup(fc));
    }

    size_t
    contact_boundary_type(const size_t face_i) const
    {
        return std::get<1>(m_faces_is_contact.at(face_i));
    }

    size_t
    contact_boundary_type(const face_type& fc) const
    {
        return contact_boundary_type(m_msh.lookup(fc));
    }

    size_t
    dirichlet_boundary_id(const size_t face_i) const
    {
        return std::get<2>(m_faces_is_dirichlet.at(face_i));
    }

    size_t
    dirichlet_boundary_id(const face_type& fc) const
    {
        return dirichlet_boundary_id(m_msh.lookup(fc));
    }

    size_t
    neumann_boundary_id(const size_t face_i) const
    {
        return std::get<2>(m_faces_is_neumann.at(face_i));
    }

    size_t
    neumann_boundary_id(const face_type& fc) const
    {
        return neumann_boundary_id(m_msh.lookup(fc));
    }

    size_t
    robin_boundary_id(const size_t face_i) const
    {
        return std::get<2>(m_faces_is_robin.at(face_i));
    }

    size_t
    robin_boundary_id(const face_type& fc) const
    {
        return robin_boundary_id(m_msh.lookup(fc));
    }

    size_t
    contact_boundary_id(const size_t face_i) const
    {
        return std::get<2>(m_faces_is_contact.at(face_i));
    }

    size_t
    contact_boundary_id(const face_type& fc) const
    {
        return robin_boundary_id(m_msh.lookup(fc));
    }

    bool
    cell_has_dirichlet_faces(const cell_type& cl) const
    {
        const auto fcs_id = faces_id(m_msh, cl);

        for (auto& fc_id : fcs_id)
        {
            if (is_dirichlet_face(fc_id))
            {
                return true;
            }
        }

        return false;
    }

    bool
    cell_has_robin_faces(const cell_type& cl) const
    {
        const auto fcs_id = faces_id(m_msh, cl);

        for (auto& fc_id : fcs_id)
        {
            if (is_robin_face(fc_id))
            {
                return true;
            }
        }

        return false;
    }

    bool
    cell_has_contact_faces(const cell_type& cl) const
    {
        const auto fcs_id = faces_id(m_msh, cl);

        for (auto& fc_id : fcs_id)
        {
            if (is_contact_face(fc_id))
            {
                return true;
            }
        }

        return false;
    }

    bool
    cell_has_neumann_faces(const cell_type& cl) const
    {
        const auto fcs_id = faces_id(m_msh, cl);

        for (auto& fc_id : fcs_id)
        {
            if (is_neumann_face(fc_id))
            {
                return true;
            }
        }

        return false;
    }

    size_t
    howmany_contact_faces(const cell_type& cl) const
    {
        const auto fcs_id = faces_id(m_msh, cl);
        size_t     nb  = 0;

        for (auto& fc_id : fcs_id)
        {
            if (is_contact_face(fc_id))
            {
                nb++;
            }
        }

        return nb;
    }

    size_t
    howmany_without_contact_faces(const cell_type& cl) const
    {
        return howmany_faces(m_msh, cl) - howmany_contact_faces(cl);
    }

    // In fact it returns the faces which have not a tag SIGNORINI
    std::vector<face_type>
    faces_without_contact(const cell_type& cl) const
    {
        const auto             fcs = faces(m_msh, cl);
        std::vector<face_type> ret;

        for (auto& fc : fcs)
        {
            if (contact_boundary_type(fc) == NO_CONTACT)
            {
                ret.push_back(fc);
            }
        }
        return ret;
    }

    // In fact it returns the faces which have not a tag SIGNORINI_CELL
    std::vector<face_type>
    faces_with_unknowns(const cell_type& cl) const
    {
        const auto             fcs = faces(m_msh, cl);
        std::vector<face_type> ret;

        for (auto& fc : fcs)
        {
            if (contact_boundary_type(fc) != SIGNORINI_CELL)
            {
                ret.push_back(fc);
            }
        }

        return ret;
    }

    // In fact it returns the faces which have a tag SIGNORINI
    std::vector<face_type>
    faces_with_contact(const cell_type& cl) const
    {
        const auto             fcs = faces(m_msh, cl);
        std::vector<face_type> ret;

        for (auto& fc : fcs)
        {
            if (contact_boundary_type(fc) != NO_CONTACT)
            {
                ret.push_back(fc);
            }
        }

        return ret;
    }

    auto
    dirichlet_boundary_func(const size_t face_i) const
    {
        if (!is_dirichlet_face(face_i))
        {
            throw std::logic_error("You want the Dirichlet function of face which is not a Dirichlet face");
        }

        const auto        func   = m_dirichlet_func.at(std::get<3>(m_faces_is_dirichlet.at(face_i)));
        const scalar_type factor = m_factor;

        auto rfunc = [ func, factor ](const point_type& p) -> auto { return priv::bnd_product(factor, func(p)); };

        return rfunc;
    }

    auto
    dirichlet_boundary_func(const face_type& fc) const
    {
        return dirichlet_boundary_func(m_msh.lookup(fc));
    }

    auto
    neumann_boundary_func(const size_t face_i) const
    {
        if (!is_neumann_face(face_i))
        {
            throw std::logic_error("You want the Neumann function of face which is not a Neumann face");
        }

        const auto        func   = m_neumann_func.at(std::get<3>(m_faces_is_neumann.at(face_i)));
        const scalar_type factor = m_factor;

        auto rfunc = [ func, factor ](const point_type& p) -> auto { return priv::bnd_product(factor, func(p)); };

        return rfunc;
    }

    auto
    neumann_boundary_func(const face_type& fc) const
    {
        return neumann_boundary_func(m_msh.lookup(fc));
    }

    auto
    robin_boundary_func(const size_t face_i) const
    {
        if (!is_robin_face(face_i))
        {
            throw std::logic_error("You want the Robin function of face which is not a Robin face");
        }

        const auto        func   = m_robin_func.at(std::get<3>(m_faces_is_robin.at(face_i)));
        const scalar_type factor = m_factor;

        auto rfunc = [ func, factor ](const point_type& p) -> auto
        {
            return disk::priv::inner_product(factor, func(p));
        };
        return rfunc;
    }

    auto
    robin_boundary_func(const face_type& fc) const
    {
        return robin_boundary_func(m_msh.lookup(fc));
    }

    auto
    contact_boundary_func(const size_t face_i) const
    {
        if (!is_contact_face(face_i))
        {
            throw std::logic_error("You want the contact function of a face which is not a conatact face");
        }

        const auto fid = std::get<3>(m_faces_is_contact.at(face_i));

        if (fid < 0)
        {
            throw std::logic_error("You want the contact function of a face which has not function");
        }

        return m_contact_func.at(fid);
    }

    auto
    contact_boundary_func(const face_type& fc) const
    {
        return contact_boundary_func(m_msh.lookup(fc));
    }

    void
    boundary_info() const
    {
        std::cout << "Number of boundary faces: " << nb_faces_boundary() << std::endl;
        std::cout << "including: " << std::endl;
        std::cout << " - Number of Dirichlet faces: " << nb_faces_dirichlet() << std::endl;
        std::cout << " - Number of Neumann faces: " << nb_faces_neumann() << std::endl;
        std::cout << " - Number of Robin faces: " << nb_faces_robin() << std::endl;
        std::cout << " - Number of Contact faces: " << nb_faces_contact() << std::endl;
    }

    size_t
    dirichlet_imposed_dofs(const size_t& face_id, const size_t face_degree) const
    {
        const size_t dimension     = mesh_type::dimension;
        const size_t num_face_dofs = priv::num_face_dofs<ScalarBoundary>::size(face_degree, dimension);

        if (is_dirichlet_face(face_id))
        {
            const size_t btype = dirichlet_boundary_type(face_id);

            return priv::imposed_dofs<ScalarBoundary>::dirichlet_imposed_dofs(btype, num_face_dofs, dimension);
        }

        return 0;
    }
};

template<typename MeshType>
using scalar_boundary_conditions = BoundaryConditions<MeshType, true>;

template<typename MeshType>
using vector_boundary_conditions = BoundaryConditions<MeshType, false>;
} // end disk
