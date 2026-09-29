/*
 * spatial_operator.hpp : the u-independent, boundary-independent bulk.
 *
 * Owns the discretization (MeshDegreeInfo, LS: cell at k+1, faces at k) and
 * caches, per cell:
 *
 *     K[ci] = 2 mu gr.second + lam dr.second + 2 mu stab      (elasticity)
 *     G[ci] = gr.first                                        (gradrec)
 *     M[ci] = cell-block mass matrix                          (optional)
 */

#pragma once

#include <vector>

#include "diskpp/methods/hho"
#include "diskpp/bases/bases.hpp"
#include "diskpp/adaptivity/adaptivity.hpp"

namespace hho_contact {

template<typename Mesh>
class spatial_operator
{
public:
    using T           = typename Mesh::coordinate_type;
    using matrix_type = disk::dynamic_matrix<T>;
    using vector_type = disk::dynamic_vector<T>;
    using cell_type   = typename Mesh::cell_type;

    const Mesh&               msh;
    disk::MeshDegreeInfo<Mesh> di;

    size_t k;            // face degree
    size_t cell_deg;     // k+1 (Lehrenfeld-Schoberl)
    size_t cbs;          // cell basis size
    T      mu, lam;

    std::vector<matrix_type> K;   // per-cell elasticity, full local dofs
    std::vector<matrix_type> G;   // per-cell gradient reconstruction
    std::vector<matrix_type> M;   // per-cell cell-block mass (empty if not requested)

    spatial_operator(const Mesh& m, size_t k_, T mu_, T lam_,
                     bool with_mass = false)
        : msh(m), di(m, k_ + 1, k_), k(k_), cell_deg(k_ + 1),
          mu(mu_), lam(lam_)
    {
        constexpr size_t DIM = Mesh::dimension;
        cbs = disk::vector_basis_size(cell_deg, DIM, DIM);

        const size_t nc = msh.cells_size();
        K.resize(nc);
        G.resize(nc);
        if (with_mass) M.resize(nc);

        size_t ci = 0;
        for (auto& cl : msh)
        {
            auto gr   = disk::make_matrix_hho_symmetric_gradrec(msh, cl, di);
            auto dr   = disk::make_hho_divergence_reconstruction(msh, cl, di);
            auto stab = disk::make_vector_hdg_stabilization(msh, cl, di);

            K[ci] = 2.0*mu*gr.second + lam*dr.second + 2.0*mu*stab;
            G[ci] = gr.first;

            if (with_mass)
            {
                auto cb = disk::make_vector_monomial_basis(msh, cl, cell_deg);
                M[ci]  = disk::make_mass_matrix(msh, cl, cb);
            }
            ci++;
        }
    }

    size_t num_cells() const { return msh.cells_size(); }

    size_t faces_dofs(const cell_type& cl) const
    {
        return disk::vector_faces_dofs(msh, cl, di);
    }

    size_t total_dofs(const cell_type& cl) const
    {
        return cbs + faces_dofs(cl);
    }
};

} // namespace hho_contact
