/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2025
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"
#include "diskpp/methods/implementation_hho/curl.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho_assemblers.hpp"
#include "mumps.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/output/silo.hpp"

template<typename Mesh>
using DM = disk::dynamic_matrix<typename Mesh::coordinate_type>;

template<typename Mesh>
using DV = disk::dynamic_vector<typename Mesh::coordinate_type>;

template<typename Mesh>
std::pair<DM<Mesh>, DM<Mesh>>
hho_mixedhigh_symlapl(const Mesh& msh,
    const typename Mesh::cell_type& cl, size_t degree)
{
    const static size_t DIM = Mesh::dimension;
    static_assert(DIM==2 or DIM==3, "Symmetric laplacian: only DIM = 2 or DIM=3");

    using scalar_type = typename Mesh::coordinate_type;
    /* Reconstruction space basis */
    auto rb = disk::make_vector_monomial_basis(msh, cl, degree+1);
    auto rbs = rb.size();

    auto fcs = faces(msh, cl);
    auto fbs = disk::vector_basis_size(degree, DIM-1, DIM);
    auto n_allfacedofs = fcs.size() * fbs;

    /* Stiffness */
    disk::dynamic_matrix<scalar_type> K =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);

    /* Nitsche contributions (consistency, symmetry, penalization) */
    disk::dynamic_matrix<scalar_type> LHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs);
    
    /* Local problem RHS */
    disk::dynamic_matrix<scalar_type> RHS =
        disk::dynamic_matrix<scalar_type>::Zero(rbs, rbs + n_allfacedofs);
}

int main(void) {

    using T = double;
    using mesh_type = disk::cartesian_mesh<T,2>;

    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    mesher.refine();

    for (auto& cl : msh) {
        hho_mixedhigh_symlapl(msh, cl, 1);
    }
}