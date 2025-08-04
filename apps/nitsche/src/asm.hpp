#pragma once

/***************************************************************
 * Boundary conditions helpers
 */

enum bc {
    none,
    dirichlet,
    neumann,
};

template<typename Mesh>
void
set_boundary(const Mesh& msh, std::vector<bc>& bcs, bc bc_type, size_t bnd)
{
    if (bcs.size() != msh.faces_size()) {
        bcs.resize( msh.faces_size() );
    }

    size_t fcnum = 0;
    for (auto& fc : faces(msh)) {
        auto bi = msh.boundary_info(fc);
        if (bi.is_boundary() and bi.id() == bnd) {
            bcs[fcnum] = bc_type;
        }
        fcnum++;
    }   
}