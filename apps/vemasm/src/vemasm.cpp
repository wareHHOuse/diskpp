#include <vector>
#include <cassert>
#include <utility>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/methods/vem"

#define TEST_MESH_TRI  "../../../meshes/2D_triangles/netgen/tri01.mesh2d"
#define TEST_MESH_QUAD "../../../meshes/2D_quads/diskpp/testmesh-2-2.quad"
#define TEST_MESH_POLY "../../../meshes/2D_hex/fvca5/hexagonal_1.typ1"

template<disk::mesh_2D Mesh>
void run_vemasm(const Mesh& msh, size_t k)
{
    disk::dof_mapper dofmap(msh, k);

    for (size_t cl_i = 0; cl_i < msh.cells_size(); cl_i++)
    {
        const auto& cl = msh[cl_i];
        std::cout << "*************" << std::endl;
        auto ptids = cl.point_ids();
        std::cout << "element point ids: ";
        for (auto& ptid : ptids)
            std::cout << ptid << " ";
        std::cout << std::endl;

        /* Vector for the whole cell l2g */
        std::vector<size_t> cl_l2g = dofmap.cell_to_global(cl_i);

        std::cout << "element l2g: ";
        for (auto& p : cl_l2g)
            std::cout << p << " ";
        std::cout << std::endl;
    }
}

int main(void)
{
    using T = double;

    std::cout << " *** POLYHEDRAL ***" << std::endl;
    disk::generic_mesh<T, 2> msh_poly;
    disk::load_mesh_fvca5_2d<T>(TEST_MESH_POLY, msh_poly);
    run_vemasm(msh_poly, 3);

    std::cout << " *** CARTESIAN ***" << std::endl;
    disk::cartesian_mesh<T, 2> msh_quad;
    disk::load_mesh_diskpp_cartesian(TEST_MESH_QUAD, msh_quad);
    run_vemasm(msh_quad, 3);

    std::cout << " *** SIMPLICIAL ***" << std::endl;
    disk::simplicial_mesh<T, 2> msh_tri;
    disk::load_mesh_netgen(TEST_MESH_TRI, msh_tri);
    run_vemasm(msh_tri, 3);

    return 0;
}