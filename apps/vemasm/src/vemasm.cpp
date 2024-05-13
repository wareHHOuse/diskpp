#include <vector>
#include <cassert>
#include <utility>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/methods/vem"

#define TEST_MESH_POLY "../../../meshes/2D_hex/fvca5/hexagonal_1.typ1"

#define TEST_MESH_QUAD "../../../meshes/2D_quads/diskpp/testmesh-2-2.quad"

template<disk::mesh_2D Mesh>
void run_vemasm(const Mesh& msh, size_t k)
{
    if (k < 1) {
        std::cout << "k must be greater than 0" << std::endl;
        return;
    }

    /* Offset of the first edge-based dof */
    size_t pts_base = msh.points_size();

    /* For each edge in the mesh (edges are ordered lexicographically)
     * compute the offset of all its edge-based dofs */
    std::vector<std::vector<size_t>> all_l2g;
    for (auto& fc : faces(msh)) {
        auto ptids = fc.point_ids();
        std::vector<size_t> l2g(k+1);
        l2g[0] = ptids[0];
        for (size_t i = 1; i < k; i++)
            l2g[i] = pts_base++;
        l2g[k] = ptids[1];
        all_l2g.push_back( std::move(l2g) );
    }

    for (auto& cl : msh)
    {
        std::cout << "*************" << std::endl;
        auto ptids = cl.point_ids();
        std::cout << "element point ids: ";
        for (auto& ptid : ptids)
            std::cout << ptid << " ";
        std::cout << std::endl;

        /* Vector for the whole cell l2g */
        std::vector<size_t> cl_l2g;

        auto fcs = faces(msh, cl);

        for (auto& fc : fcs) {
            std::cout << fc << std::endl;
        }

        auto fcs_ccw = faces_ccw(msh, cl);
        for (auto& [fc, flip] : fcs_ccw) {
            std::cout << fc << " " << std::boolalpha << flip << std::endl;
        }

        /* OK, now we build the whole element l2g */
        for (size_t i = 0; i < fcs_ccw.size(); i++) {
            auto [fc, flip] = fcs_ccw[i];
            std::cout << "  " << fc << std::endl;
            std::cout << "    ";
            auto& l2g = all_l2g[ offset(msh, fc) ];
            for (auto& ofs : l2g)
                std::cout << ofs << " ";
            std::cout << std::endl;
            /* The last node of this segment is the first of
             * the next, so we discard it. */
            if (flip)
                cl_l2g.insert(cl_l2g.end(), l2g.rbegin(), l2g.rend()-1);
            else
                cl_l2g.insert(cl_l2g.end(), l2g.begin(), l2g.end()-1);
        }
        std::cout << "element l2g: ";
        for (auto& p : cl_l2g)
            std::cout << p << " ";
        std::cout << std::endl;
    }
}

int main(void)
{
    using T = double;
    disk::generic_mesh<T, 2> msh_poly;
    disk::load_mesh_fvca5_2d<T>(TEST_MESH_POLY, msh_poly);
    run_vemasm(msh_poly, 3);

    std::cout << " *** CARTESIAN ***" << std::endl;
    disk::cartesian_mesh<T, 2> msh_quad;
    disk::load_mesh_diskpp_cartesian(TEST_MESH_QUAD, msh_quad);
    run_vemasm(msh_quad, 3);

    return 0;
}