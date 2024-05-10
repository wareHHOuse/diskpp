#include <vector>
#include <cassert>
#include <utility>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/geometry/geometry.hpp"

#define TEST_MESH "../../../meshes/2D_hex/fvca5/hexagonal_1.typ1"

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

        /* faces() gives you the faces in lex order, but we need CCW.
         * This sucks (it's n^2) but does the job. */
        std::vector<std::pair<size_t,bool>> reorder;
        for (size_t i = 0; i < ptids.size(); i++)
        {
            bool flipped = false;
            auto n0 = ptids[i];
            auto n1 = ptids[(i+1)%ptids.size()];
            if (n0 > n1) {
                std::swap(n0, n1);
                flipped = true;
            }
            /* This builds the lex-to-CCW table, it also stores if the
             * edge is flipped. */
            for (size_t j = 0; j < fcs.size(); j++) {
                auto fc_ptids = fcs[j].point_ids();
                assert(fc_ptids.size() == 2);
                if (fc_ptids[0] == n0 and fc_ptids[1] == n1)
                    reorder.push_back({j,flipped});
            }
        }

        assert(reorder.size() == fcs.size());
        std::cout << "lex-to-CCW: ";
        for (auto& r : reorder)
            std::cout << "(" << r.first << " " << std::boolalpha << r.second << ") ";
        std::cout << std::endl;

        /* OK, now we build the whole element l2g */
        for (size_t i = 0; i < fcs.size(); i++) {
            auto [ofs, flip] = reorder[i];
            auto fc = fcs[ofs];
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
    disk::generic_mesh<T, 2> msh;
    disk::load_mesh_fvca5_2d<T>(TEST_MESH, msh);
    run_vemasm(msh, 3);
}