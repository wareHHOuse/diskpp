/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <vector>
#include <set>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

struct edge {
    size_t p1;
    size_t p2;
};

std::ostream& operator<<(std::ostream& os, const edge& e) {
    os << "[" << e.p1 << ", " << e.p2 << "]";
    return os;
}

template<typename SrcMesh>
void agglomerate_by_subdomain(const SrcMesh& srcmsh,
    disk::generic_mesh<typename SrcMesh::coordinate_type, 2>& dstmsh)
{
    using src_mesh_type = SrcMesh;
    using coord_type = typename SrcMesh::coordinate_type;
    using dst_mesh_type = disk::generic_mesh<coord_type, 2>;
    using src_face_type = typename src_mesh_type::face_type;

    std::map<size_t, std::vector<edge>> elem_faces;

    for (auto& cl : srcmsh) {
        auto di = srcmsh.domain_info(cl);
        auto fcs = faces(srcmsh, cl);
        for (auto& fc : fcs) {
            auto bi = srcmsh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;

            auto [pi1, pi2] = fc.point_ids();

            elem_faces[di.tag()].push_back({pi1, pi2});
        }
    }

    for (auto& [tag, edges] : elem_faces) {
        assert(edges.size() > 2);
        for (size_t i = 0; i < edges.size(); i++) {
            for (size_t j = 0; j < edges.size(); j++) {
            if (edges[i-1].p2 == edges[i].p2)
                std::swap(edges[i].p1, edges[i].p2);
            }
        }
    }

    for (auto& [tag, edges] : elem_faces) {
        for (auto& e : edges) 
            std::cout << e << " ";
        std::cout << std::endl;
    }
}

template<typename Mesh>
void
partition_unit_square_mesh(Mesh& msh, size_t np)
{
    if (np < 2)
        return;

    auto storage = msh.backend_storage();
    for (size_t cell_i = 0; cell_i < msh.cells_size(); cell_i++) {
        auto cl = msh.cell_at(cell_i);
        auto bar = barycenter(msh, cl);
        auto domxy = bar * np;

        size_t subdom = 0;
        if constexpr (Mesh::dimension == 1)
            subdom = size_t(domxy.x());
        if constexpr (Mesh::dimension == 2)
            subdom = size_t(domxy.x()) + np*size_t(domxy.y());
        if constexpr (Mesh::dimension == 3)
            subdom = size_t(domxy.x()) + np*(size_t(domxy.y()) + np*size_t(domxy.z()));

        storage->subdomain_info[cell_i] = disk::subdomain_descriptor(subdom);
    }

    disk::make_interpartition_boundaries(msh);
}

int main(void)
{
    using T = double;
    disk::simplicial_mesh<T,2> srcmsh;
    auto mesher = make_simple_mesher(srcmsh);
    for (size_t l = 0; l < 1; l++)
        mesher.refine();
        
    partition_unit_square_mesh(srcmsh, 2);

    std::vector<double> cp;
    for (auto& cl : srcmsh) {
        auto di = srcmsh.domain_info(cl);
        cp.push_back(di.tag());
    }


    disk::generic_mesh<T, 2> dstmsh;

    agglomerate_by_subdomain(srcmsh, dstmsh);


    disk::silo_database silo;
    silo.create("agglo.silo");
    silo.add_mesh(srcmsh, "srcmesh");
    silo.add_variable("srcmesh", "partnum", cp, disk::zonal_variable_t);
    return 0;
}