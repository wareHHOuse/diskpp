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
#include <map>
#include <set>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/output/silo.hpp"

struct adjlist {
    size_t nodes[2];
    size_t used;

    adjlist() : used(0) {}

    void insert(size_t n) {
        if (used > 1)
            throw std::logic_error("Adjacency list full.");

        nodes[used++] = n;
    }
};

template<typename SrcMesh>
void agglomerate_by_subdomain(const SrcMesh& srcmsh,
    disk::generic_mesh<typename SrcMesh::coordinate_type, 2>& dstmsh)
{
    using src_mesh_type = SrcMesh;
    using coord_type = typename SrcMesh::coordinate_type;
    using dst_mesh_type = disk::generic_mesh<coord_type, 2>;
    using src_face_type = typename src_mesh_type::face_type;
    using dst_node_type = typename dst_mesh_type::node_type;
    using dst_edge_type = typename dst_mesh_type::edge_type;
    using dst_surf_type = typename dst_mesh_type::surface_type;

    using polygraph_t = std::map<size_t, adjlist>;
    std::map<size_t, polygraph_t> pgs;
    std::vector<std::optional<size_t>> compress_map;
    compress_map.resize(srcmsh.points_size());

    /* Collect all the boundary edges of the subdomains we want to
     * agglomerate in a polygon. Make a graph out of them. */
    for (auto& cl : srcmsh) {
        auto di = srcmsh.domain_info(cl);
        auto& pg = pgs[di.tag()];
        auto fcs = faces(srcmsh, cl);
        for (auto& fc : fcs) {
            auto bi = srcmsh.boundary_info(fc);
            if (not bi.is_boundary())
                continue;
            auto [pi1, pi2] = fc.point_ids();
            pg[pi1].insert(pi2);
            pg[pi2].insert(pi1);
            compress_map[pi1] = 1;
            compress_map[pi2] = 1;
        }
    }
    
    /* Make the original mesh to new mesh node mapping. Collect the
     * necessary points and node numbers. */
    auto srcstor = srcmsh.backend_storage();
    auto dststor = dstmsh.backend_storage();
    for (size_t i = 0, ci = 0; i < compress_map.size(); i++) {
        if (compress_map[i]) {
            dststor->points.push_back(srcstor->points[i]);
            dststor->nodes.push_back(typename dst_mesh_type::node_type{ci});
            compress_map[i] = ci++;
        }
    }

    /* Do a DFS to build the polygon out of its edges. This is a special
     * case of DFS as each node is guaranteed to have only two neighbours.*/
    std::map<size_t, std::vector<size_t>> paths;
    for (auto& [tag, pg] : pgs) {
        auto nvtx = pg.size();
        assert(nvtx >= 3);
        std::vector<size_t> path;
        path.reserve(nvtx);
        auto [start, adj] = *pg.begin();
        assert(adj.used == 2);
        size_t visiting = adj.nodes[0];
        path.push_back(start);
        for (size_t i = 1; i < nvtx; i++) {
            path.push_back(visiting);
            adj = pg.at(visiting);
            visiting = (adj.nodes[0] == path[i-1]) ? adj.nodes[1] : adj.nodes[0];
        }
        assert(visiting == start);
        
        /* Reverse the path if vertices are in clockwise order */
        double dir = 0.0;
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = srcstor->points[ path[i] ];
            auto p1 = srcstor->points[ path[(i+1)%path.size()] ];
            dir += (p1.x() - p0.x())*(p1.y() + p0.y());
        }

        if (dir > 0) std::reverse(path.begin(), path.end());

        /* and translate the node numbers to the new mesh numbering */
        for (auto & n : path)
            n = compress_map[ n ].value();

        /* Now put all the edges in the mesh */
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = path[i];
            auto p1 = path[(i+1)%path.size()];
            auto node1 = typename dst_node_type::id_type(p0);
            auto node2 = typename dst_node_type::id_type(p1);
            auto e = dst_edge_type(node1, node2);
            dststor->edges.push_back(e);
        }

        paths[tag] = std::move(path);
    }

    /* Sort the edges and make them unique */
    std::sort(dststor->edges.begin(), dststor->edges.end());
    auto last = std::unique(dststor->edges.begin(), dststor->edges.end());
    dststor->edges.erase(last, dststor->edges.end());

    /* and finally all the elements */
    for (auto& [tag, path] : paths) {
        std::vector<typename dst_edge_type::id_type> surfedges;
        for (size_t i = 0; i < path.size(); i++) {
            auto p0 = path[i];
            auto p1 = path[(i+1)%path.size()];
            auto node1 = typename dst_node_type::id_type(p0);
            auto node2 = typename dst_node_type::id_type(p1);
            auto e = dst_edge_type(node1, node2);
            auto eid = find_element_id(dststor->edges.begin(),
                dststor->edges.end(), e);
            surfedges.push_back(eid.second);
        }

        dst_surf_type newsurf(surfedges);
        newsurf.set_point_ids(path.begin(), path.end());
        dststor->surfaces.push_back(newsurf);
    }

    std::sort(dststor->surfaces.begin(), dststor->surfaces.end());
    //disk::mark_internal_faces(dstmsh);
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

template<disk::mesh_2D Mesh>
void dump_mesh(const Mesh& msh)
{
    for (size_t cli = 0; cli < msh.cells_size(); cli++) {
        auto cl = msh.cell_at(cli);

        std::stringstream ss;
        ss << "cell_" << cli << ".m";
        std::ofstream ofs(ss.str());
        auto bar = barycenter(msh, cl);
        ofs << "plot(" << bar.x() << "," << bar.y() << ",'*');\n";
        ofs << "hold on;";
        auto pts = points(msh, cl);
        for (size_t i = 0; i < pts.size(); i++) {
            auto p0 = pts[i];
            auto p1 = pts[(i+1)%pts.size()];
            if (i == 0)
            {
                ofs << "line([" << p0.x() << ", " << p1.x() << "], "
                            "[" << p0.y() << ", " << p1.y() << "], 'Color', 'red');\n";
            }
            else
            {
                ofs << "line([" << p0.x() << ", " << p1.x() << "], "
                            "[" << p0.y() << ", " << p1.y() << "], 'Color', 'blue');\n";
            }
        }
    }
}

int main(void)
{
    using T = double;
    
    
    disk::simplicial_mesh<T,2> srcmsh;
    auto mesher = make_simple_mesher(srcmsh);
    for (size_t l = 0; l < 1; l++)
        mesher.refine();
    
/*
    disk::generic_mesh<T,2> srcmsh;
    auto mesher = make_fvca5_hex_mesher(srcmsh);
    mesher.make_level(2);
*/
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
    silo.add_mesh(dstmsh, "dstmesh");
    silo.add_variable("srcmesh", "partnum", cp, disk::zonal_variable_t);
    
    std::map<size_t, double> areas;
    for (auto& scl : srcmsh) {
        auto tag = srcmsh.domain_info(scl).tag();
        areas[tag] += measure(srcmsh, scl);
    }

    for (auto& [tag, area] : areas)
        std::cout << tag << " " << area << std::endl;

    for (auto& dcl : dstmsh) {
        std::cout << dcl << std::endl;
        auto area = measure(dstmsh, dcl);
        auto fcs = faces(dstmsh, dcl);
        std::cout << area << " " << fcs.size() << std::endl;
        for (auto& fc : fcs)
            std::cout << normal(dstmsh, dcl, fc).transpose() << std::endl;
    }

    dump_mesh(dstmsh);

    return 0;
}