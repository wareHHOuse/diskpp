

/* If called with `base = 0`, simply make the subdomain
 * numbering zero-based. Otherwise, make the subdomain
 * numbering start from `base`. */
template<typename Mesh>
void rebase_subdomain_numbering(Mesh& msh, size_t base = 0)
{
    if (msh.cells_size() == 0)
        return;

    auto di0 = msh.domain_info(msh.cell_at(0));
    auto min_tag = di0.tag();

    for (size_t i = 1; i < msh.cells_size(); i++) {
        const auto& cl = msh.cell_at(i);
        auto di = msh.domain_info(cl);
        min_tag = std::min(min_tag, di.tag());
    }

    auto storage = msh.backend_storage();
    for (auto& di : storage->subdomain_info)
        di = disk::subdomain_descriptor( (di.tag() - min_tag) + base);
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
    disk::renumber_hypercube_boundaries(msh);
}

template<disk::mesh_2D FineMesh>
using coarse_mesh_t = disk::generic_mesh<typename FineMesh::coordinate_type, 2>;

template<typename FineMesh>
void agglomerate_by_subdomain(const FineMesh& srcmsh,
    coarse_mesh_t<FineMesh>& dstmsh)
{
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

    using src_mesh_type = FineMesh;
    using coord_type = typename FineMesh::coordinate_type;
    using dst_mesh_type = coarse_mesh_t<FineMesh>;
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
        
        /* Reverse the path if vertices are in clockwise order 
         * (uses shoelace formula).*/
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

    /* We have all the edges: sort them and make them unique */
    std::sort(dststor->edges.begin(), dststor->edges.end());
    auto last = std::unique(dststor->edges.begin(), dststor->edges.end());
    dststor->edges.erase(last, dststor->edges.end());

    using newsurf_pair_t = std::pair<size_t, dst_surf_type>;
    using newsurfs_vec_t = std::vector<newsurf_pair_t>;
    newsurfs_vec_t newsurfs;
    newsurfs.reserve(paths.size());

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

        /* We store in [tag, surf] pairs to not lose the association
         * during sorting. Later we unzip the pair and store the
         * data in the mesh. */
        dst_surf_type newsurf(surfedges);
        newsurf.set_point_ids(path.begin(), path.end());
        newsurfs.push_back({tag, newsurf});
    }

    /* Sort all the tag-surface pairs according to the surface ordering */
    std::sort(newsurfs.begin(), newsurfs.end(),
        [](const newsurf_pair_t& a, const newsurf_pair_t& b) {
            return a.second < b.second;
        }
    );

    /* Unzip the newsurf array in the appropriate storage arrays */
    dststor->surfaces.clear();
    dststor->subdomain_info.clear();
    dststor->surfaces.reserve( newsurfs.size() );
    dststor->subdomain_info.reserve( newsurfs.size() );
    for (const auto& [tag, newsurf] : newsurfs) {
        dststor->surfaces.push_back(newsurf);
        auto di = disk::subdomain_descriptor(tag);
        dststor->subdomain_info.push_back(di);
    }
    assert( dststor->surfaces.size() == newsurfs.size() );
    assert( dststor->subdomain_info.size() == newsurfs.size() );

    dststor->boundary_info.resize( dststor->edges.size() );

    /* Count the elements adjacent to each face */
    std::vector<size_t> counts(dststor->edges.size());
    for (auto& cl : dstmsh) {
        auto fcs = faces(dstmsh, cl);
        for (auto& fc : fcs) {
            counts[offset(dstmsh, fc)]++;
        }
    }

    /* All the faces adjacent to only one element are
     * boundary faces. This should preserve the original
     * boundaries, but we leave this for later. Not
     * necessary for now. */
    for ( size_t i = 0; i < counts.size(); i++) {
        if (counts[i] == 1) {
            dststor->boundary_info[i].is_boundary(true);
            dststor->boundary_info[i].id(0);
        }
    }

    disk::mark_internal_faces(dstmsh);
}

template<typename FM>
using cc2ff_t = std::map<typename coarse_mesh_t<FM>::cell_type, std::set<typename FM::face_type>>;

template<typename FineMesh>
auto make_cc2ff(const FineMesh& fmsh,
    const coarse_mesh_t<FineMesh>& cmsh)
{
    using fine_mesh_type = FineMesh;
    using coord_type = typename fine_mesh_type::coordinate_type;
    using coarse_mesh_type = coarse_mesh_t<FineMesh>;
    using fine_face_type = typename fine_mesh_type::face_type;
    using coarse_cell_type = typename coarse_mesh_type::cell_type;

    /* This maps from Coarse Cells to Fine Faces */
    cc2ff_t<fine_mesh_type> cc2ff;

    for (const auto& fcl : fmsh) {
        auto di = fmsh.domain_info(fcl);
        auto ccl = cmsh.cell_at(di.tag());
        auto ffcs = faces(fmsh, fcl);
        cc2ff[ccl].insert(ffcs.begin(), ffcs.end());
    }

    std::vector<int> occs(fmsh.faces_size());
    for (const auto& [ccl, ffcs] : cc2ff) {
        for (const auto& ffc : ffcs) {
            occs[ offset(fmsh, ffc) ]++;
        }
    }

    int min = occs[0];
    int max = occs[0];
    for (size_t i = 1; i < occs.size(); i++) {
        min = std::min(min, occs[i]);
        max = std::max(max, occs[i]);
    }

    return std::pair(cc2ff, occs);
}