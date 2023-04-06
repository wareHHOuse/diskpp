/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/* This tool makes the FVCA5 hexagonal mesh sequence, but much faster than
 * the original matlab code. */

#include <sstream>

#include "diskpp/mesh/meshgen_fvca5hex.hpp"


bool
write_typ1(const disk::fvca5::hexmesh& hm, const std::string& fn)
{
    std::cout << "Writing" << std::endl;
    std::ofstream ofs(fn);
    if ( !ofs.is_open() )
        return false;
    
    ofs << "vertices" << std::endl;
    ofs << hm.points.size() << std::endl;
    for (const auto& point : hm.points) {
        ofs << double(point.x)/hm.coord_max << '\t';
        ofs << double(point.y)/hm.coord_max << std::endl;
    }
    
    size_t num_triangles = 0;
    size_t num_quadrangles = 0;
    size_t num_pentagons = 0;
    size_t num_hexagons = 0;
    for (const auto& polygon : hm.polygons) {
        switch (polygon.points.size()) {
            case 3: num_triangles++; break;
            case 4: num_quadrangles++; break;
            case 5: num_pentagons++; break;
            case 6: num_hexagons++; break;
            default:
                throw std::logic_error("wrong number of vertices, this is a bug.");
        }
    }
    
    auto itor = hm.polygons.begin();
    ofs << "triangles" << std::endl;
    ofs << num_triangles << std::endl;
    for(; itor != hm.polygons.end() and (*itor).points.size() < 4; itor++)
    {
        const auto& points = (*itor).points;
        assert(points.size() == 3);
        for (size_t i = 0; i < points.size(); i++)
            ofs << disk::fvca5::offset(hm, points[i], true) << '\t';
        ofs << std::endl;
    }

    ofs << "quadrangles" << std::endl;
    ofs << num_quadrangles << std::endl;
    for(; itor != hm.polygons.end() and (*itor).points.size() < 5; itor++)
    {
        const auto& points = (*itor).points;
        assert(points.size() == 4);
        for (size_t i = 0; i < points.size(); i++)
            ofs << disk::fvca5::offset(hm, points[i], true) << '\t';
        ofs << std::endl;
    }
    
    ofs << "pentagons" << std::endl;
    ofs << num_pentagons << std::endl;
    for(; itor != hm.polygons.end() and (*itor).points.size() < 6; itor++)
    {
        const auto& points = (*itor).points;
        assert(points.size() == 5);
        for (size_t i = 0; i < points.size(); i++)
            ofs << disk::fvca5::offset(hm, points[i], true) << '\t';
        ofs << std::endl;
    }
    
    ofs << "hexagons" << std::endl;
    ofs << num_hexagons << std::endl;
    for(; itor != hm.polygons.end() and (*itor).points.size() < 7; itor++)
    {
        const auto& points = (*itor).points;
        assert(points.size() == 6);
        for (size_t i = 0; i < points.size(); i++)
            ofs << disk::fvca5::offset(hm, points[i], true) << '\t';
        ofs << std::endl;
    }
    
    size_t num_bndedges = 0;
    for (const auto& [boundary, edges] : hm.boundary_edges)
        num_bndedges += edges.size();
    
    ofs << "edges of the boundary" << std::endl;
    ofs << num_bndedges << std::endl;
    for (const auto& [boundary, edges] : hm.boundary_edges) {
        for (const auto& edge : edges) {
            ofs << disk::fvca5::offset(hm, edge.a, true) << '\t';
            ofs << disk::fvca5::offset(hm, edge.b, true) << std::endl;
        }
    }
    
    auto owner1 = [&](const disk::fvca5::edge& e) {
        const auto eo = hm.edge_owners.at(e);
        return eo.first.value()+1;
    };
    
    auto owner2 = [&](const disk::fvca5::edge& e) {
        auto eo = hm.edge_owners.at(e);
        if (eo.second)
            return eo.second.value()+1;
        return size_t(0);
    };
    
    ofs << "all edges" << std::endl;
    ofs << hm.edges.size() << std::endl;
    for (const auto& edge : hm.edges) {
        ofs << disk::fvca5::offset(hm, edge.a, true) << '\t';
        ofs << disk::fvca5::offset(hm, edge.b, true) << '\t' << owner1(edge);
        ofs << '\t' << owner2(edge) << std::endl;
    }
    
    return true;
}

int main(int argc, char **argv)
{
    size_t level = 1;
    std::string filename = "xxx_hexagonal_1.typ1";

    if (argc == 2) {
        level = std::stoi(argv[1]);
        std::stringstream ss;
        ss << "xxx_hexagonal_" << level << ".typ1";
        filename = ss.str();
    }

    if (argc > 2) {
        level = std::stoi(argv[1]);
        filename = argv[2];
    }

    disk::fvca5::hexmesh hm;
    disk::fvca5::make_hexmesh(hm, level);
    write_typ1(hm, filename);

    return 0;
}

