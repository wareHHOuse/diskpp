/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include <iostream>
#include <optional>
#include <fstream>
#include <vector>
#include <map>

#include <cassert>
#include <cmath>
#include <algorithm>

namespace disk::fvca5 {

enum class boundary {
    BOTTOM = 0,
    RIGHT = 1,
    TOP = 2,
    LEFT = 3
};

struct point {
    size_t x;
    size_t y;

    point(size_t px, size_t py)
        : x(px), y(py)
    {}

    void make_absolute(size_t base_x, size_t base_y)
    {
        x += base_x;
        y += base_y;
    }

    bool operator<(const point& other) const
    {
        return (x < other.x) or (x == other.x and y < other.y);
    }

    bool operator==(const point& other) const
    {
        return x == other.x and y == other.y;
    }
};

struct edge {
    point   a;
    point   b;

    edge(const point& pa, const point& pb)
        : a(pa), b(pb)
    {
        assert( not(a == b) );
        if (b < a)
            std::swap(a, b);
    }

    void make_absolute(size_t base_x, size_t base_y)
    {
        a.make_absolute(base_x, base_y);
        b.make_absolute(base_x, base_y);
    }

    bool operator<(const edge& other) const
    {
        return (a < other.a) or (a == other.a and b < other.b);
    }

    bool operator==(const edge& other) const
    {
        return a == other.a and b == other.b;
    }
};

struct polygon {
    std::vector<point>  points;
    size_t              domain_id;

    void add_point(const point& pt) {
        points.push_back(pt);
    }

    void make_absolute(size_t base_x, size_t base_y)
    {
        for (auto& point : points)
            point.make_absolute(base_x, base_y);
    }

    bool operator<(const polygon& other) const {
        if (points.size() < other.points.size())
            return true;

        if (points.size() > other.points.size())
            return false;

        return std::lexicographical_compare(points.begin(), points.end(),
            other.points.begin(), other.points.end());
    }
};

void make_main_hex(std::vector<polygon>& local_polygons)
{
    polygon hex;
    hex.add_point( point(0,2) );
    hex.add_point( point(1,1) );
    hex.add_point( point(2,1) );
    hex.add_point( point(3,2) );
    hex.add_point( point(2,3) );
    hex.add_point( point(1,3) );
    local_polygons.push_back(hex);
}

void add_bottom(std::vector<polygon>& local_polygons,
                std::map<boundary, std::vector<edge>>& bndedges)
{
    bndedges[boundary::BOTTOM].push_back( edge(point(0,0), point(3,0)) );
    bndedges[boundary::BOTTOM].push_back( edge(point(3,0), point(4,0)) );

    polygon bottom_quad;
    bottom_quad.add_point( point(0,0) );
    bottom_quad.add_point( point(3,0) );
    bottom_quad.add_point( point(2,1) );
    bottom_quad.add_point( point(1,1) );
    local_polygons.push_back( bottom_quad );
}

void add_top(std::vector<polygon>& local_polygons,
             std::map<boundary, std::vector<edge>>& bndedges)
{
    bndedges[boundary::TOP].push_back( edge(point(0,4), point(3,4)) );
    bndedges[boundary::TOP].push_back( edge(point(3,4), point(4,4)) );

    polygon top_quad;
    top_quad.add_point( point(0,4) );
    top_quad.add_point( point(1,3) );
    top_quad.add_point( point(2,3) );
    top_quad.add_point( point(3,4) );
    local_polygons.push_back( top_quad );
}

void add_top(std::vector<polygon>& local_polygons)
{
    polygon top_hex;
    top_hex.add_point( point(0,4) );
    top_hex.add_point( point(1,3) );
    top_hex.add_point( point(2,3) );
    top_hex.add_point( point(3,4) );
    top_hex.add_point( point(2,5) ); // on the neighboring tile
    top_hex.add_point( point(1,5) ); // on the neighboring tile
    local_polygons.push_back( top_hex );
}

void add_left(std::vector<polygon>& local_polygons,
              std::map<boundary, std::vector<edge>>& bndedges)
{
    bndedges[boundary::LEFT].push_back( edge(point(0,0), point(0,2)) );
    bndedges[boundary::LEFT].push_back( edge(point(0,2), point(0,4)) );

    polygon tri_lower;
    tri_lower.add_point( point(0,0) );
    tri_lower.add_point( point(1,1) );
    tri_lower.add_point( point(0,2) );
    local_polygons.push_back(tri_lower);

    polygon tri_upper;
    tri_upper.add_point( point(0,2) );
    tri_upper.add_point( point(1,3) );
    tri_upper.add_point( point(0,4) );
    local_polygons.push_back(tri_upper);
}

void add_right(std::vector<polygon>& local_polygons,
               std::map<boundary, std::vector<edge>>& bndedges)
{
    bndedges[boundary::RIGHT].push_back( edge(point(4,0), point(4,2)) );
    bndedges[boundary::RIGHT].push_back( edge(point(4,2), point(4,4)) );

    polygon penta_lower;
    penta_lower.add_point( point(4,0) );
    penta_lower.add_point( point(4,2) );
    penta_lower.add_point( point(3,2) );
    penta_lower.add_point( point(2,1) );
    penta_lower.add_point( point(3,0) );
    local_polygons.push_back(penta_lower);

    polygon penta_upper;
    penta_upper.add_point( point(4,2) );
    penta_upper.add_point( point(4,4) );
    penta_upper.add_point( point(3,4) );
    penta_upper.add_point( point(2,3) );
    penta_upper.add_point( point(3,2) );
    local_polygons.push_back(penta_upper);
}

void add_right(std::vector<polygon>& local_polygons)
{
    polygon hex_lower;
    hex_lower.add_point( point(4,0) );
    hex_lower.add_point( point(5,1) ); // on the neighboring tile
    hex_lower.add_point( point(4,2) );
    hex_lower.add_point( point(3,2) );
    hex_lower.add_point( point(2,1) );
    hex_lower.add_point( point(3,0) );
    local_polygons.push_back(hex_lower);

    polygon hex_upper;
    hex_upper.add_point( point(4,2) );
    hex_upper.add_point( point(5,3) ); // on the neighboring tile
    hex_upper.add_point( point(4,4) );
    hex_upper.add_point( point(3,4) );
    hex_upper.add_point( point(2,3) );
    hex_upper.add_point( point(3,2) );
    local_polygons.push_back(hex_upper);
}

struct hexmesh {
    size_t                                  coord_max;
    std::vector<point>                      points;
    std::vector<edge>                       edges;
    std::vector<polygon>                    polygons;
    std::map<boundary, std::vector<edge>>   boundary_edges;

    using optn = std::optional<size_t>;
    using npair = std::pair<optn, optn>;
    std::map<edge, npair>                   edge_owners;
};

void
update_geometry(hexmesh& hm, const std::vector<polygon>& polys,
    std::map<boundary, std::vector<edge>>& bndedges)
{
    hm.polygons.insert(hm.polygons.end(), polys.begin(), polys.end());
    for (const auto& [boundary, edges] : bndedges) {
        auto end = hm.boundary_edges[boundary].end();
        hm.boundary_edges[boundary].insert(end, edges.begin(), edges.end());
    }

    for (const auto& poly : polys)
    {
        const auto& points = poly.points;
        for (size_t i = 0; i < points.size(); i++)
        {
            const auto& pa = points[i];
            const auto& pb = points[(i+1)%points.size()];
            hm.edges.push_back( edge(pa, pb) );
            hm.points.push_back(pa);
        }
    }
}

void hexmesh_work(hexmesh& hm, size_t xmin, size_t xmax, size_t ymin, size_t ymax)
{
    assert( (xmax - xmin) == (ymax - ymin) );
    if ( (xmax - xmin == 4) and (ymax - ymin == 4) )
    {
        bool bottom = (ymin == 0);
        bool top = (ymax == hm.coord_max);
        bool left = (xmin == 0);
        bool right = (xmax == hm.coord_max);

        std::vector<polygon>                    local_polygons;
        std::map<boundary, std::vector<edge>>   local_boundary_edges;

        make_main_hex(local_polygons);

        if (bottom)
            add_bottom(local_polygons, local_boundary_edges);

        if (top)
            add_top(local_polygons, local_boundary_edges);
        else
            add_top(local_polygons);

        if (left)
            add_left(local_polygons, local_boundary_edges);

        if (right)
            add_right(local_polygons, local_boundary_edges);
        else
            add_right(local_polygons);

        for(auto& polygon : local_polygons)
            polygon.make_absolute(xmin, ymin);

        for(auto& [boundary,edges] : local_boundary_edges)
            for (auto& edge : edges)
                edge.make_absolute(xmin, ymin);

        update_geometry(hm, local_polygons, local_boundary_edges);

        return;
    }

    hexmesh_work(hm, xmin, (xmax+xmin)/2, ymin, (ymax+ymin)/2);
    hexmesh_work(hm, (xmax+xmin)/2, xmax, ymin, (ymax+ymin)/2);
    hexmesh_work(hm, xmin, (xmax+xmin)/2, (ymax+ymin)/2, ymax);
    hexmesh_work(hm, (xmax+xmin)/2, xmax, (ymax+ymin)/2, ymax);
}

size_t ipow(size_t base, size_t exp)
{
    size_t ret = 1;
    while (exp--)
        ret *= base;
    return ret;
}

void make_hexmesh(hexmesh& hm, size_t l)
{
    std::cout << "Meshing" << std::endl;
    hm.coord_max = ipow(2,l+2);
    std::cout << "coord_max: " << hm.coord_max << std::endl;
    hexmesh_work(hm, 0, hm.coord_max, 0, hm.coord_max);

    std::cout << "Sorting" << std::endl;
    std::sort(hm.points.begin(), hm.points.end());
    hm.points.erase(
        std::unique(hm.points.begin(), hm.points.end()), hm.points.end() );

    std::sort(hm.edges.begin(), hm.edges.end());
    hm.edges.erase(
        std::unique(hm.edges.begin(), hm.edges.end()), hm.edges.end() );

    std::sort(hm.polygons.begin(), hm.polygons.end());

    auto add_edge_owner = [&](const edge& edge, size_t polyofs) {
        if (hm.edge_owners[edge].first)
        {
            if (hm.edge_owners[edge].second)
                throw std::logic_error("more than two owners for a face.");
            else
                hm.edge_owners[edge].second = polyofs;
        }
        else
            hm.edge_owners[edge].first = polyofs;
    };

    for (size_t polyofs = 0; polyofs < hm.polygons.size(); polyofs++)
    {
        const auto& polygon = hm.polygons[polyofs];
        const auto& points = polygon.points;
        for (size_t i = 0; i < points.size(); i++)
        {
            const auto& pa = points[i];
            const auto& pb = points[(i+1)%points.size()];
            add_edge_owner( edge(pa, pb), polyofs);
        }
    }

    std::cout << "Points: " << hm.points.size() << ", Edges: ";
    std::cout << hm.edges.size() << ", Polygons: " << hm.polygons.size();
    std::cout << std::endl;
    for (auto& [boundary,edges] : hm.boundary_edges)
        std::cout << "Boundary: " << edges.size() << std::endl;
}

size_t
offset(const hexmesh& hm, const point& pt, bool one_based)
{
    auto ipt = std::lower_bound(hm.points.begin(), hm.points.end(), pt);
    if (ipt == hm.points.end() or not(*ipt == pt))
        throw std::logic_error("point not found, this is a bug.");

    if (one_based)
        return std::distance(hm.points.begin(), ipt) + 1;

    return std::distance(hm.points.begin(), ipt);
}


} //namespace disk::fvca5
