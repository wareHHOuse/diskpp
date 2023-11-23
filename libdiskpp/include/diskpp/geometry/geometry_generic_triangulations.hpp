#pragma once

#include "diskpp/common/simplicial_formula.hpp"

#include "triangle_mesher.h"

namespace disk {

template<typename T, size_t DIM>
struct triangle {
    point<T,DIM>    p0;
    point<T,DIM>    p1;
    point<T,DIM>    p2;
};

template<typename T, size_t DIM>
auto
barycenter(const triangle<T,DIM>& t)
{
    return (t.p0 + t.p1 + t.p2)/3.;
}

template<typename T, size_t DIM>
auto
measure(const triangle<T,DIM>& t)
{
    return area_triangle_kahan(t.p0, t.p1, t.p2);   
}

/* Call J. R. Shewchuk's Triangle to triangulate a mesh element */
inline
std::vector<triangle<double,2>>
triangulate_nonconvex_polygon(const generic_mesh<double,2>& msh,
    const typename generic_mesh<double,2>::cell_type& cl)
{
    auto pts = points(msh, cl);
    std::vector<double> tri_pts;
    tri_pts.reserve(2*pts.size());
    for (auto& pt : pts) {
        tri_pts.push_back(pt.x());
        tri_pts.push_back(pt.y());
    }

    std::vector<int> seglist;
    for (size_t i = 0; i < pts.size(); i++) {
        seglist.push_back(i);
        seglist.push_back((i+1)%pts.size());
    }

    struct triangulateio tio_in;
    memset(&tio_in, '\0', sizeof(struct triangulateio));
    tio_in.pointlist = tri_pts.data();
    tio_in.numberofpoints = pts.size();
    tio_in.segmentlist = seglist.data();
    tio_in.numberofsegments = pts.size();

    struct triangulateio tio_out;
    memset(&tio_out, '\0', sizeof(struct triangulateio));

    triangulate("zpQ", &tio_in, &tio_out, nullptr);

    std::vector<triangle<double,2>> ret;
    for (int i = 0; i < tio_out.numberoftriangles; i++)
    {
        int nbase = 3*i;
        int p0base = tio_out.trianglelist[nbase+0];
        assert(p0base < tio_out.numberofpoints);
        assert(2*p0base+1 < 2*tio_out.numberofpoints);
        int p1base = tio_out.trianglelist[nbase+1];
        assert(p1base < tio_out.numberofpoints);
        assert(2*p1base+1 < 2*tio_out.numberofpoints);
        int p2base = tio_out.trianglelist[nbase+2];
        assert(p2base < tio_out.numberofpoints);
        assert(2*p2base+1 < 2*tio_out.numberofpoints);
        triangle<double, 2> t;
        t.p0 = point<double,2>( tio_out.pointlist[2*p0base+0], tio_out.pointlist[2*p0base+1] );
        t.p1 = point<double,2>( tio_out.pointlist[2*p1base+0], tio_out.pointlist[2*p1base+1] );
        t.p2 = point<double,2>( tio_out.pointlist[2*p2base+0], tio_out.pointlist[2*p2base+1] );
        
        ret.push_back(t);
    }

    if (tio_out.pointlist) trifree(tio_out.pointlist);
    if (tio_out.pointattributelist) trifree(tio_out.pointattributelist);
    if (tio_out.pointmarkerlist) trifree(tio_out.pointmarkerlist);
    if (tio_out.trianglelist) trifree(tio_out.trianglelist);
    if (tio_out.triangleattributelist) trifree(tio_out.triangleattributelist);
    if (tio_out.trianglearealist) trifree(tio_out.trianglearealist);
    if (tio_out.neighborlist) trifree(tio_out.neighborlist);
    if (tio_out.segmentlist) trifree(tio_out.segmentlist);
    if (tio_out.segmentmarkerlist) trifree(tio_out.segmentmarkerlist);
    if (tio_out.holelist) trifree(tio_out.holelist);
    if (tio_out.regionlist) trifree(tio_out.regionlist);
    if (tio_out.edgelist) trifree(tio_out.edgelist);
    if (tio_out.edgemarkerlist) trifree(tio_out.edgemarkerlist);
    if (tio_out.normlist) trifree(tio_out.normlist);

    return ret;
}

template<typename T>
auto
triangulate_convex_polygon(const generic_mesh<T,2>& msh,
    const typename generic_mesh<T,2>::cell_type& cl)
{
    std::vector<triangle<T,2>> ret;

    auto pts = points(msh, cl);
    assert(pts.size() > 2);
    auto center = std::accumulate(pts.begin(), pts.end(), point<T,2>(0,0));
    center = center/T(pts.size());

    for (size_t i = 0; i < pts.size(); i++)
    {
        triangle<T,2> t;
        t.p0 = pts[i];
        t.p1 = pts[(i+1)%pts.size()];
        t.p2 = center;
        ret.push_back(t);
    }

    return ret;
}

template<typename T>
auto
triangulate_polygon(const generic_mesh<T,2>& msh,
    const typename generic_mesh<T,2>::cell_type& cl)
{
    if ( is_convex(msh, cl) )
        return triangulate_convex_polygon(msh, cl);

    return triangulate_nonconvex_polygon(msh, cl);
}

} // namespace disk