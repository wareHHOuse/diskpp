/*
 * mesh_tagging.hpp -- assign boundary ids by geometry.
 *
 * Netgen meshes carry markers from the mesh file, but the polytopal formats
 * (FVCA5, FVCA6, poly) do not: none of the loaders writes boundary_id, so
 * every boundary face comes back with id 0.  Since the boundary_conditions
 * object matches faces by msh.boundary_id(face_id) == b_id, the ids have to
 * be written into the mesh before any BC is declared.
 *
 * DiSk++ has mark_boundary_faces() and make_interface_boundaries(), but the
 * first only sets the is_boundary flag (topology, not ids) and the second
 * tags interfaces between subdomains.  Tagging by position is the missing
 * piece; this is the analogue of Gridap's face_labeling_from_vertex_filter.
 *
 * Usage (unit square, the usual 1=top 2=right 3=bottom 4=left convention):
 *
 *     disk::generic_mesh<T,2> msh;
 *     disk::load_mesh_fvca5_2d<T>("hex.typ1", msh);
 *
 *     tag_boundary_plane(msh, 1, 1, 1.0);   // y = 1
 *     tag_boundary_plane(msh, 2, 0, 1.0);   // x = 1
 *     tag_boundary_plane(msh, 3, 1, 0.0);   // y = 0
 *     tag_boundary_plane(msh, 4, 0, 0.0);   // x = 0
 *     report_boundary_tags(msh);
 *
 * The predicate is evaluated at the face BARYCENTRE, so a face is tagged
 * only if it lies on the plane as a whole (up to the tolerance).
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include <cassert>
#include <array>

#include "diskpp/mesh/mesh.hpp"

namespace disk
{

/* Tag every boundary face whose barycentre satisfies pred.  Returns how
 * many faces were tagged: the cheapest way to catch a predicate that
 * silently matches nothing.
 *
 * offset(msh, fc) is the library's own way of indexing boundary_info, so
 * this makes no assumption about face iteration order matching face ids.
 * The three-argument boundary_descriptor keeps id and tag equal and leaves
 * is_internal false, which is what a single-domain mesh wants.            */
template<typename Mesh, typename Predicate>
size_t
tag_boundary(Mesh& msh, const size_t bnd_id, Predicate pred)
{
    auto storage = msh.backend_storage();
    assert(storage->boundary_info.size() == msh.faces_size());

    size_t tagged = 0;

    for (auto& fc : faces(msh))
    {
        const auto ofs = offset(msh, fc);

        if (!storage->boundary_info.at(ofs).is_boundary())
            continue;

        if (pred(barycenter(msh, fc)))
        {
            storage->boundary_info.at(ofs) =
                boundary_descriptor(bnd_id, bnd_id, true);
            tagged++;
        }
    }

    return tagged;
}

/* Axis-aligned plane: coordinate `axis` equal to `value`, with
 * axis = 0 (x), 1 (y), 2 (z).  This covers the box-shaped domains of every
 * benchmark here.
 *
 * Coordinates are read through pt.at(axis) rather than pt.x()/y()/z():
 * point<T,2> has no z(), so a switch over the three accessors fails to
 * compile in 2D even when the z branch is unreachable at run time.        */
template<typename Mesh>
size_t
tag_boundary_plane(Mesh&                                msh,
                   const size_t                         bnd_id,
                   const size_t                         axis,
                   const typename Mesh::coordinate_type value,
                   const typename Mesh::coordinate_type tol = 1e-8)
{
    assert(axis < Mesh::dimension);

    return tag_boundary(msh, bnd_id, [axis, value, tol](const auto& pt) {
        return std::abs(pt.at(axis) - value) < tol;
    });
}


/* ==========================================================================
 * Geometry transformation
 *
 * The built-in generators produce the unit square or the unit cube.  Moving
 * the points is enough to obtain any box: faces and cells store point ids,
 * and every geometric quantity (barycenter, measure, normal, diameter) is
 * computed from the points on demand, so nothing needs rebuilding.
 *
 *     scale_mesh(msh, {8.0, 4.0, 1.0});     // unit cube -> (0,8)x(0,4)x(0,1)
 * ========================================================================*/

template<typename Mesh, typename Map>
void
transform_mesh(Mesh& msh, Map map)
{
    auto storage = msh.backend_storage();
    for (auto& pt : storage->points)
        pt = map(pt);
}

template<typename Mesh>
void
scale_mesh(Mesh& msh,
           const std::array<typename Mesh::coordinate_type,
                            Mesh::dimension>& factors)
{
    using point_type = typename Mesh::point_type;

    transform_mesh(msh, [&factors](const point_type& pt) {
        point_type q = pt;
        for (size_t d = 0; d < Mesh::dimension; d++)
            q.at(d) = pt.at(d) * factors[d];
        return q;
    });
}

/* Print how many boundary faces carry each id.  Run once after tagging: an
 * id with zero faces, or a leftover id 0, is a tagging bug that would
 * otherwise surface much later as a silently wrong problem.               */
template<typename Mesh>
void
report_boundary_tags(const Mesh& msh)
{
    auto storage = msh.backend_storage();

    std::map<size_t, size_t> counts;
    size_t n_boundary = 0;

    for (size_t i = 0; i < msh.faces_size(); i++)
    {
        const auto& bi = storage->boundary_info.at(i);
        if (!bi.is_boundary())
            continue;
        counts[bi.id()]++;
        n_boundary++;
    }

    std::cout << "  boundary faces: " << n_boundary << "   ";
    for (auto& [id, n] : counts)
        std::cout << "id " << id << ": " << n << "   ";
    std::cout << "\n";

    auto untagged = counts.find(0);
    if (untagged != counts.end())
        std::cout << "  *** " << untagged->second
                  << " boundary faces still carry id 0 (untagged) ***\n";
}

} // namespace disk