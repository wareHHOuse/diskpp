/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */
#include "core/output/hdf5_io.hpp"

namespace disk {

namespace mesh_v2 {

template<size_t DIM>
struct simplicial_storage_class;

template<>
struct simplicial_storage_class<1> {
    typedef edge<0,void>        edge_type;
    typedef node<1,void>        node_type;
};

template<>
struct simplicial_storage_class<2> {
    typedef triangle<0,void>    surface_type;
    typedef edge<1,void>        edge_type;
    typedef node<2,void>        node_type;
};

template<typename T, size_t DIM>
using simplicial_mesh_storage = mesh_storage<T, DIM, simplicial_storage_class<DIM>>;

template<typename T, size_t DIM>
using simplicial_mesh = mesh<DIM, simplicial_mesh_storage<T, DIM>>;

template<typename T>
using triangular_mesh = simplicial_mesh<T, 2>;


template<typename T>
bool save(const hdf5_context& hctx, const triangular_mesh<T>& mesh,
          const std::string& name)
{
    hid_t file_id = hctx.get_file_descriptor();

    hid_t group_id;

    herr_t status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    status = H5Gget_objinfo (file_id, "/meshes/standalone", 0, NULL);

    if (status != 0)
    {
        group_id = H5Gcreate(file_id, "meshes",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        group_id = H5Gcreate(group_id, "standalone",
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        group_id = H5Gopen(file_id, "/meshes/standalone", H5P_DEFAULT);
    }

    group_id = H5Gcreate(group_id, name.c_str(),
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //if (group_id < 0)
    //{
    //    group_id = H5Gcreate(file_id, "/meshes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //    std::cout << "GID: " << group_id << std::endl;
    //}

    return true;
}


struct child_cells_offsets
{
    std::array<size_t, 4>   offsets;
    bool                    has_childs;

    child_cells_offsets() : has_childs(false) {}
};

template<size_t DIM>
struct hierarchical_simplicial_storage_class;

template<>
struct hierarchical_simplicial_storage_class<2>
{
    typedef triangle<0, child_cells_offsets>    surface_type;
    typedef edge<1, void>                       edge_type;
    typedef node<2, void>                       node_type;
};

template<typename T, size_t DIM>
using hierarchical_simplicial_mesh_storage = mesh_storage<T, DIM,
                                    hierarchical_simplicial_storage_class<DIM>>;

template<typename T, size_t DIM>
using hierarchical_simplicial_mesh = mesh<DIM, hierarchical_simplicial_mesh_storage<T, DIM>>;

template<typename T>
using hierarchical_triangular_mesh = hierarchical_simplicial_mesh<T, 2>;



/* Return the point identifiers of a specified element. This is a catch-all function,
 * but it is likely that some elements will not have the point_identifiers() method.
 * This could happen for example in three-dimensional arbitrary polytopes, where is
 * not possible to establish a clear association between number of vertices and
 * shape.
 * In those cases this function has to be specialized accordingly.
 */
template<typename Mesh, typename Element>
auto
point_identifiers(const Mesh& msh, const Element& elem)
{
    return elem.point_identifiers();
}

template<typename Mesh, size_t EDIM, size_t CODIM, typename UserData, size_t N>
auto
points(const Mesh& msh, const polytope<EDIM, CODIM, UserData, fixed_storage_polytope<N>>& poly)
{
    typedef polytope<EDIM, CODIM, UserData, fixed_storage_polytope<N>> element_type;

    auto ptids = point_identifiers(msh, poly);
    assert( ptids.size() == N );
    std::array<typename Mesh::point_type, N> ret;

    for (size_t i = 0; i < ptids.size(); i++)
        ret[i] = *std::next(msh.points_begin(), ptids[i]);

    return ret;
}

/* Computes the measure of a triangle. Assumes that points are stored in CCW order.
 */
template<typename Mesh, size_t CODIM, typename UserData>
typename Mesh::coordinate_type
measure(const Mesh& msh, const triangle<CODIM, UserData>& t)
{
    static_assert(Mesh::dimension == 2 || Mesh::dimension == 3,
                  "Wrong mesh dimension for a triangle");

    auto pts = points(msh, t);
    assert(pts.size() == 3);

    auto v1 = (pts[1] - pts[0]).to_vector();
    auto v2 = (pts[2] - pts[1]).to_vector();

    if (Mesh::dimension == 2)
        return std::abs( (v1(0)*v2(1) - v2(0)*v1(1)) / 2.0 );
    else if (Mesh::dimension > 2)
        return v1.norm() * v2.norm() / 2.0;
    else
        throw std::logic_error("invalid dimension");
}

template<typename Mesh, size_t CODIM, typename UserData>
typename Mesh::coordinate_type
measure(const Mesh& msh, const edge<CODIM, UserData>& e)
{
    static_assert(Mesh::dimension > 0, "Wrong mesh dimension for an edge");

    auto pts = points(msh, e);
    assert(pts.size() == 2);

    return (pts[1] - pts[0]).to_vector().norm();
}

template<typename Mesh, size_t CODIM, typename UserData>
typename Mesh::point_type
barycenter(const Mesh& msh, const triangle<CODIM, UserData>& t)
{
    static_assert(Mesh::dimension == 2 || Mesh::dimension == 3,
                  "Wrong mesh dimension for a triangle");

    auto pts = points(msh, t);
    assert(pts.size() == 3);

    return (pts[0] + pts[1] + pts[2])/3.0;
}

/* Return the faces of an element. This is the case for fixed storage, so we return
 * a std::array.
 */
template<typename Mesh>
std::array<typename Mesh::face_type, polytope_traits<typename Mesh::cell_type>::storage_size>
faces(const Mesh& msh,
      typename std::enable_if<is_fixed_storage<typename Mesh::cell_type>::value,
                              const typename Mesh::cell_type>::type& cl)
{
    using return_type = std::array<typename Mesh::face_type,
                                   polytope_traits<typename Mesh::cell_type>::storage_size>;

    return_type ret;
    throw std::logic_error("You should not be here");

    return ret;
}

template<typename Mesh>
std::vector<typename Mesh::face_type>
faces(const Mesh& msh,
        typename std::enable_if<!is_fixed_storage<typename Mesh::cell_type>::value,
                                const typename Mesh::cell_type>::type& cl)
{
    throw std::logic_error("You should not be here");
}

} //namespace mesh_v2

} //namespace disk
