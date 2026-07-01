/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2016-2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 *
 * This file is copyright of the following author:
 * Guillaume Delay (C) 2020-2021         guillaume.delay@sorbonne-universite.fr
 * Sorbonne Universite
 * Laboratoire Jacques-Louis Lions (LJLL)
 *
 * Matteo Cicuttin (C) 2016-2026
 * matteo.cicuttin@polito.it
 *
 */
/*
 * The content of this file corresponds to the implementation of Lagrange basis functions.
 */

#include "diskpp/bases/bases.hpp"

//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   LAGRANGE BASES    ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

using namespace disk;

/* Generic template for Lagrange bases. */
template<typename MeshType, typename Element>
struct Lagrange_scalar_basis
{
    static_assert(sizeof(MeshType) == -1, "Lagrange_scalar_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "Lagrange_scalar_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_scalar_Lagrange_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
    return Lagrange_scalar_basis<MeshType, ElementType>(msh, elem, degree);
}


/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Lagrange_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 2>     gradient_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>     function_type;

  private:
    std::vector<point_type>       vertices;
    size_t                        basis_degree, basis_size;

#ifdef POWER_CACHE
    mutable std::vector<scalar_type> power_cache;
#endif

  public:
    Lagrange_scalar_basis(const mesh_type& msh, const cell_type& cl, size_t degree)
    {
        if( degree > 4 )
            throw std::invalid_argument("degree > 4 not yet supported");
        basis_degree = degree;

        // store the vertices
        auto pts = points(msh, cl);
        assert( pts.size() == 3);
        vertices.push_back( pts[0] );
        vertices.push_back( pts[1] );
        vertices.push_back( pts[2] );
        // vertices = points(msh, cl);

        basis_size = scalar_basis_size(degree, 2);
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        if(basis_degree == 0)
            ret(0) = 1.0;
        else if(basis_degree == 1)
        {
            ret = bar_coord(pt);
        }
        else if(basis_degree == 2)
        {
            auto bar_c = bar_coord(pt);
            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
                ret[i] = bar_c[i] * (2.0*bar_c[i] - 1.0);

            // 3-5 : mid-points (12,02,01)
            ret[3] = 4.0*bar_c[1]*bar_c[2];
            ret[4] = 4.0*bar_c[0]*bar_c[2];
            ret[5] = 4.0*bar_c[0]*bar_c[1];
        }
        else if(basis_degree == 3)
        {
            auto bar_c = bar_coord(pt);
            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
                ret[i] = 0.5*bar_c[i] * (3.0*bar_c[i] - 1.0) * (3.0*bar_c[i] - 2.0);

            // 3-8 : face-points (112,122,002,022,001,011)
            ret[3] = 4.5*bar_c[1]*bar_c[2]*(3.0*bar_c[1]-1.0);
            ret[4] = 4.5*bar_c[1]*bar_c[2]*(3.0*bar_c[2]-1.0);
            ret[5] = 4.5*bar_c[0]*bar_c[2]*(3.0*bar_c[0]-1.0);
            ret[6] = 4.5*bar_c[0]*bar_c[2]*(3.0*bar_c[2]-1.0);
            ret[7] = 4.5*bar_c[0]*bar_c[1]*(3.0*bar_c[0]-1.0);
            ret[8] = 4.5*bar_c[0]*bar_c[1]*(3.0*bar_c[1]-1.0);

            // 9 : center of mass
            ret[9] = 27.0*bar_c[0]*bar_c[1]*bar_c[2];
        }
        else // degree == 4
        {
            auto bar_c = bar_coord(pt);
            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
                ret[i] = (1.0/6.0)*bar_c[i] * (4.0*bar_c[i] - 1.0)
                    * (4.0*bar_c[i] - 2.0) * (4.0*bar_c[i] - 3.0);

            // 3-5 : mid-points (12,02,01)
            ret[3] = 4.0*bar_c[1]*bar_c[2]*(4.0*bar_c[1]-1.0)*(4.0*bar_c[2]-1.0);
            ret[4] = 4.0*bar_c[0]*bar_c[2]*(4.0*bar_c[0]-1.0)*(4.0*bar_c[2]-1.0);
            ret[5] = 4.0*bar_c[0]*bar_c[1]*(4.0*bar_c[0]-1.0)*(4.0*bar_c[1]-1.0);

            // 6-11 : face-points (1112,1222,0002,0222,0001,0111)
            ret[6] = (8.0/3.0)*bar_c[1]*bar_c[2]*(4.0*bar_c[1]-1.0)*(4.0*bar_c[1]-2.0);
            ret[7] = (8.0/3.0)*bar_c[2]*bar_c[1]*(4.0*bar_c[2]-1.0)*(4.0*bar_c[2]-2.0);
            ret[8] = (8.0/3.0)*bar_c[0]*bar_c[2]*(4.0*bar_c[0]-1.0)*(4.0*bar_c[0]-2.0);
            ret[9] = (8.0/3.0)*bar_c[2]*bar_c[0]*(4.0*bar_c[2]-1.0)*(4.0*bar_c[2]-2.0);
            ret[10] = (8.0/3.0)*bar_c[0]*bar_c[1]*(4.0*bar_c[0]-1.0)*(4.0*bar_c[0]-2.0);
            ret[11] = (8.0/3.0)*bar_c[1]*bar_c[0]*(4.0*bar_c[1]-1.0)*(4.0*bar_c[1]-2.0);

            // 12-14 : others (0012,0112,0122)
            ret[12] = 32.0 * bar_c[0]*bar_c[1]*bar_c[2]*(4.0*bar_c[0]-1.0);
            ret[13] = 32.0 * bar_c[0]*bar_c[1]*bar_c[2]*(4.0*bar_c[1]-1.0);
            ret[14] = 32.0 * bar_c[0]*bar_c[1]*bar_c[2]*(4.0*bar_c[2]-1.0);
        }
        return ret;
    }

    gradient_type
    eval_gradients(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 2);

        if(basis_degree == 0)
        {
            ret(0,0) = 0.0;
            ret(0,1) = 0.0;
        }
        else if(basis_degree == 1)
        {
            ret = bar_coord_grad(pt);
        }
        else if(basis_degree == 2)
        {
            auto bar_c = bar_coord(pt);
            auto bar_c_g = bar_coord_grad(pt);

            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
            {
                ret(i,0) = (4.0*bar_c[i]-1.0) * bar_c_g(i,0);
                ret(i,1) = (4.0*bar_c[i]-1.0) * bar_c_g(i,1);
            }
            // 3-5 : mid-points (12,02,01)
            for(size_t j = 0; j < 2; j++)
            {
                ret(3,j) = 4.0*( bar_c[1]*bar_c_g(2,j) + bar_c[2]*bar_c_g(1,j) );
                ret(4,j) = 4.0*( bar_c[0]*bar_c_g(2,j) + bar_c[2]*bar_c_g(0,j) );
                ret(5,j) = 4.0*( bar_c[0]*bar_c_g(1,j) + bar_c[1]*bar_c_g(0,j) );
            }
        }
        else if(basis_degree == 3)
        {
            auto bar_c = bar_coord(pt);
            auto bar_c_g = bar_coord_grad(pt);

            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
            {
                T coeff = 0.5 * (27.0*bar_c[i]*bar_c[i] - 18.0 * bar_c[i] + 2.0);
                ret(i,0) = coeff * bar_c_g(i,0);
                ret(i,1) = coeff * bar_c_g(i,1);
            }
            // 3-8 : face-points (112,122,002,022,001,011)
            for(size_t j = 0; j < 2; j++)
            {
                ret(3,j) = 4.5*( (6.0*bar_c[1]-1.0)*bar_c[2]*bar_c_g(1,j)
                                 + (3.0*bar_c[1]-1.0)*bar_c[1]*bar_c_g(2,j) );
                ret(4,j) = 4.5*( (6.0*bar_c[2]-1.0)*bar_c[1]*bar_c_g(2,j)
                                 + (3.0*bar_c[2]-1.0)*bar_c[2]*bar_c_g(1,j) );
                ret(5,j) = 4.5*( (6.0*bar_c[0]-1.0)*bar_c[2]*bar_c_g(0,j)
                                 + (3.0*bar_c[0]-1.0)*bar_c[0]*bar_c_g(2,j) );
                ret(6,j) = 4.5*( (6.0*bar_c[2]-1.0)*bar_c[0]*bar_c_g(2,j)
                                 + (3.0*bar_c[2]-1.0)*bar_c[2]*bar_c_g(0,j) );
                ret(7,j) = 4.5*( (6.0*bar_c[0]-1.0)*bar_c[1]*bar_c_g(0,j)
                                 + (3.0*bar_c[0]-1.0)*bar_c[0]*bar_c_g(1,j) );
                ret(8,j) = 4.5*( (6.0*bar_c[1]-1.0)*bar_c[0]*bar_c_g(1,j)
                                 + (3.0*bar_c[1]-1.0)*bar_c[1]*bar_c_g(0,j) );
            }
            // 9 : center of mass
            for(size_t j = 0; j < 2; j++)
            {
                ret(9,j) = 27.0 * (bar_c[1] * bar_c[2] * bar_c_g(0,j)
                                   + bar_c[0] * bar_c[1] * bar_c_g(2,j)
                                   + bar_c[2] * bar_c[0] * bar_c_g(1,j) );
            }
        }
        else // degree == 4
        {
            auto bar_c = bar_coord(pt);
            auto bar_c_g = bar_coord_grad(pt);
            // 0-2 : vertices
            for(size_t i = 0; i < 3; i++)
            {
                T coeff = (1.0/6.0) * (256.0*bar_c[i]*bar_c[i]*bar_c[i] - 288.0*bar_c[i]*bar_c[i]
                                       + 88.0*bar_c[i] - 6.0);
                ret(i,0) = coeff * bar_c_g(i,0);
                ret(i,1) = coeff * bar_c_g(i,1);
            }
            // 3-5 : mid-points (12,02,01)
            for(size_t j = 0; j < 2; j++)
            {
                ret(3,j) =
                    4.0*bar_c[1]*(32.0*bar_c[1]*bar_c[2] - 4.0*(2.0*bar_c[2]+bar_c[1]) + 1.0)*bar_c_g(2,j)
                    + 4.0*bar_c[2]*(32.0*bar_c[2]*bar_c[1] - 4.0*(2.0*bar_c[1]+bar_c[2])+1.0)*bar_c_g(1,j);

                ret(4,j) =
                    4.0*bar_c[0]*(32.0*bar_c[0]*bar_c[2] - 4.0*(2.0*bar_c[2]+bar_c[0]) + 1.0)*bar_c_g(2,j)
                    + 4.0*bar_c[2]*(32.0*bar_c[2]*bar_c[0] - 4.0*(2.0*bar_c[0]+bar_c[2])+1.0)*bar_c_g(0,j);

                ret(5,j) =
                    4.0*bar_c[0]*(32.0*bar_c[0]*bar_c[1] - 4.0*(2.0*bar_c[1]+bar_c[0]) + 1.0)*bar_c_g(1,j)
                    + 4.0*bar_c[1]*(32.0*bar_c[1]*bar_c[0] - 4.0*(2.0*bar_c[0]+bar_c[1])+1.0)*bar_c_g(0,j);
            }
            // 6-11 : face-points (1112,1222,0002,0222,0001,0111)
            for(size_t j = 0; j < 2; j++)
            {
                ret(6,j) =
                    (16.0/3.0) * bar_c[2] * (24.0*bar_c[1]*bar_c[1] -12.0*bar_c[1] + 1.0)*bar_c_g(1,j)
                    +(16.0/3.0) * bar_c[1] * (8.0*bar_c[1]*bar_c[1] -6.0*bar_c[1] + 1.0)*bar_c_g(2,j);

                ret(7,j) =
                    (16.0/3.0) * bar_c[1] * (24.0*bar_c[2]*bar_c[2] -12.0*bar_c[2] + 1.0)*bar_c_g(2,j)
                    +(16.0/3.0) * bar_c[2] * (8.0*bar_c[2]*bar_c[2] -6.0*bar_c[2] + 1.0)*bar_c_g(1,j);

                ret(8,j) =
                    (16.0/3.0) * bar_c[2] * (24.0*bar_c[0]*bar_c[0] -12.0*bar_c[0] + 1.0)*bar_c_g(0,j)
                    +(16.0/3.0) * bar_c[0] * (8.0*bar_c[0]*bar_c[0] -6.0*bar_c[0] + 1.0)*bar_c_g(2,j);

                ret(9,j) =
                    (16.0/3.0) * bar_c[0] * (24.0*bar_c[2]*bar_c[2] -12.0*bar_c[2] + 1.0)*bar_c_g(2,j)
                    +(16.0/3.0) * bar_c[2] * (8.0*bar_c[2]*bar_c[2] -6.0*bar_c[2] + 1.0)*bar_c_g(0,j);

                ret(10,j) =
                    (16.0/3.0) * bar_c[1] * (24.0*bar_c[0]*bar_c[0] -12.0*bar_c[0] + 1.0)*bar_c_g(0,j)
                    +(16.0/3.0) * bar_c[0] * (8.0*bar_c[0]*bar_c[0] -6.0*bar_c[0] + 1.0)*bar_c_g(1,j);

                ret(11,j) =
                    (16.0/3.0) * bar_c[0] * (24.0*bar_c[1]*bar_c[1] -12.0*bar_c[1] + 1.0)*bar_c_g(1,j)
                    +(16.0/3.0) * bar_c[1] * (8.0*bar_c[1]*bar_c[1] -6.0*bar_c[1] + 1.0)*bar_c_g(0,j);
            }
            // 12-14 : others (0012,0112,0122)
            for(size_t j = 0; j < 2; j++)
            {
                ret(12,j) = 32.0*(4.0*bar_c[0]-1.0)
                    *(bar_c[0]*bar_c[1]*bar_c_g(2,j)
                      + bar_c[1]*bar_c[2]*bar_c_g(0,j) + bar_c[2]*bar_c[0]*bar_c_g(1,j))
                    + 128.0 * bar_c[0]*bar_c[1]*bar_c[2]*bar_c_g(0,j);

                ret(13,j) = 32.0*(4.0*bar_c[1]-1.0)
                    *(bar_c[0]*bar_c[1]*bar_c_g(2,j)
                      + bar_c[1]*bar_c[2]*bar_c_g(0,j) + bar_c[2]*bar_c[0]*bar_c_g(1,j))
                    + 128.0 * bar_c[0]*bar_c[1]*bar_c[2]*bar_c_g(1,j);

                ret(14,j) = 32.0*(4.0*bar_c[2]-1.0)
                    *(bar_c[0]*bar_c[1]*bar_c_g(2,j)
                      + bar_c[1]*bar_c[2]*bar_c_g(0,j) + bar_c[2]*bar_c[0]*bar_c_g(1,j))
                    + 128.0 * bar_c[0]*bar_c[1]*bar_c[2]*bar_c_g(2,j);
            }
        }
        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }


    // barycentric coordinates
    function_type
    bar_coord(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        auto pts = vertices;
        auto x0 = pts[0].x(); auto y0 = pts[0].y();
        auto x1 = pts[1].x(); auto y1 = pts[1].y();
        auto x2 = pts[2].x(); auto y2 = pts[2].y();

        auto m = (x1*y2 - y1*x2 - x0*(y2 - y1) + y0*(x2 - x1));

        ret(0) = (x1*y2 - y1*x2 - pt.x() * (y2 - y1) + pt.y() * (x2 - x1)) / m;
        ret(1) = (x2*y0 - y2*x0 + pt.x() * (y2 - y0) - pt.y() * (x2 - x0)) / m;
        ret(2) = (x0*y1 - y0*x1 - pt.x() * (y1 - y0) + pt.y() * (x1 - x0)) / m;

        return ret;
    }

    // gradients of the barycentric coordinates
    gradient_type
    bar_coord_grad(const point_type& pt) const
    {
        gradient_type ret = gradient_type::Zero(basis_size, 2);

        auto pts = vertices;
        auto x0 = pts[0].x(); auto y0 = pts[0].y();
        auto x1 = pts[1].x(); auto y1 = pts[1].y();
        auto x2 = pts[2].x(); auto y2 = pts[2].y();

        auto m = (x1*y2 - y1*x2 - x0*(y2 - y1) + y0*(x2 - x1));

        ret(0,0) = (y1 - y2) / m;
        ret(1,0) = (y2 - y0) / m;
        ret(2,0) = (y0 - y1) / m;
        ret(0,1) = (x2 - x1) / m;
        ret(1,1) = (x0 - x2) / m;
        ret(2,1) = (x1 - x0) / m;

        return ret;
    }
};

/* Specialization for 2D meshes, faces */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Lagrange_scalar_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::face>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::face            face_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>     function_type;

  private:
    point_type  face_ref, base;
    scalar_type face_h;
    size_t      basis_degree, basis_size;

#ifdef POWER_CACHE
    mutable std::vector<scalar_type> power_cache;
#endif

  public:
    Lagrange_scalar_basis(const mesh_type& msh, const face_type& fc, size_t degree)
    {
        if( degree > 3 )
            throw std::invalid_argument("degree > 3 not yet supported");
        basis_degree = degree;
        basis_size   = degree + 1;

        const auto pts = points(msh, fc);
        face_ref = pts[0];
        face_h       = diameter(msh, fc);
        base = pts[1] - pts[0];
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size);

        // pos on the face
        const auto v   = base.to_vector();
        const auto t   = (pt - face_ref).to_vector();
        const auto dot = v.dot(t);
        const auto pos = dot/(face_h*face_h); // 0 -> pts[0], 1 -> pts[1]

        if(basis_degree == 0)
            ret(0) = 1.0;
        else if(basis_degree == 1)
        {
            ret(0) = - pos + 1.0; // 0
            ret(1) = pos;         // 1
        }
        else if(basis_degree == 2)
        {
            ret(0) = 2.0*pos*pos - 3.0*pos + 1.0; // 0
            ret(1) = -4.0*pos*pos + 4.0*pos;      // 0.5
            ret(2) = 2.0*pos*pos - pos;           // 1
        }
        else  // degree == 3
        {
            ret(0) = -0.5*(3*pos-1)*(3*pos-2)*(pos-1); // 0
            ret(1) = 4.5*pos*(3*pos-2)*(pos-1);        // 1/3
            ret(2) = -4.5*pos*(3*pos-1)*(pos-1);       // 2/3
            ret(3) = 0.5*pos*(3*pos-1)*(3*pos-2);      // 1
        }
        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};

//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////    VECTOR BASES     ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/* Generic template for bases. */
template<typename MeshType, typename Element>
struct Lagrange_vector_basis
{
    static_assert(sizeof(MeshType) == -1, "Lagrange_vector_basis: not suitable for the requested kind of mesh");
    static_assert(sizeof(Element) == -1,
                  "Lagrange_vector_basis: not suitable for the requested kind of element");
};

/* Basis 'factory'. */
template<typename MeshType, typename ElementType>
auto
make_vector_Lagrange_basis(const MeshType& msh, const ElementType& elem, size_t degree)
{
    return Lagrange_vector_basis<MeshType, ElementType>(msh, elem, degree);
}


/* Specialization for 2D meshes, cells */
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
class Lagrange_vector_basis<Mesh<T, 2, Storage>, typename Mesh<T, 2, Storage>::cell>
{

  public:
    typedef Mesh<T, 2, Storage>                 mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::point_type      point_type;
    typedef Eigen::Matrix<scalar_type, 2, 2>           gradient_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 2>     function_type;
    typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>     divergence_type;

  private:
    size_t basis_degree, basis_size;

    typedef Lagrange_scalar_basis<mesh_type, cell_type>    scalar_basis_type;
    scalar_basis_type                                      scalar_basis;

  public:
    Lagrange_vector_basis(const mesh_type& msh, const cell_type& cl, size_t degree) :
      scalar_basis(msh, cl, degree)
    {
        basis_degree = degree;
        basis_size   = 2 * scalar_basis.size();
    }

    function_type
    eval_functions(const point_type& pt) const
    {
        function_type ret = function_type::Zero(basis_size, 2);

        const auto phi = scalar_basis.eval_functions(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(2 * i, 0)     = phi(i);
            ret(2 * i + 1, 1) = phi(i);
        }
        return ret;
    }

    eigen_compatible_stdvector<gradient_type>
    eval_gradients(const point_type& pt) const
    {
        eigen_compatible_stdvector<gradient_type> ret;
        ret.reserve(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            const Eigen::Matrix<scalar_type, 1, 2> dphi_i = dphi.row(i);
            gradient_type                   g;

            g        = gradient_type::Zero();
            g.row(0) = dphi_i;
            ret.push_back(g);

            g        = gradient_type::Zero();
            g.row(1) = dphi_i;
            ret.push_back(g);
        }
        assert(ret.size() == basis_size);
        return ret;
    }

    divergence_type
    eval_divergences(const point_type& pt) const
    {
        divergence_type ret = divergence_type::Zero(basis_size);

        const function_type dphi = scalar_basis.eval_gradients(pt);

        for (size_t i = 0; i < scalar_basis.size(); i++)
        {
            ret(2 * i)     = dphi(i, 0);
            ret(2 * i + 1) = dphi(i, 1);
        }

        return ret;
    }

    size_t
    size() const
    {
        return basis_size;
    }

    size_t
    degree() const
    {
        return basis_degree;
    }
};
