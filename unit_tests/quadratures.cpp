/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Karol Cascavita (2023)
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

 /* Test for 1D raw monomial bases: interpolate a sine, measure the error
  * on the interpolation and its gradient. */

#include <iostream>
#include <fstream>
#include <functional>
#include <silo.h>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/output/silo.hpp"
using namespace disk;
using namespace disk::basis;


double function_unit(const point<double, 3>& pt)
{
    return 1.0;
};

double function_3d(const point<double, 3>& pt)
{
    return pt.x()* pt.y()* pt.z();
};

double function_2d (const point<double, 2>& pt)
{
    return (1.0 - pt.x () * pt.x ()) * (1.0 - pt.y () * pt.y ());
}


template<typename Mesh>
void test_quadratures (const std::string& name, const Mesh& msh,
        const std::function<typename Mesh::coordinate_type (typename Mesh::point_type)> func,
        const size_t& degree,
        const std::vector< typename Mesh::coordinate_type>& ref)
{
    using mesh_type = Mesh;
    using T = typename mesh_type::coordinate_type;
    using point_type = typename mesh_type::point_type;

    auto passfail = [](bool pass) {
        return pass ? "[PASS]" : "[FAIL]";
        };

    auto cl = *msh.cells_begin ();
    auto fcs = faces (msh, cl);
    std::vector<bool> pass (fcs.size (), false);

    std::cout << __FILE__ << std::endl;
    size_t i = 0;
    bool pass_test = true;
    for (auto& fc : fcs)
    {
        auto res = integrate (msh, fc, 2 * degree + 2, func);

        pass[i] = std::abs (res - ref[i]) < 1.e-12;
        //std::cout << std::setprecision (16) << res << std::endl;
        pass_test = pass[i] && pass_test;
        i++;
    }
    auto res_cell = integrate(msh, cl, 2 * degree + 2, func);
    pass[i] = std::abs (res_cell - ref[fcs.size()]) < 1.e-12;

    std::cout << " * "<< name<< " ...................." << passfail (pass_test) << std::endl;
    std::cout << " * Cell  : "<< res_cell << "  "<<  passfail (pass[fcs.size()]) << std::endl;

    for (auto i = 0; i < fcs.size(); i++)
    {
        std::cout << " * Face[" << i << "] :  " << passfail (pass[i]) << std::endl;
    }
}

template<typename T>
void test_triangles(const size_t& degree)
{
    std::vector<T> ref_tria = { 2. / 3., 0., std::sqrt(2.0) * (8.0 / 15.0),
                2.0/9.0, 0.};
    simplicial_mesh<T, 2>  msh_tria;
    make_single_element_mesh(msh_tria, { 0,0 }, { 1,0 }, { 1,1 });
    test_quadratures("Triangles  ", msh_tria, function_2d, degree, ref_tria);

}

template<typename T>
void test_quadrangles(const size_t& degree)
{
    T x0 (2.0), hx(1.0);
    T y0 (3.0), hy(1.5);

    auto y1 = y0 + hy;
    auto x1 = x0 + hx;
    auto x0_cub = x0 * x0 * x0;
    auto x1_cub = x1 * x1 * x1;
    auto y0_cub = y0 * y0 * y0;
    auto y1_cub = y1 * y1 * y1;

    auto int_vert  = ((y1 - (1. / 3.) * y1_cub) - (y0 - (1. / 3.) * y0_cub));
    auto int_horiz = ((x1 - (1. / 3.) * x1_cub) - (x0 - (1. / 3.) * x0_cub));

    std::vector<T> ref_quads = { (1.0 - y0 * y0) * int_horiz,
        (1.0 - x0 * x0) * int_vert,
        (1.0 - x1 * x1) * int_vert,
        (1.0 - y1 * y1) * int_horiz,
        int_horiz * int_vert };

    cartesian_mesh<T, 2> msh_quads;
    make_single_element_mesh(msh_quads, { x0, y0 }, hx, hy);
    test_quadratures("Quadrangles", msh_quads, function_2d, degree, ref_quads);
}


template<typename T>
void test_tetras(const size_t& degree)
{
    T x0 (0.0), hx(3.0);
    T y0 (0.0), hy(4.0);
    T z0 (0.0), hz(5.0);

    auto x1 = x0 + hx;
    auto y1 = y0 + hy;
    auto z1 = z0 + hz;

    std::vector<T> ref_tetra = { 0., 0., 0., 0., 10. };
    simplicial_mesh<T, 3> msh_tetra;
    make_single_element_mesh(msh_tetra,
        { x1, y0, z0 }, { x0, y1, z0 }, {x0, y0, z1}, { x0, y0, z0});
    test_quadratures("Tetrahedras", msh_tetra, function_unit, degree, ref_tetra);
}


template<typename T>
void test_hexas(const size_t& degree)
{
    T hx(1.0), hy(1.0), hz(1.8);
    point<T, 3> pt(1.5, 2.0, 3.0);
    auto y2 = (pt.y() + hy) * (pt.y() + hy) - pt.y() * pt.y();
    auto x2 = (pt.x() + hx) * (pt.x() + hx) - pt.x() * pt.x();
    auto z2 = (pt.z() + hz) * (pt.z() + hz) - pt.z() * pt.z();

    std::vector<T> ref_hexa = { 0.25 * y2 * z2 * pt.x(),
                                0.25 * y2 * z2 * (pt.x() + hx) ,
                                0.25 * x2 * y2 * pt.z(),
                                0.25 * x2 * y2 * (pt.z() + hz) ,
                                0.25 * x2 * z2 * pt.y(),
                                0.25 * x2 * z2 * (pt.y() + hy),
                                0.125 * x2 * y2 * z2};

    cartesian_mesh<T, 3> msh_hexa;
    make_single_element_mesh(msh_hexa, {pt.x(), pt.y(), pt.z()}, hx, hy, hz);
    test_quadratures("Hexahedras ", msh_hexa, function_3d, degree, ref_hexa);
}


int main(void)
{
    using T = double;
    size_t degree = 6;

    std::cout << "Quadratures Tests: " << std::endl;

    test_triangles<T>(degree);

    test_quadrangles<T>(degree);

    test_tetras<T>(degree);

    test_hexas<T>(degree);
}
