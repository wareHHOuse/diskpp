/*
 *       /\         DISK++, a template library for DIscontinuous SKeletal
 *      /__\        methods.
 *     /_\/_\
 *    /\    /\      Matteo Cicuttin (C) 2016, 2017, 2018
 *   /__\  /__\     matteo.cicuttin@enpc.fr
 *  /_\/_\/_\/_\    École Nationale des Ponts et Chaussées - CERMICS
 *
 * This file is copyright of the following authors:
 * Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
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

#include "common.hpp"

template<typename T>
bool test_barycenter(const point<T,1>& pt)
{
    T ULP_max = 2;
    
    if ( almost_equal(pt.x(), 0.5, ULP_max) )
        return true;
    
    std::cout << "  pt.x() not accurate: " << std::setprecision(16) << pt.x() << std::endl;
    return false;
}

template<typename T>
bool test_barycenter(const point<T,2>& pt)
{
    T ULP_max = 2;
    bool success = true;
    
    if ( !almost_equal(pt.x(), 0.5, ULP_max) )
    {
        std::cout << "  pt.x() not accurate: ";
        std::cout << std::setprecision(16) << pt.x() << std::endl;
        success = false;
    }
    
    if ( !almost_equal(pt.y(), 0.5, ULP_max) )
    {
        std::cout << "  pt.y() not accurate: ";
        std::cout << std::setprecision(16) << pt.y() << std::endl;
        success = false;
    }
    
    return success;
}

template<typename T>
bool test_barycenter(const point<T,3>& pt)
{
    T ULP_max = 2;
    bool success = true;
    
    if ( !almost_equal(pt.x(), 0.5, ULP_max) )
    {
        std::cout << "  pt.x() not accurate: ";
        std::cout << std::setprecision(16) << pt.x() << std::endl;
        success = false;
    }
    
    if ( !almost_equal(pt.y(), 0.5, ULP_max) )
    {
        std::cout << "  pt.y() not accurate: ";
        std::cout << std::setprecision(16) << pt.y() << std::endl;
        success = false;
    }
    
    if ( !almost_equal(pt.z(), 0.5, ULP_max) )
    {
        std::cout << "  pt.z() not accurate: ";
        std::cout << std::setprecision(16) << pt.z() << std::endl;
        success = false;
    }
    
    return success;
}

template<typename Mesh>
bool
test_barycenter(const Mesh& msh)
{
    typedef typename Mesh::point_type point_type;
    typedef typename Mesh::coordinate_type T;
    point_type  tot_bar;
    T           tot_meas;
    
    for (auto& cl : msh)
    {
        auto meas = measure(msh, cl);
        auto bar = barycenter(msh,cl);
        
        tot_bar = tot_bar + meas*bar;
        tot_meas += meas;
    }
    
    test_barycenter(tot_bar);
    
    return true;
}

template<typename T>
void test(const T& meshes)
{
    size_t num_meshes = meshes.size();
    size_t cur_mesh = 1;
    for (auto& msh : meshes)
    {
        std::cout << cyan << " Mesh " << cur_mesh << "/";
        std::cout << num_meshes << nocolor << std::endl;
        test_barycenter(msh);
        cur_mesh++;
    }
}

int main(void)
{
    std::cout << bold << yellow << "Testing triangular meshes (generic)" << reset << std::endl;
    test( get_triangle_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing triangular meshes (netgen)" << reset << std::endl;
    test( get_triangle_netgen_meshes<double>() );
    std::cout << bold << yellow << "Testing quad meshes" << reset << std::endl;
    test( get_quad_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing polygonal meshes" << reset << std::endl;
    test( get_polygonal_generic_meshes<double>() );
    std::cout << bold << yellow << "Testing tetrahedral meshes" << reset << std::endl;
    test( get_tetrahedra_netgen_meshes<double>() );
}
