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

#include "core/loaders/loader.hpp"
#include "geometry/geometry.hpp"
#include "common.hpp"

#include "mesh/rational.hpp"

point<rational<int>,2>
barycenter(const std::vector<point<rational<int>,2>>& mp,
           size_t a, size_t b, size_t c)
{
    return (mp[a] + mp[b] + mp[c])/3;
}

rational<int>
area(const std::vector<point<rational<int>,2>>& mp,
     size_t a, size_t b, size_t c)
{
    return abs(det(mp[b] - mp[a], mp[c] - mp[a])/2);
}

template<typename Mesh>
void x(const Mesh& msh)
{
    using R = rational<int>;
    using P = point<R,2>;
    
    std::vector<P> mp;
    mp.resize(14);
    mp[0]  = P( R(0,10), R(0,10) );
    mp[1]  = P( R(1,10), R(0,10) );
    mp[2]  = P( R(3,10), R(0,10) );
    mp[3]  = P( R(6,10), R(0,10) );
    mp[4]  = P( R(6,10), R(2,10) );
    mp[5]  = P( R(4,10), R(3,10) );
    mp[6]  = P( R(2,10), R(5,10) );
    mp[7]  = P( R(6,10), R(5,10) );
    mp[8]  = P( R(4,10), R(8,10) );
    mp[9]  = P( R(0,10), R(10,10) );
    mp[10] = P( R(1,10), R(10,10) );
    mp[11] = P( R(2,10), R(10,10) );
    mp[12] = P( R(4,10), R(10,10) );
    mp[13] = P( R(6,10), R(10,10) );
    
    for (auto& cl : msh)
    {
        auto ptids = cl.point_ids();
        
        auto rbar = barycenter(mp, ptids[0], ptids[1], ptids[2]);
        auto rmeas = area(mp, ptids[0], ptids[1], ptids[2]);
        
        auto bar = barycenter(msh, cl);
        auto meas = measure(msh, cl);
        
        double ULP_max = 1;
        bool success = true;
        
        if ( !almost_equal(bar.x(), double(rbar.x()), ULP_max) )
        {
            std::cout << "  pt.x() not accurate: ";
            std::cout << std::setprecision(16) << bar.x() << " ";
            std::cout << rbar.x() << std::endl;
            success = false;
        }
        
        if ( !almost_equal(bar.y(), double(rbar.y()), ULP_max) )
        {
            std::cout << "  pt.y() not accurate: ";
            std::cout << std::setprecision(16) << bar.y() << " ";
            std::cout << rbar.y() << std::endl;
            success = false;
        }
    }
}

int main(void)
{
    using T = double;
    
    typedef disk::simplicial_mesh<T, 2>  mesh_type;

    mesh_type msh;
    disk::netgen_mesh_loader<T, 2> loader;

    if (!loader.read_mesh("../../../diskpp/meshes/tests/pattern1_scale1.mesh2d"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }

    loader.populate_mesh(msh);

    x(msh);
    return 0;
}
