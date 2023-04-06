/*
 *       /\        Matteo Cicuttin (C) 2016-2020
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

/* This DiSk++ tutorial shows how to use GMSH meshes.
 */



#include <iostream>

// For the mesh data structure
#include "diskpp/mesh/mesh.hpp"

// For the loaders and related helper functions
#include "diskpp/loaders/loader.hpp"

#include "diskpp/output/silo.hpp"

using disk::cells;
using disk::faces;

template<typename Mesh>
void run(Mesh& msh)
{
    /* Iterate on the mesh cells */
    for (auto& cl : cells(msh))
    {
        auto di = msh.domain_info(cl);
        std::cout << "Cell " << cl << ": domain id is " << di.tag() << std::endl;
    }

    /* iterate on faces */
    for (auto& fc : faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        std::cout << "Face " << fc << ": ";
        std::cout << "is boundary = " << bi.is_boundary();
        if (bi.is_boundary())
        {
            std::cout << ", bnd id = " << bi.tag();
            
            if (bi.is_internal())
                std::cout << ", internal boundary";
        }
        std::cout << std::endl;
    }

    disk::silo_database silo;
    silo.create("mesh.silo");
    silo.add_mesh(msh, "mesh");
    silo.close();
}

int main(int argc, const char **argv)
{
    if (argc < 2)
    {
        std::cout << "Please specify a mesh" << std::endl;
        return 1;
    }

    const char *mesh_filename = argv[1];

    using T = double;

    /* GMSH default extension for models is .geo and for meshes is .msh.
     * It is therefore not possible to distinguish programmatically if a
     * mesh is 2D/3D simplicial or polyhedral. You need then to change the
     * extensions as follows. */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial 2D" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        
        run(msh);

        return 0;
    }
    
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH simplicial 3D" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        run(msh);
        
        return 0;
    }
    
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH poly 3D" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);
        
        run(msh);

        return 0;
    }


    return 0;
}
