/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <iostream>

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "common/eigen.hpp"



int main(int argc, char **argv)
{
    using RealType = double;

    //const char *msh_fn = "/Users/matteo/Desktop/MsHHO meshes/mesh6.mesh2d";
    const char *msh_fn = "../multiscale/monoscale_hho_mesh.mesh2d";
    const char *sol_fn = "../diffusion/solution_monoscale_k0.bin";

    typedef disk::simplicial_mesh<RealType, 2>      mesh_type;
    typedef typename mesh_type::cell                cell_type;
    //typedef typename mesh_type::face                face_type;
    typedef typename mesh_type::scalar_type         scalar_type;
    
    mesh_type msh;
    
    disk::netgen_mesh_loader<RealType, 2> loader;
    loader.verbose(true);
    if (!loader.read_mesh(msh_fn))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);
    
    

    std::ifstream ifs(sol_fn, std::ifstream::binary);
    
    size_t sol_num_elements, cell_basis_deg, face_basis_deg;
    
    ifs.read(reinterpret_cast<char *>(&sol_num_elements), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&cell_basis_deg), sizeof(size_t));
    ifs.read(reinterpret_cast<char *>(&face_basis_deg), sizeof(size_t));
    
    if (sol_num_elements != msh.cells_size())
    {
        std::cout << "Solution has a different number of elements than the mesh (";
        std::cout << sol_num_elements << " vs. " << msh.cells_size() << ")" << std::endl;
        return 1;
    }
    
    //std::cout << "Solution has " << sol_num_elements << " elements" << std::endl;
    //std::cout << "Cell basis degree: " << cell_basis_deg << std::endl;
    //std::cout << "Face basis degree: " << face_basis_deg << std::endl;
    
    
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    //typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;
    
    cell_basis_type     cell_basis_k1(cell_basis_deg+1);
    
    size_t local_dofs_size = cell_basis_k1.size();
    size_t dofs_vec_size = msh.cells_size() * local_dofs_size;
    dynamic_vector<scalar_type> cell_dofs = dynamic_vector<scalar_type>::Zero(dofs_vec_size);
    
    for (size_t i = 0; i < dofs_vec_size; i++)
        ifs.read(reinterpret_cast<char *>(&cell_dofs(i)), sizeof(scalar_type));
    
    ifs.close();
    
    std::ofstream ofs("pippo.dat");
    
    size_t elemnum = 0;
    for (auto& cl : msh)
    {
        auto offset_start = elemnum * local_dofs_size;
        dynamic_vector<scalar_type> local_dofs;
        local_dofs = cell_dofs.block(offset_start, 0, local_dofs_size, 1);
        
        auto bar = barycenter(msh, cl);

        auto phi = cell_basis_k1.eval_functions(msh, cl, bar);
        auto val = local_dofs.dot(phi);
        
        elemnum++;
        
        ofs << bar.x() << " " << bar.y() << " " << val << std::endl;
    }
    
    ofs.close();
    
    return 0;

}




