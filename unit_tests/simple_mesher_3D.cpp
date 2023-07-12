/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Karol Cascavita (C) 2023
 *
 */

/* Test primitive operations like barycenter(), measure() etc. */

#include <iostream>
#include <fstream>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/geometry/geometry.hpp"

#define MEAS_THRESH 1e-15
#define BARY_THRESH 1e-15

#include "diskpp/output/silo.hpp"

template<typename Mesh>
int export_mesh_to_silo(Mesh& msh, const char* silo_filename, const std::vector<double>& en)
{
    disk::silo_database silo_db;
    silo_db.create(silo_filename);
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<double> ens("en", en);
    silo_db.add_variable("mesh", ens);

    silo_db.close();

    return 0;
}


template<typename Mesh>
bool test(Mesh& msh,  size_t num_ref)
{
    size_t num_tetras = 6 * std::pow(8, num_ref);
    size_t num_fcs_per_boundary = 2 * std::pow(4, num_ref);

    if (msh.cells_size() != num_tetras)
        return false;

    // Volume verification
    for (auto& cl : msh)
        if (measure(msh, cl) != (1.0 /num_tetras))
            return false;

    // Number of boundary faces
    std::array<size_t, 6> counter{};
    for (auto& fc : boundary_faces(msh))
    {
        auto bi = msh.boundary_info(fc);
        counter.at(bi.tag()) += 1;
    }

    for (auto& c : counter)
        if (c != num_fcs_per_boundary)
            return false;

    return true;
}


int main(void)
{
    using T = double;

    size_t num_ref = 2;


    std::cout << "Simplicial 3D" << std::endl;
    disk::simplicial_mesh<T, 3> msh;
    auto mesher = make_simple_mesher(msh);
    for(size_t r = 0; r < num_ref; r++)
        mesher.refine();

    /*
    std::vector<double> elemnum;
    size_t n = 0;
    for (auto& cl : msh)
        elemnum.push_back(n++);
    export_mesh_to_silo(msh, "tetras.silo", elemnum);
    */
    return test(msh, num_ref);
}