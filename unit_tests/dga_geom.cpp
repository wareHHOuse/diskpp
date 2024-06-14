/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2024
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include "diskpp/methods/dga"
#include "diskpp/mesh/meshgen.hpp"

int main(int argc, char **argv)
{
    using T = double;
    using mesh_type = disk::simplicial_mesh<T, 3>;

    /* Test domain is [0,1]^3 */
    mesh_type msh;
    auto mesher = make_simple_mesher(msh);
    for (int i = 0; i < 2; i++)
        mesher.refine();

    for (auto& cl : msh)
    {
        Eigen::Matrix<T,3,3> vol1 = Eigen::Matrix<T,3,3>::Zero();
        auto pev = primal_edge_vectors(msh, cl);
        auto dav = dual_area_vectors(msh, cl);
        for (size_t i = 0; i < 6; i++)
            vol1 += pev[i]*dav[i].transpose();

        Eigen::Matrix<T,3,3> vol2 = Eigen::Matrix<T,3,3>::Zero();
        auto pav = primal_area_vectors(msh, cl);
        auto dev = dual_edge_vectors(msh, cl);
        for (size_t i = 0; i < 4; i++)
            vol2 += pav[i]*dev[i].transpose();

        T expected = volume_unsigned(msh, cl);
        T computed = vol1.array().sum() + vol2.array().sum();
        T error = std::abs(computed - 6*expected);

        if (error > 1e-16) {
            std::cout << "Error too big in " << cl << ": ";
            std::cout << error << std::endl;
            return 1;
        }
    }

    return 0;
}