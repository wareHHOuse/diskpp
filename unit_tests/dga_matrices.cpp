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

    disk::vec3<T> field = { 1.0, 1.0, 1.0 };

    T ee = 0.0;
    T me = 0.0;
    for (auto& cl : msh)
    {
        /* Make the edge matrix, project the test field
         * on the edges of the current element to obtain
         * voltages across edges, then compute electric
         * energy. Accumulate in ee. */
        auto Me = edge_matrix(msh, cl, 1.0);
        Eigen::Matrix<T, 6, 1> Us;
        auto pev = primal_edge_vectors(msh, cl);
        for (size_t i = 0; i < 6; i++)
            Us(i) = field.dot(pev[i]);
        ee += Us.dot(Me*Us);

        /* Make the face matrix, project the test field
         * on the faces of the current element to obtain
         * fluxes across faces, then compute magnetic
         * energy. Accumulate in me. */
        auto Mf = face_matrix(msh, cl, 1.0);
        Eigen::Matrix<T, 4, 1> Phis;
        auto pav = primal_area_vectors(msh, cl);
        for (size_t i = 0; i < 4; i++)
            Phis(i) = field.dot(pav[i]);
        me += Phis.dot(Mf*Phis);
    }

    T ee_err = ee - field.dot(field);
    T me_err = me - field.dot(field); 

    std::cout << "Electric energy: " << ee << ", error: " << ee_err << std::endl;
    std::cout << "Magnetic energy: " << me << ", error: " << me_err << std::endl;

    if ( ee_err > 1e-16 or me_err > 1e-16)
        return 1;

    return 0;
}