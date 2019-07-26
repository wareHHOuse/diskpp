#include <iostream>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "geometry/geometry.hpp"
#include "loaders/loader.hpp"
#include "methods/hho"
#include "solvers/solver.hpp"

#include "../tests/common.hpp"

using namespace disk;


template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
oldstab(const Mesh&                                                     msh,
        const typename Mesh::cell_type&                                 cl,
        const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction,
        const hho_degree_info&                                          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    // Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)

    // Step 1: compute \pi_T^k p_T^k v (third term).
    const matrix_type M1    = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type M2    = mass_mat.block(0, 1, cbs, rbs - 1);
    matrix_type       proj1 = -M1.ldlt().solve(M2 * reconstruction);

    assert(M2.cols() == reconstruction.rows());

    // Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
    proj1.block(0, 0, cbs, cbs) += matrix_type::Identity(cbs, cbs);

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(cbs + num_faces * fbs, cbs + num_faces * fbs);

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        const auto fc = fcs[face_i];
        const auto hf = diameter(msh, fc);
        auto       fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix.block(0, 1, fbs, rbs - 1);

        matrix_type       proj2 = piKF.solve(MR1 * reconstruction);
        proj2.block(0, cbs + face_i * fbs, fbs, fbs) -= matrix_type::Identity(fbs, fbs);

        // Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
        const matrix_type MR2   = face_trace_matrix.block(0, 0, fbs, cbs);
        const matrix_type proj3 = piKF.solve(MR2 * proj1);
        const matrix_type BRF   = proj2 + proj3;

        data += BRF.transpose() * face_mass_matrix * BRF / hf;
    }

    return data;
}

template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
newstab(const Mesh&                                                     msh,
        const typename Mesh::cell_type&                                 cl,
        const Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>& reconstruction,
        const hho_degree_info&                                          di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto recdeg = di.reconstruction_degree();
    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();

    const auto rbs = scalar_basis_size(recdeg, Mesh::dimension);
    const auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

    const auto cb = make_scalar_monomial_basis(msh, cl, recdeg);

    const matrix_type mass_mat = make_mass_matrix(msh, cl, cb);

    const matrix_type projT = mass_mat.block(0, 0, cbs, cbs);
    const matrix_type Reval = mass_mat.block(0, 1, cbs, rbs - 1);
    const matrix_type Re = Reval * reconstruction;
    matrix_type       pR = projT.ldlt().solve(Re);

    matrix_type rdiff = matrix_type::Zero(rbs, reconstruction.cols());
    rdiff.block(1,0,rbs-1,reconstruction.cols()) = -reconstruction;
    assert(rdiff.cols() == pR.cols());
    rdiff.block(0,0,cbs,pR.cols()) += pR.block(0,0,cbs,pR.cols());
    //rdiff(0,0) = 1;

    const auto fcs       = faces(msh, cl);
    const auto num_faces = fcs.size();

    matrix_type data = matrix_type::Zero(/*cbs + num_faces * */fbs, cbs + num_faces * fbs);
/*
    std::cout << "pR" << std::endl;
    std::cout << pR << std::endl;
    std::cout << "rdiff" << std::endl;
    std::cout << rdiff << std::endl;
*/

    // Step 3: project on faces (eqn. 21)
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        const auto fc = fcs[face_i];
        const auto hf = 1.;//diameter(msh, fc);
        auto       fb = make_scalar_monomial_basis(msh, fc, facdeg);

        matrix_type face_mass_matrix  = make_mass_matrix(msh, fc, fb);
        matrix_type face_trace_matrix = matrix_type::Zero(fbs, rbs);

        const auto face_quadpoints = integrate(msh, fc, recdeg + facdeg);
        for (auto& qp : face_quadpoints)
        {
            assert(qp.weight() > 0);
            const auto        f_phi   = fb.eval_functions(qp.point());
            const auto        c_phi   = cb.eval_functions(qp.point());
            face_trace_matrix += priv::outer_product(priv::inner_product(qp.weight(), f_phi), c_phi);
        }

        /* Build v_F - v_T term */

        LLT<matrix_type> piKF;
        piKF.compute(face_mass_matrix);

        matrix_type T = face_trace_matrix.block(0, 0, fbs, cbs);

        // Step 3a: \pi_F^k( v_F - p_T^k v )
        const matrix_type MR1 = face_trace_matrix.block(0, 0, fbs, rbs);

        matrix_type RRR = MR1*rdiff;

        matrix_type       BRF = piKF.solve(RRR);


        //data -= BRF.transpose() * face_mass_matrix * BRF / hf;

        
        matrix_type MT = piKF.solve(T);
        /*
        data.block(0,0,cbs,cbs) += T.transpose() * MT / hf;
        data.block(cbs + face_i * fbs, cbs + face_i * fbs, fbs, fbs) += face_mass_matrix/hf;
        data.block(0, cbs + face_i * fbs, cbs, fbs) += -T.transpose()/hf;
        data.block(cbs + face_i * fbs, 0, fbs, cbs) += -T/hf;
        
        data -= RRR.transpose() * BRF / hf;
        */

        matrix_type D = matrix_type::Zero(fbs, cbs + num_faces * fbs);
        D.block(0,0,fbs,cbs) = -MT;
        D.block(0,cbs+face_i*fbs,fbs,fbs) = matrix_type::Identity(fbs, fbs);

        data += D;
    }

    return data;
}

template<typename Mesh>
void
newstab_test(const Mesh& msh, size_t degree)
{
    using T = typename Mesh::coordinate_type;
    hho_degree_info hdi(degree);

    auto f = make_scalar_testing_data(msh);

    for (auto& cl : msh)
    {
        std::cout << cl << std::endl;
        auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = oldstab(msh, cl, gr.first, hdi);
        auto stab2   = newstab(msh, cl, gr.first, hdi);
    
        auto proj = disk::project_function(msh, cl, hdi, f, 2);
        
        Matrix<T, Dynamic, Dynamic> s = stab2 * proj;

        std::cout << proj.transpose() << std::endl;

        std::cout << s << std::endl;

        /*
        std::cout << " ** CELL: **" << cl << std::endl;
        std::cout << "OLD STAB" << std::endl;
        std::cout << stab << std::endl;
        std::cout << "NEW STAB" << std::endl;
        std::cout << stab2 << std::endl;
        */
        break;
    }
}

int main(int argc, char **argv)
{
    using T = double;

    if (argc != 3)
    {
        std::cout << "Please specify file name." << std::endl;
        return 1;
    }

    char *mesh_filename = argv[1];

    /* DiSk++ cartesian 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: DiSk++ Cartesian 2D" << std::endl;
        auto msh = disk::load_cartesian_2d_mesh<T>(mesh_filename);
        newstab_test(msh, atoi(argv[2]));
        return 0;
    }
}
