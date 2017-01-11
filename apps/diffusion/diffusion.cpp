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
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <map>

#include "colormanip.h"

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
    #include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

struct assembly_info
{
    size_t  linear_system_size;
    double  time_gradrec, time_statcond, time_stab;
};

struct solver_info
{
    double  time_solver;
};

struct postprocess_info
{
    double  time_postprocess;
};

template<typename Mesh>
class diffusion_solver
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::diffusion_like_static_condensation_bq<bqdata_type>    statcond_type;
    typedef disk::assembler<mesh_type, face_basis_type, face_quadrature_type> assembler_type;

    size_t m_cell_degree, m_face_degree;

    bqdata_type     m_bqd;

    typename assembler_type::sparse_matrix_type     m_system_matrix;
    typename assembler_type::vector_type            m_system_rhs, m_system_solution;

    const mesh_type& m_msh;

    std::vector<dynamic_vector<scalar_type>>        m_postprocess_data;

    bool m_verbose;

public:
    diffusion_solver(const mesh_type& msh, size_t degree, int l = 0)
        : m_msh(msh), m_verbose(false)
    {
        if ( l < -1 or l > 1)
        {
            if (degree > 0)
            {
                std::cout << "'l' should be -1, 0 or 1. Reverting to 0." << std::endl;
                l = 0;
            }
            else
            {
                std::cout << "'l' should be 0 or 1. Reverting to 0." << std::endl;
                l = 0;
            }
        }

        m_cell_degree = degree + l;
        m_face_degree = degree;

        m_bqd = bqdata_type(m_cell_degree, m_face_degree);
    }

    bool    verbose(void) const     { return m_verbose; }
    void    verbose(bool v)         { m_verbose = v; }

    template<typename LoadFunction, typename BoundaryConditionFunction>
    assembly_info
    assemble(const LoadFunction& lf, const BoundaryConditionFunction& bcf)
    {
        auto gradrec    = gradrec_type(m_bqd);
        auto stab       = stab_type(m_bqd);
        auto statcond   = statcond_type(m_bqd);
        auto assembler  = assembler_type(m_msh, m_face_degree);

        assembly_info ai;
        bzero(&ai, sizeof(ai));

        timecounter tc;

        for (auto& cl : m_msh)
        {
            tc.tic();
            gradrec.compute(m_msh, cl);
            tc.toc();
            ai.time_gradrec += tc.to_double();

            tc.tic();
            stab.compute(m_msh, cl, gradrec.oper);
            tc.toc();
            ai.time_stab += tc.to_double();

            tc.tic();
            auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            auto scnp = statcond.compute(m_msh, cl, loc, cell_rhs);
            tc.toc();
            ai.time_statcond += tc.to_double();

            assembler.assemble(m_msh, cl, scnp);
        }

        assembler.impose_boundary_conditions(m_msh, bcf);
        assembler.finalize(m_system_matrix, m_system_rhs);

        ai.linear_system_size = m_system_matrix.nonZeros();
        return ai;
    }

    solver_info
    solve(void)
    {
#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

        solver_info si;

        size_t systsz = m_system_matrix.rows();
        size_t nnz = m_system_matrix.nonZeros();

        if (verbose())
        {
            std::cout << "Starting linear solver..." << std::endl;
            std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
            std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;
        }

        timecounter tc;

        tc.tic();
        solver.analyzePattern(m_system_matrix);
        solver.factorize(m_system_matrix);
        m_system_solution = solver.solve(m_system_rhs);
        tc.toc();
        si.time_solver = tc.to_double();

        return si;
    }

    template<typename LoadFunction>
    postprocess_info
    postprocess(const LoadFunction& lf)
    {
        auto gradrec    = gradrec_type(m_bqd);
        auto stab       = stab_type(m_bqd);
        auto statcond   = statcond_type(m_bqd);

        size_t fbs = m_bqd.face_basis.size();

        postprocess_info pi;

        m_postprocess_data.reserve(m_msh.cells_size());

        timecounter tc;
        tc.tic();
        for (auto& cl : m_msh)
        {
            auto fcs = faces(m_msh, cl);
            auto num_faces = fcs.size();

            dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto fc = fcs[face_i];
                auto eid = find_element_id(m_msh.faces_begin(), m_msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                auto face_id = eid.second;

                dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
                xF = m_system_solution.block(face_id * fbs, 0, fbs, 1);
                xFs.block(face_i * fbs, 0, fbs, 1) = xF;
            }

            gradrec.compute(m_msh, cl);
            stab.compute(m_msh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(m_msh, cl, lf, m_cell_degree);
            dynamic_vector<scalar_type> x = statcond.recover(m_msh, cl, loc, cell_rhs, xFs);
            m_postprocess_data.push_back(x);
        }
        tc.toc();

        pi.time_postprocess = tc.to_double();

        return pi;
    }

    template<typename AnalyticalSolution>
    scalar_type
    compute_l2_error(const AnalyticalSolution& as)
    {
        scalar_type err_dof = 0.0;

        disk::projector_bq<bqdata_type> projk(m_bqd);

        size_t i = 0;
        for (auto& cl : m_msh)
        {
            auto x = m_postprocess_data.at(i++);
            dynamic_vector<scalar_type> true_dof = projk.compute_cell(m_msh, cl, as);
            dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
            dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
            err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
        }

        return sqrt(err_dof);
    }

    void
    plot_solution(const std::string& filename)
    {
        std::ofstream ofs(filename);

        size_t cell_i = 0;
        for (auto& cl : m_msh)
        {
            auto x = m_postprocess_data.at(cell_i++);
            auto qps = m_bqd.cell_quadrature.integrate(m_msh, cl);
            for (auto& qp : qps)
            {
                auto phi = m_bqd.cell_basis.eval_functions(m_msh, cl, qp.point());

                scalar_type pot = 0.0;
                for (size_t i = 0; i < m_bqd.cell_basis.range(0, m_cell_degree).size(); i++)
                    pot += phi[i] * x(i);

                auto tp = qp.point();
                for (size_t i = 0; i < mesh_type::dimension; i++)
                    ofs << tp[i] << " ";
                ofs << pot << std::endl;
            }
        }

        ofs.close();
    }
};

#if 0

template<typename MeshType, typename LoaderType>
void
test_mesh_format(const std::vector<std::string>& paths,
                 size_t runs, size_t mindeg, size_t maxdeg,
                 const std::string& output_basename)
{

    auto f = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return M_PI * M_PI * sin(p.x() * M_PI);
    };

    auto sf = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return sin(p.x() * M_PI);
    };

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        std::stringstream ss;
        ss << output_basename << "_degree_" << i << ".txt";

        std::ofstream ofs(ss.str());

        size_t mesh_index = 1;
        for (auto& tsp : paths)
        {
            MeshType    msh;
            LoaderType  loader;

            if (!loader.read_mesh(tsp))
            {
                std::cout << "Problem loading mesh." << std::endl;
                return;
            }
            loader.populate_mesh(msh);


            struct timing_data td;
            td.time_gradrec     = 0.0;
            td.time_stab        = 0.0;
            td.time_statcond    = 0.0;
            td.time_solver      = 0.0;

            for (size_t run = 0; run < runs; run++)
            {
                struct timing_data rtd;
                rtd = test_diffusion(msh, f, sf, i, "plot.dat");

                td.time_gradrec     += rtd.time_gradrec;
                td.time_stab        += rtd.time_stab;
                td.time_statcond    += rtd.time_statcond;
                td.time_solver      += rtd.time_solver;
                td.solved_dofs      = rtd.solved_dofs;
            }

            td.time_gradrec /= double(runs);
            td.time_stab /= double(runs);
            td.time_statcond /= double(runs);
            td.time_solver /= double(runs);

            ofs << td.solved_dofs << " ";
            ofs << td.time_gradrec << " ";
            ofs << td.time_stab << " ";
            ofs << td.time_statcond << " ";
            ofs << td.time_solver << std::endl;

            std::cout << "Degree:                  " << i << std::endl;
            std::cout << "DOFs:                    " << td.solved_dofs << std::endl;
            std::cout << "Gradient reconstruction: " << td.time_gradrec << std::endl;
            std::cout << "Stabilization:           " << td.time_stab << std::endl;
            std::cout << "Static condensation:     " << td.time_statcond << std::endl;
            std::cout << "Solver:                  " << td.time_solver << std::endl;

        }

        ofs.close();
    }
}

#endif


template<typename MeshType, typename LoaderType>
bool
verify_convergence(const std::vector<std::string>& paths,
                   size_t mindeg, size_t maxdeg)
{
    typedef typename MeshType::scalar_type scalar_type;

    auto f = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return M_PI * M_PI * sin(p.x() * M_PI);
    };

    auto sf = [](const point<typename MeshType::scalar_type, MeshType::dimension>& p) -> auto {
        return sin(p.x() * M_PI);
    };

    bool success = true;

    for (size_t i = mindeg; i <= maxdeg; i++)
    {
        scalar_type expected_rate = i+2;

        std::vector<std::pair<scalar_type, scalar_type>> errdiams;

        for (auto& tsp : paths)
        {
            MeshType    msh;
            LoaderType  loader;

            if (!loader.read_mesh(tsp))
            {
                std::cout << "Problem loading mesh." << std::endl;
                return false;
            }
            loader.populate_mesh(msh);

            diffusion_solver<MeshType> dp(msh, i);

            dp.assemble(f, sf);
            dp.solve();
            dp.postprocess(f);
            auto error = dp.compute_l2_error(sf);
            auto diam = average_diameter(msh);

            errdiams.push_back( std::make_pair(diam, error) );
        }

        std::cout << "Convergence rates for k = " << i << ":   ";

        bool pass       = true;
        bool warning    = false;
        bool high, low, ok;
        for (size_t i = 1; i < errdiams.size(); i++)
        {
            auto d = log2(errdiams[i-1].first/errdiams[i].first);
            auto e = log2(errdiams[i-1].second/errdiams[i].second);
            auto rate = e/d;

            ok   = (std::abs(expected_rate - rate) < 0.4); /* Test passed */
            low  = ((expected_rate - rate) > 0.2); /* a bit too low, warn */
            high = ((rate - expected_rate) > 0.2); /* a bit too high, warn */

            if (low)    std::cout << magenta;
            if (high)   std::cout << cyan;
            std::cout << std::fixed << std::setprecision(3) << rate << "  ";
            if (low or high)
            {
                std::cout << reset;
                warning = true;
            }
        }

        std::string             okfail = "[\x1b[31;1mFAIL\x1b[0m]";
        if (ok && not warning)  okfail = "[\x1b[32;1m OK \x1b[0m]";
        if (ok && warning)      okfail = "[\x1b[33;1mWARN\x1b[0m]";

        std::cout << okfail << std::endl;

        success &= ok;
    }

    return success;
}


enum test_type
{
    TEST_VERIFY_CONVERGENGE,
    TEST_MEASURE_TIMES
};

void test_triangles_specialized()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
    paths.push_back("../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");

    typedef disk::simplicial_mesh<double, 2>      MeshType;
    typedef disk::netgen_mesh_loader<double, 2>   LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "triangle_spec");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_triangles_generic()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_1.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_triangles/fvca5/mesh1_4.typ1");

    typedef disk::generic_mesh<double, 2>       MeshType;
    typedef disk::fvca5_mesh_loader<double, 2>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "triangle_gen");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_hexagons_generic()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    paths.push_back("../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");

    typedef disk::generic_mesh<double, 2>       MeshType;
    typedef disk::fvca5_mesh_loader<double, 2>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexagons_gen");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}


void test_hexahedra_specialized()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

    typedef disk::cartesian_mesh<double, 3>         MeshType;
    typedef disk::cartesian_mesh_loader<double, 3>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexahedra_spec");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_hexahedra_generic()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

    typedef disk::generic_mesh<double, 3>       MeshType;
    typedef disk::fvca6_mesh_loader<double, 3>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexahedra_gen");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_tetrahedra_specialized()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet1.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet2.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet3.mesh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/netgen/fvca6_tet4.mesh");

    typedef disk::simplicial_mesh<double, 3>         MeshType;
    typedef disk::netgen_mesh_loader<double, 3>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "tetrahedra_spec");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_tetrahedra_generic()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.1.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.2.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.3.msh");
    paths.push_back("../../../diskpp/meshes/3D_tetras/fvca6/tet.4.msh");

    typedef disk::generic_mesh<double, 3>       MeshType;
    typedef disk::fvca6_mesh_loader<double, 3>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "tetrahedra_gen");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}

void test_polyhedra_generic()
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
    //paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
    //paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

    typedef disk::generic_mesh<double, 3>       MeshType;
    typedef disk::fvca6_mesh_loader<double, 3>  LoaderType;

    //test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "polyhedra_gen");

    verify_convergence<MeshType, LoaderType>(paths, 0, 3);
}


#if 0
int main(int argc, char **argv)
{
    test_triangles_specialized();
    #if 0
    using RealType = double;
    typedef disk::cartesian_mesh<RealType, 3>   mesh_type;
    mesh_type msh;

    disk::cartesian_mesh_loader<RealType, 3> loader;
    if (!loader.read_mesh("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex"))
    {
        std::cout << "Problem loading mesh." << std::endl;
        return 1;
    }
    loader.populate_mesh(msh);

    auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
        return M_PI * M_PI * sin(p.x() * M_PI);
        //return 1.0;
    };

    auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
        return sin(p.x() * M_PI);
        //return -p.x() * p.x() * 0.5;
    };

    diffusion_solver<mesh_type> solver(msh, 0);

    solver.assemble(f, sf);
    solver.solve();
    solver.postprocess(f);
    std::cout << solver.compute_l2_error(sf) << std::endl;
    solver.plot_solution("test.dat");
    #endif
}

#endif

int main(int argc, char **argv)
{
    std::cout << bold << underline << "Triangles specialized" << reset << std::endl;
    test_triangles_specialized();

    std::cout << bold << underline << "Triangles generic" << reset << std::endl;
    test_triangles_generic();

    std::cout << bold << underline << "Hexagons" << reset << std::endl;
    test_hexagons_generic();

    std::cout << bold << underline << "Hexahedra specialized" << reset << std::endl;
    test_hexahedra_specialized();

    std::cout << bold << underline << "Hexahedra generic" << reset << std::endl;
    test_hexahedra_generic();

    std::cout << bold << underline << "Tetrahedra specialized" << reset << std::endl;
    test_tetrahedra_specialized();

    std::cout << bold << underline << "Tetrahedra generic" << reset << std::endl;
    test_tetrahedra_generic();

    std::cout << bold << underline << "Polyhedra" << reset << std::endl;
    test_polyhedra_generic();
}


#if 0



int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int     degree          = 1;
    int     l               = 0;
    int     elems_1d        = 8;
    bool    submesh_flag    = false;
    int ch;

    while ( (ch = getopt(argc, argv, "k:n:sl:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree < 0)
                {
                    std::cout << "Degree must be positive. Falling back to 1." << std::endl;
                    degree = 1;
                }
                break;

            case 'n':
                elems_1d = atoi(optarg);
                if (elems_1d < 0)
                {
                    std::cout << "Num of elems must be positive. Falling back to 8." << std::endl;
                    elems_1d = 8;
                }
                break;

            case 's':
                submesh_flag = true;
                break;

            case 'l':
                l = atoi(optarg);
                if (l < -1 or l > 1)
                {
                    std::cout << "l can be -1, 0 or 1. Falling back to 0." << std::endl;
                    l = 0;
                }
                break;

            case 'h':
            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;


    if (argc == 0)
    {
        std::cout << "Running 1D test simulation" << std::endl;

        typedef disk::generic_mesh<RealType, 1>  mesh_type;

        mesh_type msh;
        disk::uniform_mesh_loader<RealType, 1> loader(0,1,elems_1d);
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            //return 2.;
            return M_PI * M_PI * sin(p.x() * M_PI);
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            //return - p.x() * p.x();
            return sin(p.x() * M_PI);
        };

        //test_gradrec(msh, degree);
        test_diffusion(msh, f, sf, degree, "plot.dat");

        return 0;
    }

    filename = argv[0];

    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 6. * p.x();
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return - p.x() * p.x() * p.x();
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 6. * p.x();
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return - p.x() * p.x() * p.x();
        };

        if (submesh_flag)
        {
            /*
            disk::submesher<mesh_type> sm;

            size_t count = 0;
            for (auto& cl : msh)
            {
                std::cout << "Solving on element " << count << std::endl;
                std::stringstream ss;
                ss << "plot_element_" << count << ".dat";
                auto submesh = sm.generate_mesh(msh, cl);
                test_diffusion(submesh, f, sf, degree, ss.str());

                ss.str("");
                ss << "element_" << count << ".m";
                dump_to_matlab(submesh, ss.str());

                count++;
            }
            */
            disk::multiscale_local_problem<mesh_type> mlp(degree);

            for (auto& cl : msh)
            {
                mlp.assemble(msh, cl);
            }

        }
        else
        {
            test_diffusion(msh, f, sf, degree, "plot.dat");
        }
    }

    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<RealType, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        double m = 0.0;
        for (auto& cl : msh)
            m += measure(msh, cl);
        std::cout << m << std::endl;

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    if (std::regex_match(filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto f = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return M_PI * M_PI * sin(p.x() * M_PI);
            //return 1.0;
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
            //return -p.x() * p.x() * 0.5;
        };

        test_diffusion(msh, f, sf, degree, "plot.dat");
        //test_gradrec(msh, degree);
    }

    return 0;
}

#endif
