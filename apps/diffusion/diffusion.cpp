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

#include <map>

#include "../../config.h"

#ifdef HAVE_SOLVER_WRAPPERS
    #include "agmg/agmg.hpp"
#endif

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#if 0
template<typename MeshType>
void
test_gradrec(MeshType& msh, size_t degree)
{
    typedef MeshType                                   mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    auto f = [](const point<scalar_type, mesh_type::dimension>& p) -> auto {
        return sin(p.x() * M_PI);
        //return p.x();
    };

    disk::gradient_reconstruction<mesh_type,
                                        cell_basis_type,
                                        cell_quadrature_type,
                                        face_basis_type,
                                        face_quadrature_type> gradrec(degree);

    disk::projector<mesh_type,
                    cell_basis_type,
                    cell_quadrature_type,
                    face_basis_type,
                    face_quadrature_type> projk(degree);

    cell_basis_type cbk1(degree+1);

    cell_quadrature_type cq(2*degree+2);

    std::ofstream ofs("grad.dat");
    for (auto& cl : msh)
    {
        gradrec.compute(msh, cl);
        dynamic_vector<scalar_type> pk = projk.compute_whole(msh, cl, f);
        dynamic_vector<scalar_type> pk1 = gradrec.oper * pk;

        auto qps = cq.integrate(msh, cl);
        for (auto& qp : qps)
        {
            auto tp = qp.point();
            auto phi = cbk1.eval_functions(msh, cl, tp);
            scalar_type val = pk(0);
            for (size_t i = 1; i < phi.size(); i++)
                val += pk1(i-1) * phi[i];

            for (size_t i = 0; i < MeshType::dimension; i++)
                ofs << tp[i] << " ";
            ofs << val << " " << std::abs(val - f(tp)) << std::endl;
        }
    }
    ofs.close();
}
#endif

#if 0
template<typename MeshType>
class diffusion_solver
{
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef quadrature<mesh_type, cell_type>    cell_quadrature_type;
    typedef quadrature<mesh_type, face_type>    face_quadrature_type;

    typedef scaled_monomial_scalar_basis<mesh_type, cell_type>  cell_basis_type;
    typedef scaled_monomial_scalar_basis<mesh_type, face_type>  face_basis_type;

    typedef gradient_reconstruction<mesh_type, cell_basis_type,
                                    cell_quadrature_type, face_basis_type,
                                    face_quadrature_type> gradrec_type;

    typedef diffusion_like_stabilization<mesh_type, cell_basis_type,
                                         cell_quadrature_type, face_basis_type,
                                         face_quadrature_type> stab_type;

    typedef diffusion_like_static_condensation<mesh_type, cell_basis_type,
                                               cell_quadrature_type,
                                               face_basis_type,
                                               face_quadrature_type> statcond_type;

    typedef assembler<mesh_type, face_basis_type, face_quadrature_type> assembler_type;


    gradrec_type        m_gradrec;
    stab_type           m_stab;
    statcond_type       m_statcond;
    assembler_type      m_assembler;

    size_t              m_degree;

    template<typename Function, typename BoundaryConditions, typename MaterialParamFunction>
    void assemble_problem(const Function& load,
                          const BoundaryConditions& bcs,
                          const MaterialParamFunction& mpf)
    {
        timecounter_new tc;

        /* ASSEMBLE PROBLEM */
        std::cout << "Assembling..." << std::endl;
        tc.tic();
        for (auto& cl : msh)
        {
            auto mp = mpf(barycenter(msh, cl));

            gradrec.compute(msh, cl, mp);
            stab.compute(msh, cl, gradrec.oper);
            auto cell_rhs = compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, degree);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            auto scnp = statcond.compute(msh, cl, loc, cell_rhs);
            assembler.assemble(msh, cl, scnp);
        }

        std::cout << "  Imposing boundary conditions" << std::endl;
        assembler.impose_boundary_conditions(msh, solution);
        std::cout << "  Finalizing" << std::endl;
        assembler.finalize();
        tc.toc();
        std::cout << "Assembly total time: " << tc << " seconds." << std::endl;
    }

public:
    diffusion_solver() : m_degree(1)
    {}

    diffusion_solver(size_t degree) : m_degree(degree)
    {
        m_gradrec   = gradrec_type(m_degree);
        m_stab      = stab_type(m_degree);
        m_statcond  = statcond_type(m_degree);
        m_assembler = assembler_type(m_degree);
    }

};

#endif


struct timing_data
{
    size_t  mesh_elems;
    size_t  cell_degree, face_degree;
    size_t  solved_dofs;
    double  mesh_diameter;
    double  time_gradrec, time_statcond, time_stab, time_solver;
    double  l2_error;
};


template<typename MeshType, typename Function, typename Solution>
timing_data
test_diffusion(MeshType& msh,               /* handle to the mesh */
               const Function& load,        /* rhs */
               const Solution& solution,    /* solution of the problem */
               size_t degree,               /* degree of the method */
               const std::string& outfile)  /* output file name */
{
    typedef MeshType                                   mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    struct timing_data td;

    /*

    disk::gradient_reconstruction<mesh_type,
                                  cell_basis_type,
                                  cell_quadrature_type,
                                  face_basis_type,
                                  face_quadrature_type> gradrec(degree);

    */

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    int l = 0;
    size_t cell_degree = degree + l;
    size_t face_degree = degree;

    td.cell_degree = cell_degree;
    td.face_degree = face_degree;
    td.mesh_elems = msh.cells_size();

    std::cout << "Running HHO with cell degree " << cell_degree << " and face degree ";
    std::cout << face_degree << std::endl;

    bqdata_type bqd(cell_degree, face_degree);

    disk::gradient_reconstruction_bq<bqdata_type> gradrec(bqd);

    /*
    disk::diffusion_like_stabilization<mesh_type,
                                       cell_basis_type,
                                       cell_quadrature_type,
                                       face_basis_type,
                                       face_quadrature_type> stab(degree);
    */

    disk::diffusion_like_stabilization_bq<bqdata_type> stab(bqd);

    /*
    disk::diffusion_like_static_condensation<mesh_type,
                                             cell_basis_type,
                                             cell_quadrature_type,
                                             face_basis_type,
                                             face_quadrature_type> statcond(degree);
    */

    disk::diffusion_like_static_condensation_bq<bqdata_type> statcond(bqd);

    disk::assembler<mesh_type,
                    face_basis_type,
                    face_quadrature_type> assembler(msh, face_degree);

    timecounter_new tc;

    td.time_gradrec = 0.0;
    td.time_stab = 0.0;
    td.time_statcond = 0.0;
    td.time_solver = 0.0;

    /* ASSEMBLE PROBLEM */
    std::cout << "Assembling..." << std::endl;
    tc.tic();
    for (auto& cl : msh)
    {
        timecounter tc_detail;

        tc_detail.tic();
        gradrec.compute(msh, cl);
        tc_detail.toc();
        td.time_gradrec += tc_detail.to_double();

        tc_detail.tic();
        stab.compute(msh, cl, gradrec.oper);
        tc_detail.toc();
        td.time_stab += tc_detail.to_double();

        tc_detail.tic();
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, cell_degree);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        auto scnp = statcond.compute(msh, cl, loc, cell_rhs);
        tc_detail.toc();
        td.time_statcond += tc_detail.to_double();

        assembler.assemble(msh, cl, scnp);
    }

    assembler.impose_boundary_conditions(msh, solution);
    assembler.finalize();
    tc.toc();
    std::cout << "Assembly total time: " << tc << " seconds." << std::endl;

    td.solved_dofs = assembler.matrix.rows();

    return td;

    /* SOLVE */
    tc.tic();

#ifdef HAVE_INTEL_MKL
    Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

    size_t systsz = assembler.matrix.rows();
    size_t nnz = assembler.matrix.nonZeros();

    td.solved_dofs = systsz;

    std::cout << "Starting linear solver..." << std::endl;
    std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
    std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

    solver.analyzePattern(assembler.matrix);
    solver.factorize(assembler.matrix);
    dynamic_vector<scalar_type> X = solver.solve(assembler.rhs);

    tc.toc();
    std::cout << "Solver time: " << tc << " seconds." << std::endl;

    td.time_solver = tc.to_double();

    scalar_type diam = 0.0;
    scalar_type err_dof = 0.0;
    scalar_type err_fun = 0.0;

    //std::ofstream ofs(outfile);

    /*
    disk::projector<mesh_type,
                    cell_basis_type,
                    cell_quadrature_type,
                    face_basis_type,
                    face_quadrature_type> projk(degree);

    */
    disk::projector_bq<bqdata_type> projk(bqd);

    //cell_basis_type         cell_basis(degree);
    //cell_quadrature_type    cell_quadrature(2*degree);
    //face_basis_type         face_basis(degree);
    size_t fbs = bqd.face_basis.size();

    for (auto& cl : msh)
    {
        diam = std::max(diameter(msh, cl), diam);
        auto fcs = faces(msh, cl);
        auto num_faces = fcs.size();

        dynamic_vector<scalar_type> xFs = dynamic_vector<scalar_type>::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_id = eid.second;

            dynamic_vector<scalar_type> xF = dynamic_vector<scalar_type>::Zero(fbs);
            xF = X.block(face_id * fbs, 0, fbs, 1);
            xFs.block(face_i * fbs, 0, fbs, 1) = xF;
        }

        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, cell_degree);
        dynamic_vector<scalar_type> x = statcond.recover(msh, cl, loc, cell_rhs, xFs);

#if 0
        auto qps = bqd.cell_quadrature.integrate(msh, cl);
        for (auto& qp : qps)
        {
            auto phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());

            scalar_type pot = 0.0;
            for (size_t i = 0; i < bqd.cell_basis.range(0, cell_degree).size(); i++)
                pot += phi[i] * x(i);

            //auto potr = solution(qp.point());

            //scalar_type diff = 0.0;
            //diff = (pot - potr) * (pot - potr) * qp.weight();
            //std::cout << pot << " " << potr << " " << qp.weight() << " " << diff << std::endl;

            //err_fun += diff;

            auto tp = qp.point();
            for (size_t i = 0; i < MeshType::dimension; i++)
                ofs << tp[i] << " ";
            ofs << pot << " " << std::abs(pot - solution(tp)) << std::endl;
        }
#endif

        dynamic_vector<scalar_type> true_dof = projk.compute_cell(msh, cl, solution);
        dynamic_vector<scalar_type> comp_dof = x.block(0,0,true_dof.size(), 1);
        dynamic_vector<scalar_type> diff_dof = (true_dof - comp_dof);
        err_dof += diff_dof.dot(projk.cell_mm * diff_dof);
    }

    //ofs.close();

    std::cout << "Mesh diameter: " << diam << std::endl;
    std::cout << "L2-norm error, dof:   " << std::sqrt(err_dof) << std::endl;
    td.l2_error = std::sqrt(err_dof);

    return td;
}


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

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "triangle_spec");
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

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "triangle_gen");
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

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexagons_gen");
}

void test_hexahedra_specialized()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-2-2-2.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-4-4-4.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-8-8-8.hex");
    paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-16-16-16.hex");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/diskpp/testmesh-32-32-32.hex");

    typedef disk::cartesian_mesh<double, 3>         MeshType;
    typedef disk::cartesian_mesh_loader<double, 3>  LoaderType;

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexahedra_spec");
}

void test_hexahedra_generic()
{
    size_t runs = 5;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_2x2x2.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_4x4x4.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_8x8x8.msh");
    paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_16x16x16.msh");
    //paths.push_back("../../../diskpp/meshes/3D_hexa/fvca6/hexa_32x32x32.hex");

    typedef disk::generic_mesh<double, 3>       MeshType;
    typedef disk::fvca6_mesh_loader<double, 3>  LoaderType;

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "hexahedra_gen");
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

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "tetrahedra_spec");
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

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "tetrahedra_gen");
}

void test_polyhedra_generic()
{
    size_t runs = 2;

    std::vector<std::string> paths;
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_10.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_20.msh");
    paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_30.msh");
    //paths.push_back("../../../diskpp/meshes/3D_general/fvca6/dbls_40.msh");

    typedef disk::generic_mesh<double, 3>       MeshType;
    typedef disk::fvca6_mesh_loader<double, 3>  LoaderType;

    test_mesh_format<MeshType, LoaderType>(paths, runs, 0, 3, "polyhedra_gen");
}

int main(int argc, char **argv)
{
    //test_triangles_specialized();
    //test_triangles_generic();
    //test_hexagons_generic();
    //test_hexahedra_specialized();
    //test_hexahedra_generic();
    //test_tetrahedra_specialized();
    //test_tetrahedra_generic();
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
