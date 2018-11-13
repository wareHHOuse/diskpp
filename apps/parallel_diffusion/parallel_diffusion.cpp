/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
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

#include <iostream>
#include <regex>
#include <unistd.h>

#include <map>

#include "loaders/loader.hpp"
#include "hho/hho.hpp"

#include "timecounter.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <mpi.h>

class MPI_Context
{
    int m_rank, m_size;

public:
    MPI_Context(int *argc, char ***argv)
    {
        MPI_Init(argc, argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    }

    ~MPI_Context()
    {
        MPI_Finalize();
    }

    int comm_rank() const { return m_rank; }
    int comm_size() const { return m_size; }
};

template<typename Mesh>
class proc_mesh_partition
{
    int     m_comm_rank, m_comm_size;
    Mesh    m_msh;

public:
    proc_mesh_partition(const Mesh& msh, const MPI_Context& ctx)
    {
        m_comm_rank   = ctx.comm_rank();
        m_comm_size   = ctx.comm_size();
        m_msh       = msh;
    }

    auto begin() const
    {
        size_t rem_elem = m_msh.cells_size() % m_comm_size;
        if (rem_elem == 0)
        {
            size_t num_elem = m_msh.cells_size() / m_comm_size;
            return std::next(m_msh.cells_begin(), num_elem*m_comm_rank);
        }
        else
        {
            size_t num_elem = m_msh.cells_size() / (m_comm_size-1);
            return std::next(m_msh.cells_begin(), num_elem*m_comm_rank);
        }
    }

    auto end() const
    {
        size_t rem_elem = m_msh.cells_size() % m_comm_size;
        if (rem_elem == 0)
        {
            size_t num_elem = m_msh.cells_size() / m_comm_size;
            return std::next(m_msh.cells_begin(), num_elem*(m_comm_rank+1));
        }
        else
        {
            size_t num_elem = m_msh.cells_size() / (m_comm_size-1);
            if (m_comm_rank == m_comm_size - 1)
                return m_msh.cells_end();

            return std::next(m_msh.cells_begin(), num_elem*(m_comm_rank+1));
        }
    }
};

template<typename MeshType, typename Function, typename Solution>
void
test_diffusion(MeshType& msh,               /* handle to the mesh */
               const Function& load,        /* rhs */
               const Solution& solution,    /* solution of the problem */
               size_t degree,               /* degree of the method */
               const MPI_Context& ctx)
{
    proc_mesh_partition<MeshType> pmp(msh, ctx);

    typedef MeshType                                   mesh_type;
    typedef typename mesh_type::coordinate_type            scalar_type;
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


    double local_sum = 0.0, rank0_sum;
    for (auto itor = pmp.begin(); itor != pmp.end(); itor++)
        local_sum += measure(msh, *itor);

    MPI_Reduce(&local_sum, &rank0_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << "Local sum: " << local_sum << std::endl;
    if (ctx.comm_rank() == 0)
        std::cout << "Computed sum: " << rank0_sum << std::endl;


    int l = 0;
    size_t cell_degree = degree + l;
    size_t face_degree = degree;

    std::cout << "Running HHO with cell degree " << cell_degree << " and face degree ";
    std::cout << face_degree << std::endl;

    bqdata_type bqd(cell_degree, face_degree);

    disk::gradient_reconstruction_bq<bqdata_type> gradrec(bqd);
    disk::diffusion_like_stabilization_bq<bqdata_type> stab(bqd);
    disk::diffusion_like_static_condensation_bq<bqdata_type> statcond(bqd);

    //disk::assembler<mesh_type,
    //                face_basis_type,
    //                face_quadrature_type> assembler(msh, face_degree);

    timecounter_new tc;

    /* ASSEMBLE PROBLEM */
    std::cout << "Assembling..." << std::endl;

    tc.tic();
    for (auto itor = pmp.begin(); itor != pmp.end(); itor++)
    {
        auto cl = *itor;
        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        auto cell_rhs = disk::compute_rhs<cell_basis_type, cell_quadrature_type>(msh, cl, load, cell_degree);
        dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
        auto scnp = statcond.compute(msh, cl, loc, cell_rhs);
        //assembler.assemble(msh, cl, scnp);
    }

    //assembler.impose_boundary_conditions(msh, solution);
    //assembler.finalize();
    tc.toc();
    std::cout << "Assembly total time: " << tc << " seconds." << std::endl;

}

int main(int argc, char **argv)
{
    MPI_Context ctx(&argc, &argv);

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
            return M_PI * M_PI * sin(p.x() * M_PI);
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
        };

        test_diffusion(msh, f, sf, degree, ctx);

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
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
        };

        test_diffusion(msh, f, sf, degree, ctx);
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
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
        };

        if (submesh_flag)
        {
            //disk::multiscale_local_problem<mesh_type> mlp(degree);

            for (auto& cl : msh)
            {
                //mlp.assemble(msh, cl);
            }

        }
        else
        {
            test_diffusion(msh, f, sf, degree, ctx);
        }
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
        };

        test_diffusion(msh, f, sf, degree, ctx);
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
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
        };

        test_diffusion(msh, f, sf, degree, ctx);
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
        };

        auto sf = [](const point<RealType, mesh_type::dimension>& p) -> auto {
            return sin(p.x() * M_PI);
        };

        test_diffusion(msh, f, sf, degree, ctx);
    }

    return 0;
}
