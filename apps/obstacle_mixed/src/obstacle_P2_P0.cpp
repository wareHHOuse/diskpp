/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2016-2026
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 *
 * This file is copyright of the following author:
 * Guillaume Delay (C) 2026           guillaume.delay@sorbonne-universite.fr
 * Sorbonne Universite
 * Laboratoire Jacques-Louis Lions (LJLL)
 *
 */
/*
 * The content of this file corresponds to solving an obstacle problem with P2 elements
 * in the cells, P1 elements on the faces and
 * using P0 Lagrange multipliers.
 * This code does not use Lagrange basis and then can consider more general meshes.
 */

#include "gnuplot_output.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/methods/hho"

// For the mesh data structure
#include "diskpp/mesh/mesh.hpp"

// For the loaders and related helper functions
#include "diskpp/loaders/loader.hpp"

#include "diskpp/solvers/direct_solvers.hpp"

using namespace disk;
using disk::cells;
using disk::faces;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   ASSEMBLERS   ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template<typename Mesh>
class contact_assembler
{
    using T = typename Mesh::coordinate_type;
    typedef disk::BoundaryConditions<Mesh, true> boundary_type;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;
    hho_degree_info         di;
    std::vector<Triplet<T>> triplets;
    bool                    use_bnd;
    std::vector< Matrix<T, Dynamic, Dynamic> > loc_LHS;
    std::vector< Matrix<T, Dynamic, 1> >       loc_RHS, loc_N;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces, system_size;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T> LHS;
    vector_type     RHS;

    contact_assembler(const Mesh& msh, hho_degree_info hdi)
        : di(hdi), use_bnd(false)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        auto num_cells = msh.cells_size();
        loc_LHS.resize( num_cells );
        loc_RHS.resize( num_cells );

        const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
        system_size    = cbs * num_cells + fbs * num_other_faces + num_cells;

        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);

        // initialize the local constraint matrices (loc_N)
        loc_N.resize( num_cells );
        for (auto& cl : msh)
        {
            const auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto cell_offset = offset(msh, cl);
            vector_type mean_cb = vector_type::Zero(cbs); // mean value of the basis functions
            const auto qps = integrate(msh, cl, hdi.cell_degree());
            for (auto& qp : qps)
            {
                auto t_phi = cb.eval_functions( qp.point() );
                mean_cb += qp.weight() * t_phi;
            }
            // store this vector
            loc_N.at( cell_offset ) = mean_cb;
        }
    }

    void
    set_loc_mat(const Mesh&                     msh,
                const typename Mesh::cell_type& cl,
                const matrix_type&              lhs,
                const vector_type&              rhs)
    {
        auto cell_offset = offset(msh, cl);
        loc_LHS.at( cell_offset ) = lhs;
        loc_RHS.at( cell_offset ) = rhs;
    }

    template<typename Function>
    void
    assemble_mat(const Mesh&                     msh,
                 const typename Mesh::cell_type& cl,
                 const matrix_type&              lhs,
                 const Function&                 dirichlet_bf)
    {
        if(use_bnd)
            throw std::invalid_argument("contact_assembler: you have to use boundary type");

        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        const auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs + fcs.size() * fbs);

        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs;

        for (size_t i = 0; i < cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        vector_type dirichlet_data = vector_type::Zero(cbs + fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = offset(msh, fc);
            const auto face_LHS_offset = msh.cells_size() * cbs + compress_table.at(face_offset) * fbs;

            const bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                dirichlet_data.block(cbs + face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, fb, dirichlet_bf, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
            }
        }

    } // assemble_mat()


    template<typename Function>
    void
    assemble_rhs(const Mesh&                     msh,
                 const typename Mesh::cell_type& cl,
                 const matrix_type&              lhs,
                 const vector_type&              rhs,
                 const Function&                 dirichlet_bf)
    {
        if(use_bnd)
            throw std::invalid_argument("contact_assembler: you have to use boundary type");

        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        const auto fcs = faces(msh, cl);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs + fcs.size() * fbs);

        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs;

        for (size_t i = 0; i < cbs; i++)
                asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        vector_type dirichlet_data = vector_type::Zero(cbs + fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = offset(msh, fc);
            const auto face_LHS_offset = msh.cells_size() * cbs + compress_table.at(face_offset) * fbs;

            const bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                dirichlet_data.block(cbs + face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, fb, dirichlet_bf, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( !asm_map[j].assemble() )
                    RHS[ asm_map[i] ] -= lhs(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS[ asm_map[i] ] += rhs(i);
        }

    } // assemble_rhs()

    // init : set no contact constraint and assemble matrix
    // (first iteration matrix)
    template<typename Function>
    void
    init(const Mesh&                     msh,
         const Function&                 dirichlet_bf)
    {
        // assemble all local contributions for Laplacian part
        for (auto& cl : msh)
        {
            auto cell_offset = offset(msh, cl);
            assemble_mat(msh, cl, loc_LHS.at(cell_offset), dirichlet_bf);
            assemble_rhs(msh, cl, loc_LHS.at(cell_offset), loc_RHS.at(cell_offset), dirichlet_bf);
        }
        // assemble constraints (no contact)
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        auto mult_offset = cbs * msh.cells_size() + fbs * num_other_faces;
        for(size_t i = 0; i < msh.cells_size(); i++)
        {
            triplets.push_back( Triplet<T>(mult_offset + i, mult_offset + i, 1.0) );
        }

        // end assembly
        finalize();
    }

    // update_mat : assemble matrix according to the previous iteration solution
    template<typename Function>
    void
    update_mat(const Mesh&                     msh,
               const vector_type&              prev_sol,
               const Function&                 dirichlet_bf)
    {
        // assemble all local contributions for Laplacian part
        for (auto& cl : msh)
        {
            auto cell_offset = offset(msh, cl);
            assemble_mat(msh, cl, loc_LHS.at(cell_offset), dirichlet_bf);
        }
        // assemble constraints
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        auto mult_offset = cbs * msh.cells_size() + fbs * num_other_faces;

        for(auto& cl : msh)
        {
            auto cell_offset = offset(msh, cl);
            auto u_T = prev_sol.block(cell_offset*cbs, 0, cbs, 1);

            auto mean_cb = loc_N.at( cell_offset ); // mean value of the basis functions
            // T mean_u = priv::inner_product(mean_cb , u_T); // mean value of u_T over T
            T mean_u = 0.;
            for(size_t i=0; i<cbs; i++)
                mean_u += mean_cb(i,0) * u_T(i,0);
            auto sol_mult = prev_sol(mult_offset + cell_offset);

            if(mean_u < sol_mult)
            {
                // impose the solution to be of null-mean-value in this cell
                for(size_t i = 0; i < cbs; i++)
                {
                    triplets.push_back( Triplet<T>(mult_offset + cell_offset, cell_offset*cbs + i, mean_cb[i]) );
                }
            }
            else
            {
                // impose the multiplier to be zero in this cell
                triplets.push_back( Triplet<T>(mult_offset + cell_offset, mult_offset + cell_offset, 1.0) );
            }

            // coupling term
            for(size_t i = 0; i < cbs; i++)
            {
                triplets.push_back( Triplet<T>(cell_offset*cbs+i, mult_offset + cell_offset, -mean_cb[i]) );
            }
        }

        // end assembly
        finalize();
    }


    bool
    stop(const Mesh&         msh,
         const vector_type&  sol)
    {
        T TOL = 1e-14;

        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        auto mult_offset = cbs * msh.cells_size() + fbs * num_other_faces;

        bool ret = true;
        for(auto& cl : msh) {
            auto cell_offset = offset(msh, cl);
            auto u_T = sol.block(cell_offset*cbs, 0, cbs, 1);

            auto mean_cb = loc_N.at( cell_offset ); // mean value of the basis functions
            T mean_u = 0.;
            for(size_t i=0; i<cbs; i++)
                mean_u += mean_cb(i,0) * u_T(i,0);
            auto sol_mult = sol(mult_offset + cell_offset);

            if(mean_u < -TOL || sol_mult < -TOL)
            {
                ret = false;
                break;
            }
        }

        return ret;
    }



    template<typename Function>
    vector_type
    take_u(const Mesh& msh, const typename Mesh::cell_type& cl,
    const vector_type& solution, const Function& dirichlet_bf)
    {
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        auto cell_offset        = offset(msh, cl);
        auto cell_SOL_offset    = cell_offset * cbs;

        vector_type ret = vector_type::Zero(cbs + num_faces*fbs);
        ret.block(0, 0, cbs, 1) = solution.block(cell_SOL_offset, 0, cbs, 1);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
                return msh.is_boundary(fc);
            };

            bool dirichlet = is_dirichlet(fc);

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(cbs + face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = msh.cells_size() * cbs + compress_table.at(face_offset)*fbs;
                ret.block(cbs + face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }


    vector_type
    take_mult(const Mesh& msh, const typename Mesh::cell_type& cl,
              const vector_type& solution)
    {
        auto celdeg = di.cell_degree();
        auto facdeg = di.face_degree();
        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        auto mult_offset = cbs * msh.cells_size() + fbs * num_other_faces;
        auto cell_offset        = offset(msh, cl);
        auto cell_SOL_offset    = mult_offset + cell_offset;

        vector_type multT = solution.block(cell_SOL_offset, 0, 1, 1);

        return multT;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};
template<typename Mesh>
auto make_assembler_P2_P0(const Mesh& msh, hho_degree_info hdi)
{
    return contact_assembler<Mesh>(msh, hdi);
}


//////////
// Stuff for static condensation
//////////

template<typename Mesh, typename T>
auto
make_contact_SC(const Mesh&                                                      msh,
                const typename Mesh::cell_type&                                  cl,
                const hho_degree_info&                                           hdi,
                const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& D,
                const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& N_T)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto facdeg         = hdi.face_degree();
    const auto celdeg         = hdi.cell_degree();
    const auto num_face_dofs  = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto num_cell_dofs  = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fcs            = faces(msh, cl);
    const auto num_faces      = fcs.size();
    const auto num_faces_dofs = num_faces * num_face_dofs;
    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);
    assert(D.rows()   == 1);
    assert(D.cols()   == num_cell_dofs + 1);
    if ((rhs.size() != num_cell_dofs) && (rhs.size() != num_total_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);
    const matrix_type K_FT = lhs.bottomLeftCorner(num_faces_dofs, num_cell_dofs);
    const matrix_type K_FF = lhs.bottomRightCorner(num_faces_dofs, num_faces_dofs);

    assert(K_TT.cols() == num_cell_dofs);
    assert(K_TT.cols() + K_TF.cols() == lhs.cols());
    assert(K_TT.rows() + K_FT.rows() == lhs.rows());
    assert(K_TF.rows() + K_FF.rows() == lhs.rows());
    assert(K_FT.cols() + K_FF.cols() == lhs.cols());

    const vector_type cell_rhs  = rhs.head(num_cell_dofs);
    vector_type       faces_rhs = vector_type::Zero(num_faces_dofs);

    if (rhs.size() == num_total_dofs)
    {
        faces_rhs = rhs.tail(num_faces_dofs);
    }

    const matrix_type D_TT = D.topLeftCorner(1, num_cell_dofs);
    const matrix_type D_lam = D.topRightCorner(1, 1);

    const auto K_TT_ldlt = K_TT.ldlt();
    if (K_TT_ldlt.info() != Eigen::Success)
    {
        throw std::invalid_argument("static condensation: K_TT is not positive definite");
    }

    const matrix_type AL = K_TT_ldlt.solve(K_TF);
    const vector_type bL = K_TT_ldlt.solve(cell_rhs);

    const auto ID = matrix_type::Identity(num_cell_dofs, num_cell_dofs);
    const auto K_TT_inv = K_TT_ldlt.solve(ID);

    const auto E = D_TT * K_TT_inv * N_T + D_lam;
    const auto E_inv = E.inverse();

    const matrix_type E2 = E_inv * D_TT;
    const matrix_type E2_bis = N_T * E2;
    const matrix_type E3 = K_TT_ldlt.solve(E2_bis);

    const matrix_type AC = K_FF - K_FT * AL + K_FT * E3 * AL;
    const vector_type bC = faces_rhs - K_FT * bL + K_FT * E3 * bL;

    return std::make_pair(AC, bC);
}


///////////////////
template<typename Mesh, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
contact_static_decondensation(const Mesh&                                                      msh,
                              const typename Mesh::cell_type&                                  cl,
                              const hho_degree_info&                                           hdi,
                              const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& lhs,
                              const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              rhs,
                              const typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& D,
                              const typename Eigen::Matrix<T, Eigen::Dynamic, 1>& N_T,
                              const typename Eigen::Matrix<T, Eigen::Dynamic, 1>&              solF)
{
    using matrix_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    const auto facdeg         = hdi.face_degree();
    const auto celdeg         = hdi.cell_degree();
    const auto num_face_dofs  = scalar_basis_size(facdeg, Mesh::dimension - 1);
    const auto num_cell_dofs  = scalar_basis_size(celdeg, Mesh::dimension);
    const auto fcs            = faces(msh, cl);
    const auto num_faces      = fcs.size();
    const auto num_faces_dofs = num_faces * num_face_dofs;
    const auto num_total_dofs = num_cell_dofs + num_faces_dofs;

    assert(lhs.rows() == lhs.cols());
    assert(lhs.cols() == num_total_dofs);

    if ((rhs.size() < num_cell_dofs))
    {
        throw std::invalid_argument("static condensation: incorrect size of the rhs");
    }

    const matrix_type K_TT = lhs.topLeftCorner(num_cell_dofs, num_cell_dofs);
    const matrix_type K_TF = lhs.topRightCorner(num_cell_dofs, num_faces_dofs);

    const matrix_type D_TT = D.topLeftCorner(1, num_cell_dofs);
    const matrix_type D_lam = D.topRightCorner(1, 1);

    vector_type uF = solF.head(num_faces_dofs);

    const auto K_TT_ldlt = K_TT.ldlt();
    if (K_TT_ldlt.info() != Eigen::Success)
    {
        throw std::invalid_argument("static condensation: K_TT is not positive definite");
    }

    const auto ID = matrix_type::Identity(num_cell_dofs, num_cell_dofs);
    const auto K_TT_inv = K_TT_ldlt.solve(ID);

    const auto E_inv = (D_TT * K_TT_inv * N_T + D_lam).inverse();
    const matrix_type E2 = E_inv * D_TT;

    const vector_type multT = - E2 * K_TT.ldlt().solve(rhs.head(num_cell_dofs) - K_TF * uF);

    const vector_type solT = K_TT_ldlt.solve(rhs.head(num_cell_dofs) - K_TF * uF + N_T * multT);

    vector_type ret          = vector_type::Zero(num_total_dofs + 1);
    ret.head(num_cell_dofs)  = solT;
    ret.block(num_cell_dofs, 0, num_faces_dofs, 1) = uF;
    ret.block(num_total_dofs, 0, 1, 1) = multT;

    return ret;
}

//////////
/* Assembler using static condensation
 * This algorithm is described in Algorithm 2 from the work
 * Jad Dabaghi, Guillaume Delay
 * A unified framework for high-order numerical discretizations of variational inequalities
 * Computers and Mathematics with Applications 92 (2021) 62-75
 */
//////////

template<typename Mesh>
class contact_condensed_assembler
{
    using T = typename Mesh::coordinate_type;
    typedef disk::BoundaryConditions<Mesh, true> boundary_type;

    std::vector<size_t>     compress_table;
    std::vector<size_t>     expand_table;
    hho_degree_info         di;
    std::vector<Triplet<T>> triplets;
    bool                    use_bnd;
    std::vector< Matrix<T, Dynamic, Dynamic> > loc_LHS;
    std::vector< Matrix<T, Dynamic, 1> >       loc_RHS, loc_N;
    std::vector< bool >     active_constr;

    size_t num_all_faces, num_dirichlet_faces, num_other_faces, system_size;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    SparseMatrix<T> LHS;
    vector_type     RHS;

    contact_condensed_assembler(const Mesh& msh, hho_degree_info hdi)
        : di(hdi), use_bnd(false)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        num_all_faces       = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces     = num_all_faces - num_dirichlet_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            const auto fc = *std::next(msh.faces_begin(), i);
            if (!is_dirichlet(fc))
            {
                compress_table.at(i)               = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        auto num_cells = msh.cells_size();
        loc_LHS.resize( num_cells );
        loc_RHS.resize( num_cells );

        const auto fbs = scalar_basis_size(hdi.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
        system_size = fbs * num_other_faces;

        active_constr.resize(num_cells);

        LHS = SparseMatrix<T>(system_size, system_size);
        RHS = vector_type::Zero(system_size);

        // initialize the local constraint matrices (loc_N)
        loc_N.resize( num_cells );
        for (auto& cl : msh)
        {
            const auto cb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
            auto cell_offset = offset(msh, cl);
            vector_type mean_cb = vector_type::Zero(cbs); // mean value of the basis functions
            const auto qps = integrate(msh, cl, hdi.cell_degree());
            for (auto& qp : qps)
            {
                auto t_phi = cb.eval_functions( qp.point() );
                mean_cb += qp.weight() * t_phi;
            }
            // store this vector
            loc_N.at( cell_offset ) = mean_cb;
        }
    }

    void
    set_loc_mat(const Mesh&                     msh,
                const typename Mesh::cell_type& cl,
                const matrix_type&              lhs,
                const vector_type&              rhs)
    {
        auto cell_offset = offset(msh, cl);
        loc_LHS.at( cell_offset ) = lhs;
        loc_RHS.at( cell_offset ) = rhs;
        active_constr.at( cell_offset ) = false;
    }

    template<typename Function>
    void
    assemble_contrib(const Mesh&                     msh,
                     const typename Mesh::cell_type& cl,
                     const matrix_type&              lhs,
                     const vector_type&              rhs,
                     const vector_type&              solF,
                     const Function&                 dirichlet_bf)
    {
        if(use_bnd)
            throw std::invalid_argument("contact_condensed_assembler: you have to use boundary type");

        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        auto cell_offset = offset(msh, cl);
        auto mean_cb = loc_N.at(cell_offset);

        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);
        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        const auto fcs = faces(msh, cl);

        matrix_type D = matrix_type::Zero(1, cbs+1); // matrix of constraints in the present cell
        if( active_constr.at( cell_offset ) )
        {
            for(size_t i = 0; i < cbs; i++)
                D(0,i) = mean_cb(i);
        }
        else
        {
            D(0,cbs) = 1.;
        }

        auto loc_sol = contact_static_decondensation(msh, cl, di, loc_LHS.at( cell_offset ), loc_RHS.at( cell_offset ), D, mean_cb, solF);

        vector_type solT = loc_sol.head(cbs);
        T multT = loc_sol(cbs + fcs.size() * fbs);

        D = matrix_type::Zero(1, cbs + 1);

        T mean_u = 0.;
        for(size_t i = 0; i < cbs; i++)
            mean_u += mean_cb(i,0) * solT(i,0);

        if(mean_u < multT)
        {
            // impose the solution to be of null-mean-value in this cell
            active_constr.at( cell_offset ) = true;
            for(size_t i = 0; i < cbs; i++)
            {
                D(0,i) = mean_cb(i);
            }
        }
        else
        {
            active_constr.at( cell_offset ) = false;
            D(0,cbs) = 1.;
        }

        auto SC = make_contact_SC(msh, cl, di, lhs, rhs, D, mean_cb);
        matrix_type lhs_sc = SC.first;
        vector_type rhs_sc = SC.second;

        std::vector<assembly_index> asm_map;
        asm_map.reserve(fcs.size() * fbs);

        vector_type dirichlet_data = vector_type::Zero(fcs.size()*fbs);

        for (size_t face_i = 0; face_i < fcs.size(); face_i++)
        {
            const auto fc              = fcs[face_i];
            const auto face_offset     = offset(msh, fc);
            const auto face_LHS_offset = compress_table.at(face_offset) * fbs;

            const bool dirichlet = is_dirichlet(fc);

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());
                dirichlet_data.block(face_i * fbs, 0, fbs, 1) =
                  project_function(msh, fc, fb, dirichlet_bf, di.face_degree());
            }
        }

        for (size_t i = 0; i < lhs_sc.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs_sc.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs_sc(i,j)) );
                else
                    RHS[ asm_map[i] ] -= lhs_sc(i,j) * dirichlet_data(j);
            }
        }

        for (size_t i = 0; i < rhs_sc.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;
            RHS[ asm_map[i] ] += rhs_sc(i);
        }

    } // assemble_contrib()


    // init : set no contact constraint and assemble matrix
    // (matrix for the first iteration)
    template<typename Function>
    void
    init(const Mesh&                     msh,
         const Function&                 dirichlet_bf)
    {
        const auto facdeg  = di.face_degree();
        const auto fbs = scalar_basis_size(facdeg, Mesh::dimension - 1);

        // assemble all local contributions for Laplacian part
        for (auto& cl : msh)
        {
            // init solution with no contact
            auto num_faces = howmany_faces(msh, cl);
            auto num_dofs = num_faces * fbs;
            auto cell_offset = offset(msh, cl);

            vector_type solF = vector_type::Zero(num_dofs);
            for(size_t i=0; i<num_dofs; i++)
                solF(i) = 1.0;

            assemble_contrib(msh, cl, loc_LHS.at(cell_offset), loc_RHS.at(cell_offset),solF, dirichlet_bf);
        }

        // end assembly
        finalize();
    }


    template<typename Function>
    vector_type
    get_solF(const Mesh& msh, const typename Mesh::cell_type& cl,
             const vector_type& solution, const Function& dirichlet_bf)
    {
        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto fcs = faces(msh, cl);

        auto num_faces = fcs.size();

        vector_type ret = vector_type::Zero(num_faces*fbs);

        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];

            auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
                return msh.is_boundary(fc);
            };

            bool dirichlet = is_dirichlet(fc);

            if (dirichlet)
            {
                auto fb = make_scalar_monomial_basis(msh, fc, di.face_degree());

                matrix_type mass = make_mass_matrix(msh, fc, fb, di.face_degree());
                vector_type rhs = make_rhs(msh, fc, fb, dirichlet_bf, di.face_degree());
                ret.block(face_i*fbs, 0, fbs, 1) = mass.ldlt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = compress_table.at(face_offset)*fbs;
                ret.block(face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    // update_mat : assemble matrix according to the previous iteration solution
    template<typename Function>
    void
    update_mat(const Mesh&                     msh,
               const vector_type&              prev_sol,
               const Function&                 dirichlet_bf)
    {
        // clear RHS
        RHS = vector_type::Zero(system_size);

        // assemble all local contributions for Laplacian part
        for (auto& cl : msh)
        {
            auto solF = get_solF(msh, cl, prev_sol, dirichlet_bf);

            auto cell_offset = offset(msh, cl);
            assemble_contrib(msh, cl, loc_LHS.at(cell_offset), loc_RHS.at(cell_offset),solF, dirichlet_bf);
        }

        finalize();
    }

    template<typename Function>
    bool
    stop(const Mesh&         msh,
         const vector_type&  sol, const Function& dirichlet_bf)
    {
        T TOL = 1e-14;

        const auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension - 1);
        const auto cbs = scalar_basis_size(di.cell_degree(), Mesh::dimension);

        // test the cells
        for (auto& cl : msh)
        {
            auto cell_offset = offset(msh, cl);
            auto mean_cb = loc_N.at(cell_offset);
            matrix_type D = matrix_type::Zero(1, cbs+1); // matrix of constraints in the present cell
            if( active_constr.at( cell_offset ) )
            {
                for(size_t i = 0; i < cbs; i++)
                    D(0,i) = mean_cb(i);
            }
            else
            {
                D(0,cbs) = 1.;
            }

            auto solF = get_solF(msh, cl, sol, dirichlet_bf);

            auto loc_sol = contact_static_decondensation(msh, cl, di, loc_LHS.at( cell_offset ), loc_RHS.at( cell_offset ), D, mean_cb, solF);

            const auto fcs = faces(msh, cl);
            vector_type solT = loc_sol.head(cbs);
            T multT = loc_sol(cbs + fcs.size() * fbs);

            T mean_u = 0.;
            for(size_t i = 0; i < cbs; i++)
                mean_u += mean_cb(i,0) * solT(i,0);

            if(mean_u < -TOL || multT < -TOL)
                return false;
        }

        return true;
    }



    template<typename Function>
    vector_type
    take_u(const Mesh& msh, const typename Mesh::cell_type& cl,
    const vector_type& solution, const Function& dirichlet_bf)
    {
        auto solF = get_solF(msh, cl, solution, dirichlet_bf);

        auto cell_offset        = offset(msh, cl);
        auto celdeg = di.cell_degree();
        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);

        auto mean_cb = loc_N.at(cell_offset);
        matrix_type D = matrix_type::Zero(1, cbs+1); // matrix of constraints in the present cell
        if( active_constr.at( cell_offset ) )
        {
            for(size_t i = 0; i < cbs; i++)
                D(0,i) = mean_cb(i);
        }
        else
        {
            D(0,cbs) = 1.;
        }

        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto num_faces = howmany_faces(msh, cl);

        vector_type full_sol = contact_static_decondensation(msh, cl, di, loc_LHS.at( cell_offset ), loc_RHS.at( cell_offset ), D, mean_cb, solF);

        return full_sol.head(cbs + num_faces * fbs);
    }

    template<typename Function>
    vector_type
    take_mult(const Mesh& msh, const typename Mesh::cell_type& cl,
              const vector_type& solution, const Function& dirichlet_bf)
    {
        auto solF = get_solF(msh, cl, solution, dirichlet_bf);

        auto cell_offset        = offset(msh, cl);
        auto celdeg = di.cell_degree();
        auto cbs = scalar_basis_size(celdeg, Mesh::dimension);
        const auto cb = make_scalar_monomial_basis(msh, cl, celdeg);

        auto mean_cb = loc_N.at(cell_offset);
        matrix_type D = matrix_type::Zero(1, cbs+1); // matrix of constraints in the present cell
        if( active_constr.at( cell_offset ) )
        {
            for(size_t i = 0; i < cbs; i++)
                D(0,i) = mean_cb(i);
        }
        else
        {
            D(0,cbs) = 1.;
        }

        auto full_sol = contact_static_decondensation(msh, cl, di, loc_LHS.at( cell_offset ), loc_RHS.at( cell_offset ), D, mean_cb, solF);

        auto facdeg = di.face_degree();
        auto fbs = scalar_basis_size(di.face_degree(), Mesh::dimension-1);
        auto num_faces = howmany_faces(msh, cl);

        auto multT = full_sol.block( cbs + num_faces * fbs, 0, 1, 1);
        return multT;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();

        dump_sparse_matrix(LHS, "diff.dat");
    }

    size_t num_assembled_faces() const
    {
        return num_other_faces;
    }

};
template<typename Mesh>
auto make_condensed_assembler_P2_P0(const Mesh& msh, hho_degree_info hdi)
{
    return contact_condensed_assembler<Mesh>(msh, hdi);
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CONTACT SOLVER   ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/////// test_info -> for error output
template<typename T>
class test_info {
public:
    test_info() : H1_u(0.) , L2_u(0.) , L2_mult(0.) , nb_dof(0), nb_Newton_iter(0) {}
    T H1_u; // H1-error for u
    T L2_u; // L2-error for u
    T L2_mult; // L2-error for mult
    size_t nb_dof; // number of degrees of freedom of the linear system
    size_t nb_Newton_iter; // number of iterations for the Newton loop
    T h_max; // maximal cell diameter
};

/////// Contact solver
template<typename Mesh>
test_info<typename Mesh::coordinate_type>
run_contact_solver(const Mesh& msh, size_t degree)
{
    using T = typename Mesh::coordinate_type;
    using point_type = typename Mesh::point_type;

    hho_degree_info hdi(degree+1, degree);

#if 1
    auto rhs_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return 8.0 * R2 - 16.0 * r2;
        else
            return - 8.0 * R2;
    };
    auto sol_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return (r2 - R2) * (r2 - R2);
        else
            return 0.0;
    };
    auto sol_grad = [](const point_type& pt) -> auto {
        Matrix<T, 1, 2> ret;
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
        {
            T coeff = 2.0*2.0*(r2 - R2);
            ret(0) =  coeff * x1;
            ret(1) =  coeff * y1;
        }
        else
        {
            ret(0) = 0.0;
            ret(1) = 0.0;
        }
        return ret;
    };
    auto mult_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return 0.0;
        else
            return 8.0 * R2;
    };
#else
    auto rhs_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return -24.0*(r2-R2)*(r2-R2)*(r2-R2)*(r2-R2)*(6.0*r2-R2);
        else
            return -r2*(R2-r2)*(R2-r2)*(R2-r2);
    };
    auto sol_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return (r2 - R2) * (r2 - R2) * (r2 - R2) * (r2 - R2) * (r2 - R2) * (r2 - R2);
        else
            return 0.0;
    };
    auto sol_grad = [](const point_type& pt) -> auto {
        Matrix<T, 1, 2> ret;
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
        {
            T coeff = 2.0*6.0*(r2 - R2) * (r2 - R2) * (r2 - R2) * (r2 - R2) * (r2 - R2);
            ret(0) =  coeff * x1;
            ret(1) =  coeff * y1;
        }
        else
        {
            ret(0) = 0.0;
            ret(1) = 0.0;
        }

        return ret;
    };
    auto mult_fun = [](const point_type& pt) -> T {
        auto x1 = pt.x() - 0.5;
        auto y1 = pt.y() - 0.5;
        auto r2 = x1*x1 + y1*y1;
        auto R2 = 1.0 / 9.0;
        if(r2 > R2)
            return 0.0;
        else
            return r2*(R2-r2)*(R2-r2)*(R2-r2);
    };
#endif

    auto assembler_sc = make_condensed_assembler_P2_P0(msh, hdi);
    auto assembler = make_assembler_P2_P0(msh, hdi);
    test_info<double> TI;

    bool scond = true; // static condensation

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_vector_hho_gradrec(msh, cl, hdi);
        auto stab   = make_scalar_hdg_stabilization(msh, cl, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;
        if(scond)
            assembler_sc.set_loc_mat(msh, cl, A, rhs);
        else
        {
            assembler.set_loc_mat(msh, cl, A, rhs);
        }
    }


    if(scond)
        std::cout << "end assembly : nb dof = " << assembler_sc.RHS.size() << std::endl;
    else
        std::cout << "end assembly : nb dof = " << assembler.RHS.size() << std::endl;


    size_t systsz, nnz;
    dynamic_vector<T> sol;
    size_t Newton_iter = 0;
    bool stop_loop = false;
    size_t max_Newton_iter = 300;

    // Newton loop
    while(!stop_loop and Newton_iter < max_Newton_iter)
    {
        if(scond)
        {
            if(Newton_iter == 0)
                assembler_sc.init(msh, sol_fun);
            else
                assembler_sc.update_mat(msh, sol, sol_fun);
            systsz = assembler_sc.LHS.rows();
        }
        else
        {
            if(Newton_iter == 0)
                assembler.init(msh, sol_fun);
            else
                assembler.update_mat(msh, sol, sol_fun);
            systsz = assembler.LHS.rows();
        }
        sol = dynamic_vector<T>::Zero(systsz);

        if(scond)
        {
            std::cout << "Running solver..." << std::flush;
            auto status = disk::solvers::sparse_lu(assembler_sc.LHS, assembler_sc.RHS,
                                                   sol, disk::solvers::direct_solver::sparselu);
            if (status != disk::solvers::direct_solver_status::ok) {
                std::cout << "LU factorization failed" << std::endl;
                return TI;
            }
            std::cout << "done" << std::endl;
        }
        else
        {
            std::cout << "Running solver..." << std::flush;
            auto status = disk::solvers::sparse_lu(assembler.LHS, assembler.RHS,
                                                   sol, disk::solvers::direct_solver::sparselu);
            // auto status = disk::solvers::sparse_lu(assembler.LHS, assembler.RHS, sol);
            if (status != disk::solvers::direct_solver_status::ok) {
                std::cout << "LU factorization failed" << std::endl;
                return TI;
            }
            std::cout << "done" << std::endl;
        }

        Newton_iter++;
        std::cout << "end Newton iter nb " << Newton_iter << std::endl;

        if(scond)
        {
            stop_loop = assembler_sc.stop(msh, sol, sol_fun);
        }

        else
            stop_loop = assembler.stop(msh, sol);
    } // Newton loop

    std::cout << "Start post-process" << std::endl;

    T u_H1_error = 0.0;
    T u_L2_error = 0.0;
    T mult_L2_error = 0.0;

    postprocess_output<T>  postoutput;

    auto uT_gp  = std::make_shared< gnuplot_output_object<T> >("uT.dat");
    auto multT_gp  = std::make_shared< gnuplot_output_object<T> >("multT.dat");

    for (auto& cl : msh)
    {
        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto cbs = cb.size();

        Eigen::Matrix<T, Eigen::Dynamic, 1> realsol = project_function(msh, cl, cb, sol_fun, 2);
        Eigen::Matrix<T, Eigen::Dynamic, 1> fullsol, mult_sol;
        auto gr     = make_vector_hho_gradrec(msh, cl, hdi);
        auto stab   = make_scalar_hdg_stabilization(msh, cl, hdi);
        auto rhs    = make_rhs(msh, cl, cb, rhs_fun);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = gr.second + stab;

        if(scond)
        {
            fullsol = assembler_sc.take_u(msh, cl, sol, sol_fun);
            mult_sol = assembler_sc.take_mult(msh, cl, sol, sol_fun);
        }
        else
        {
            fullsol = assembler.take_u(msh, cl, sol, sol_fun);
            mult_sol = assembler.take_mult(msh, cl, sol);
        }

        auto cell_dofs = fullsol.head( cb.size() );

        // errors
        const auto celdeg = hdi.cell_degree();
        const auto qps = integrate(msh, cl, 2*celdeg);
        for (auto& qp : qps)
        {
            auto grad_ref = sol_grad( qp.point() );
            auto t_dphi = cb.eval_gradients( qp.point() );
            Matrix<T, 1, 2> grad = Matrix<T, 1, 2>::Zero();

            for (size_t i = 0; i < cbs; i++ )
                grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);

            // H1-error
            u_H1_error += qp.weight() * (grad_ref - grad).dot(grad_ref - grad);

            // L2-error
            auto t_phi = cb.eval_functions( qp.point() );
            T v = cell_dofs.dot( t_phi );
            u_L2_error += qp.weight() * (sol_fun(qp.point()) - v) * (sol_fun(qp.point()) - v);

            // mult-L2-error
            T mult = mult_sol(0);
            T mult_ex = mult_fun(qp.point());
            mult_L2_error += qp.weight() * (mult_ex - mult) * (mult_ex - mult);
        }

        // gnuplot output for cells
        auto pts = points(msh, cl);
        for(size_t i=0; i < pts.size(); i++)
        {
            T sol_uT = cell_dofs.dot( cb.eval_functions( pts[i] ) );
            uT_gp->add_data( pts[i], sol_uT );
            // T sol_multT = mult_sol.dot( cb.eval_functions(pts[i]) );
            // multT_gp->add_data( pts[i], sol_multT );
        }
        multT_gp->add_data( barycenter(msh,cl), mult_sol(0) );
    }

    postoutput.add_object(uT_gp);
    postoutput.add_object(multT_gp);
    postoutput.write();

    std::cout << "ended run : H1-error is " << std::sqrt(u_H1_error) << std::endl;
    std::cout << "            L2-error is " << std::sqrt(u_L2_error) << std::endl;
    std::cout << "            mult-L2-error is " << std::sqrt(mult_L2_error) << std::endl;

    TI.nb_dof = systsz;
    TI.H1_u = std::sqrt(u_H1_error);
    TI.L2_u = std::sqrt(u_L2_error);
    TI.L2_mult = std::sqrt(mult_L2_error);
    TI.nb_Newton_iter = Newton_iter;
    TI.h_max = average_diameter(msh);

    return TI;
}

//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   TESTS AUTO   //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


template<typename T>
void
tests_auto_2d()
{
    typedef disk::simplicial_mesh<T, 2>  mesh_type;
    // typedef disk::generic_mesh<double, 2> mesh_type;

    // list of mesh files
    std::vector<std::string> meshes;
    meshes.push_back("../../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d");
    meshes.push_back("../../../../diskpp/meshes/2D_triangles/netgen/tri02.mesh2d");
    meshes.push_back("../../../../diskpp/meshes/2D_triangles/netgen/tri03.mesh2d");
    meshes.push_back("../../../../diskpp/meshes/2D_triangles/netgen/tri04.mesh2d");
    meshes.push_back("../../../../diskpp/meshes/2D_triangles/netgen/tri05.mesh2d");

    // meshes.push_back("../../../../diskpp/meshes/2D_hex/fvca5/hexagonal_1.typ1");
    // meshes.push_back("../../../../diskpp/meshes/2D_hex/fvca5/hexagonal_2.typ1");
    // meshes.push_back("../../../../diskpp/meshes/2D_hex/fvca5/hexagonal_3.typ1");
    // meshes.push_back("../../../../diskpp/meshes/2D_hex/fvca5/hexagonal_4.typ1");
    // meshes.push_back("../../../../diskpp/meshes/2D_hex/fvca5/hexagonal_5.typ1");


    size_t nb_meshes = meshes.size();

    // list of export files
    std::vector<std::string> files;
    files.push_back("./test_k0.txt");
    files.push_back("./test_P2_P0.txt");
    files.push_back("./test_k2.txt");
    files.push_back("./test_k3.txt");

    // we test degree 1
    for(int degree=1; degree <= 1; degree++)
    {
        std::cout << " WORKING WITH k = " << degree << std::endl;

        // open the output file
        std::ofstream file;
        file.open (files.at(degree), std::ios::in | std::ios::trunc);
        if (!file.is_open())
            throw std::logic_error("file not open");

        // init the file
        file << "N\tL2_u\tH1_u\tL2_mult\tdof\tnb_Newton_iter\th" << std::endl;

        // we test all the meshes in the list
        for(size_t i=0; i < nb_meshes; i++)
        {
            mesh_type msh;

            disk::netgen_mesh_loader<T, 2> loader;
            // disk::fvca5_mesh_loader<T, 2> loader;

            if( !loader.read_mesh(meshes.at(i)) )
                std::cout << "error loading mesh !" << std::endl;
            loader.populate_mesh(msh);

            // test this mesh
            auto TI = run_contact_solver(msh, degree);

            // write the results in the file
            file << i+1 << "\t" << TI.L2_u << "\t" << TI.H1_u << "\t"
                 << TI.L2_mult
                 << "\t" << TI.nb_dof << "\t" << TI.nb_Newton_iter
                 << "\t" << TI.h_max
                 << std::endl;
        }

        // close the file
        file.close();
    }
}



//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////   MAIN   /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/* run main with :
   ./obstacle_P2_P0
*/
int main(int argc, char **argv)
{
    tests_auto_2d<double>();
    return 0;
}
