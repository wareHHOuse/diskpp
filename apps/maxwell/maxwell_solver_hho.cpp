/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020, 2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "bases/bases.hpp"
#include "quadratures/quadratures.hpp"
#include "methods/hho"
#include "methods/implementation_hho/curl.hpp"
#include "core/loaders/loader.hpp"
#include "output/silo.hpp"
#include "solvers/solver.hpp"
#include "solvers/mumps.hpp"
#include "contrib/sol2/include/sol/sol.hpp"

#define EPS_0   (8.8541878128e-12)
#define MU_0    (4e-7*M_PI)

template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
class maxwell_assembler_condensed
{
    using CoordT = typename Mesh::coordinate_type;

    Mesh                                msh;
    disk::hho_degree_info               chdi;

    std::vector<Triplet<ScalT>>         triplets;
    std::vector<bool>                   is_dirichlet;

    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;

    const size_t INVALID_OFFSET = (size_t) ~0;

    size_t face_basis_size() const
    {
        if (Mesh::dimension == 3)
            return disk::vector_basis_size(chdi.face_degree(), 2, 2);
        if (Mesh::dimension == 2)
            return disk::scalar_basis_size(chdi.face_degree(), 1);
        return 0;
    }

    /* Get the offset of a face in the linear system */
    size_t get_system_offset(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(compress_table.size() == msh.faces_size());
        auto cnum = compress_table[face_num];
        assert(cnum != INVALID_OFFSET);
        auto fbs = face_basis_size();
        return cnum*fbs;
    }

    /* Determine if a face should be assembled */
    bool is_in_system(const Mesh& msh,
        const typename Mesh::face_type& fc)
    {
        auto bi = msh.boundary_info(fc);
        auto face_num = offset(msh, fc);
        assert(face_num < msh.faces_size());
        assert(is_dirichlet.size() == msh.faces_size());
        return not (bi.is_boundary() and is_dirichlet[face_num]);
    }

public:
    
    SparseMatrix<ScalT>         LHS;
    Matrix<ScalT, Dynamic, 1>   RHS;

    size_t                  syssz;
    size_t                  sysfcs;

    maxwell_assembler_condensed(const Mesh& pmsh,
                      const disk::hho_degree_info& pchdi,
                      const std::vector<bool> p_is_dirichlet)
        : msh(pmsh), chdi(pchdi), is_dirichlet(p_is_dirichlet)
    {
        auto fbs = face_basis_size();

        using face_type = typename Mesh::face_type;
        auto in_system = [&](const face_type& fc) -> bool {
            auto ofs = offset(msh, fc);
            return not (msh.is_boundary(fc) and is_dirichlet.at(ofs));
        };

        sysfcs = std::count_if(msh.faces_begin(), msh.faces_end(), in_system);
        syssz = fbs*sysfcs;

        LHS = SparseMatrix<ScalT>(syssz, syssz);
        RHS = Matrix<ScalT, Dynamic, 1>::Zero(syssz);

        compress_table.resize( msh.faces_size(), INVALID_OFFSET);
        expand_table.resize( sysfcs );

        size_t face_i = 0;
        size_t compressed_ofs = 0;
        for (auto& fc : faces(msh))
        {
            assert(compressed_ofs <= face_i);
            if ( is_in_system(msh, fc) )
            {
                assert(face_i < compress_table.size());
                compress_table[face_i] = compressed_ofs;
                assert(compressed_ofs < expand_table.size());
                expand_table[compressed_ofs] = face_i;
                compressed_ofs++;
            }

            face_i++;
        }

        assert(face_i == msh.faces_size());
    }

    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<ScalT, Dynamic, Dynamic>& lhsc,
             const Matrix<ScalT, Dynamic, 1>& rhs)
    {
        auto fbs = face_basis_size();

        auto fcs = faces(msh, cl);
        for (size_t fi = 0; fi < fcs.size(); fi++)
        {
            if ( not is_in_system(msh, fcs[fi]) )
                continue;

            auto cofsi = get_system_offset(msh, fcs[fi]);
            for (size_t fj = 0; fj < fcs.size(); fj++)
            {
                if ( not is_in_system(msh, fcs[fj]) ) 
                    continue;

                auto cofsj = get_system_offset(msh, fcs[fj]);

                for (size_t i = 0; i < fbs; i++)
                {
                    for(size_t j = 0; j < fbs; j++)
                    {
                        auto lofsi = fi*fbs;
                        auto lofsj = fj*fbs;
                        triplets.push_back( Triplet<ScalT>(cofsi+i, cofsj+j, lhsc(lofsi+i, lofsj+j)) );
                    }
                }
            }

            RHS.segment(cofsi, fbs) += rhs.segment(fi*fbs, fbs);
        }
    }

    void finalize(void)
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        std::cout << "Solving for " << syssz << " DOFs. ";
        std::cout << "Matrix has " << LHS.nonZeros() << " nonzeros." << std::endl; 
    }

    disk::dynamic_vector<ScalT>
    get_expanded_solution(const Mesh& msh, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();

        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*msh.faces_size() );

        for (size_t i = 0; i < sysfcs; i++)
        {
            auto in_offset = i*fbs;
            auto out_offset = expand_table.at(i)*fbs;
            ret.segment(out_offset, fbs) = sol.segment(in_offset, fbs);
        }

        return ret;
    }

    disk::dynamic_vector<ScalT>
    get_element_dofs(const Mesh& msh, const typename Mesh::cell& cl, disk::dynamic_vector<ScalT>& sol)
    {
        auto fbs = face_basis_size();
        auto fcs = faces(msh, cl);
        disk::dynamic_vector<ScalT> ret = disk::dynamic_vector<ScalT>::Zero( fbs*fcs.size() );

        for (size_t i = 0; i < fcs.size(); i++)
        {
            auto& fc = fcs[i];
            if ( not is_in_system(msh, fc) )
                continue;

            auto ofs = get_system_offset(msh, fc);
            ret.segment(i*fbs, fbs) = sol.segment(ofs, fbs);
        }

        return ret;
    }
};

template<typename Mesh>
class maxwell_hho_solver
{
};

int main(int argc, const char **argv)
{
    sol::state lua;
    lua.open_libraries(sol::lib::math, sol::lib::base);

}