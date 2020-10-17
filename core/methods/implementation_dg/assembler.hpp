/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group  
 */

#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace disk {

template<typename Mesh>
class discontinuous_galerkin_assembler
{
    size_t cbs;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using scal_type = typename mesh_type::coordinate_type;
    using matrix_type = Eigen::Matrix<scal_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<scal_type, Eigen::Dynamic, 1>;
    using triplet_type = Triplet<scal_type>;
    using T = scal_type;

    std::vector<triplet_type> triplets;

public:
    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    size_t                  syssz;

    discontinuous_galerkin_assembler(const Mesh& msh, size_t pcbs)
    {
        cbs = pcbs;
        syssz = cbs * msh.cells_size();

        LHS = SparseMatrix<T>(syssz, syssz);
        RHS = Matrix<T, Dynamic, 1>::Zero(syssz);
    }

    void
    assemble(const mesh_type& msh, const cell_type& clp, const cell_type& clm,
             const matrix_type& A)
    {
        auto ofs_i = cbs * offset(msh, clp);
        auto ofs_j = cbs * offset(msh, clm);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets.push_back( triplet_type(ofs_i+i, ofs_j+j, A(i,j)) );
    }

    void
    assemble(const mesh_type& msh, const cell_type& cl, const matrix_type& A,
             const vector_type& b)
    {
        auto ofs = cbs * offset(msh, cl);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets.push_back( triplet_type(ofs+i, ofs+j, A(i,j)) );

        RHS.segment(ofs, cbs) += b;
    }

    void
    finalize()
    {
        LHS.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
    }
};

template<typename Mesh>
auto
make_discontinuous_galerkin_assembler(const Mesh& msh, size_t cbs)
{
    return discontinuous_galerkin_assembler<Mesh>(msh, cbs);
}

template<typename Mesh>
class discontinuous_galerkin_eigenvalue_assembler
{
    size_t cbs;

    using mesh_type = Mesh;
    using cell_type = typename mesh_type::cell_type;
    using scal_type = typename mesh_type::coordinate_type;
    using matrix_type = Eigen::Matrix<scal_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type = Eigen::Matrix<scal_type, Eigen::Dynamic, 1>;
    using triplet_type = Triplet<scal_type>;
    using T = scal_type;

    std::vector<triplet_type> triplets_gK;
    std::vector<triplet_type> triplets_gM;

public:
    SparseMatrix<T>         gK;
    SparseMatrix<T>         gM;

    size_t                  syssz;

    discontinuous_galerkin_eigenvalue_assembler(const Mesh& msh, size_t pcbs)
    {
        cbs = pcbs;
        syssz = cbs * msh.cells_size();

        gK = SparseMatrix<T>(syssz, syssz);
        gM = SparseMatrix<T>(syssz, syssz);
    }

    void
    assemble(const mesh_type& msh, const cell_type& clp, const cell_type& clm,
             const matrix_type& A)
    {
        auto ofs_i = cbs * offset(msh, clp);
        auto ofs_j = cbs * offset(msh, clm);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets_gK.push_back( triplet_type(ofs_i+i, ofs_j+j, A(i,j)) );
    }

    void
    assemble(const mesh_type& msh, const cell_type& cl, const matrix_type& A)
    {
        auto ofs = cbs * offset(msh, cl);

        for (size_t i = 0; i < cbs; i++)
            for (size_t j = 0; j < cbs; j++)
                triplets_gM.push_back( triplet_type(ofs+i, ofs+j, A(i,j)) );
    }

    void
    finalize()
    {
        gK.setFromTriplets(triplets_gK.begin(), triplets_gK.end());
        triplets_gK.clear();
        gM.setFromTriplets(triplets_gM.begin(), triplets_gM.end());
        triplets_gM.clear();
    }
};

template<typename Mesh>
auto
make_discontinuous_galerkin_eigenvalue_assembler(const Mesh& msh, size_t cbs)
{
    return discontinuous_galerkin_eigenvalue_assembler<Mesh>(msh, cbs);
}

}
