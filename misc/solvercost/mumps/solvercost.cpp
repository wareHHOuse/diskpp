/*
 * Copyright (C) 2017, Matteo Cicuttin - matteo.cicuttin@enpc.fr
 * École Nationale des Ponts et Chaussées, Université Paris EST
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <mpi.h>

#include "mumps.hpp"

using scalar_type   = double;
template<typename T>
using triplet_type  = Eigen::Triplet<T>;
using spmat_type    = Eigen::SparseMatrix<scalar_type>;
using vector_type   = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

template<typename T>
triplet_type<T>
read_matrix_file_entry(std::ifstream& ifs)
{
    int row, col;
    T   value;

    ifs >> row >> col >> value;

    return triplet_type<T>(row, col, value);
}

template<typename T>
T
read_rhs_file_entry(std::ifstream& ifs)
{
    T   value;
    ifs >> value;
    return value;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " <matrix file> <rhs file>";
        std::cout << std::endl;
        return 1;
    }

    /* Read matrix */
    std::ifstream ifsm(argv[1]);
    if ( !ifsm.is_open() )
    {
        std::cout << "Unable to open " << argv[1] << std::endl;
        return 1;
    }

    std::vector<triplet_type<scalar_type>> triplets;
    size_t matrix_size, matrix_nnz;

    ifsm >> matrix_nnz >> matrix_size;
    triplets.reserve(matrix_nnz);
    for (size_t i = 0; i < matrix_nnz; i++)
        triplets.push_back( read_matrix_file_entry<scalar_type>(ifsm) );

    ifsm.close();

    /* Read rhs */
    std::ifstream ifsr(argv[2]);
    if ( !ifsr.is_open() )
    {
        std::cout << "Unable to open " << argv[2] << std::endl;
        return 1;
    }

    vector_type rhs;
    size_t rhs_size;
    ifsr >> rhs_size;
    rhs.resize(rhs_size);

    for (size_t i = 0; i < rhs_size; i++)
        rhs(i) = read_rhs_file_entry<scalar_type>(ifsr);

    ifsr.close();

    spmat_type  matrix(matrix_size, matrix_size);
    matrix.setFromTriplets(triplets.begin(), triplets.end());

    mumps_solver<scalar_type> mumps;

    vector_type sol = mumps.solve(matrix, rhs);

    MPI_Finalize();

    return 0;
}
