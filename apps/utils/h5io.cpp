#include <iostream>

#include "diskpp/common/eigen.hpp"

int main(void)
{
    disk::sparse_matrix<double> M;
    M.resize(10,10);
    M.setIdentity();

    disk::save_to_hdf5(M, "matrix.h5");

    return 0;
}
