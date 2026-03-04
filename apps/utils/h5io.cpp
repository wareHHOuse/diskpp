#include <iostream>

#include "diskpp/common/eigen.hpp"

int main(void)
{
    disk::sparse_matrix<double> Ms;
    Ms.resize(10,10);
    Ms.setIdentity();

    disk::save_to_hdf5(Ms, "sparse.h5");

    disk::sparse_matrix<double> Ns;
    disk::load_from_hdf5(Ns, "sparse.h5");

    std::cout << Ns << std::endl;

    disk::dynamic_matrix<double> Md = disk::dynamic_matrix<double>::Ones(4,5);
    disk::save_to_hdf5(Md, "dense.h5");

    disk::dynamic_matrix<double> Nd;
    disk::load_from_hdf5(Nd, "dense.h5");

    std::cout << Nd << std::endl;
    return 0;
}
