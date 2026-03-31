#include <iostream>
#include <vector>
#include "diskpp/common/eigen.hpp"

int main(void)
{
    using T = double;
    disk::sparse_matrix<T> Alhs;
    std::cout << "loading matrix\n"; 
    disk::load_from_hdf5(Alhs, "/home/matteo/tmp/Alhs.h5");

    std::cout << "loading vector\n";
    std::vector<T> rhs_v;
    H5Easy::File file("/home/matteo/tmp/RHS.h5", H5Easy::File::ReadOnly);
    rhs_v = file.getDataSet("/densematrix").read<std::vector<T>>();
    disk::dynamic_vector<T> RHS = disk::dynamic_vector<T>::Zero(rhs_v.size());
    for (size_t i = 0; i < rhs_v.size(); i++) {
        RHS(i) = rhs_v[i];
    }

    Eigen::PardisoLU<Eigen::SparseMatrix<T>> lu(Alhs);
    disk::dynamic_vector<T> sol = lu.solve(RHS);

    std::cout << sol << std::endl;

    return 0;
}