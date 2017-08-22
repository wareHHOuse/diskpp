#include "element_polytopes.hpp"

struct ud_a
{
    int x;
    int y;
};

struct ud_b
{
    double x;
    int y;
};

int main(void)
{

    disk::surface<0, ud_a, disk::limited_storage_polytope<5>> s1({4,5,6});
    s1.domain_id(42);
    std::cout << s1 << std::endl;
    std::cout << sizeof(s1) << std::endl;

    disk::surface<0, void, disk::limited_storage_polytope<5>> s2(s1);
    std::cout << s2 << std::endl;
    std::cout << sizeof(s2) << std::endl;

    disk::surface<0, ud_b, disk::limited_storage_polytope<5>> s3(s2);
    std::cout << s3 << std::endl;
    std::cout << sizeof(s3) << std::endl;

    disk::surface<1, void, disk::limited_storage_polytope<5>> s4({4,5,6});
    s4.boundary_id(42);
    std::cout << s4 << std::endl;
    std::cout << sizeof(s4) << std::endl;

    disk::surface<2, void, disk::limited_storage_polytope<3>> s5({4,5,6});
    std::cout << s5 << std::endl;
    std::cout << sizeof(s5) << std::endl;

    disk::node<1, void> n({4});
    std::cout << n << std::endl;

    return 0;
}

