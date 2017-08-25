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

    disk::surface<0, ud_a, disk::dynamic_storage_polytope> s1({4,5,6});
    s1.domain_id(42);
    std::cout << s1 << std::endl;
    std::cout << sizeof(s1) << std::endl;

    disk::surface<0, void, disk::dynamic_storage_polytope> s2(s1);
    std::cout << s2 << std::endl;
    std::cout << sizeof(s2) << std::endl;

    return 0;
}

