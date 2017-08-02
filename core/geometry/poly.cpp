#include "element_polytopes.hpp"

struct ud
{
    int x;
    int y;
};

int main(void)
{
    disk::simplicial_element<2> poly({1,2,3});

    std::cout << poly << std::endl;
    //std::cout << poly.user_data.x << std::endl;

    std::cout << sizeof(poly) << std::endl;

    return 0;
}