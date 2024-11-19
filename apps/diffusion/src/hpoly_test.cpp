#include <iostream>
#include <complex>

#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/bases/bases_new.hpp"

int main(void)
{
    using T = double;
    using mesh_type = disk::generic_mesh<T, 2>;
    
    auto num_faces = 4;

    mesh_type msh;
    using point_type = typename mesh_type::point_type;

    disk::make_single_element_mesh(msh, 1.0, num_faces);

    auto tr = [](const point_type& pt) {
        return point_type(pt.x()+5.0, pt.y()+4.0);
    };

    msh.transform(tr);

    auto poly = [&](const point_type& pt) {
        auto n = num_faces-1;
        auto rho = std::hypot(pt.x(), pt.y());
        auto theta = std::atan2(pt.y(), pt.x());
        auto z = std::pow(rho,n)*std::exp( std::complex<T>(0,n*theta) );
        return real(z);
    };

    auto cl = *msh.cells_begin();

    double sint_poly = disk::basis::integrate(msh, cl, num_faces-1, poly);
    double lint_poly = 0.0;

    auto area = measure(msh, cl);
    auto perimeter = 0.0;
    auto fcs = faces(msh, cl);
    for (auto& fc : fcs) {
        perimeter += measure(msh, fc);
        lint_poly += disk::basis::integrate(msh, fc, num_faces-1, poly);
    }

    std::cout << "Area: " << area << ", Perimeter: " << perimeter << std::endl;
    std::cout << "Surface integral: " << sint_poly << ", line integral: " << lint_poly << std::endl;
    std::cout << sint_poly/area - lint_poly/perimeter << std::endl;

    return 0;
}