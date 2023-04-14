#pragma once

#include "diskpp/mesh/point.hpp"

namespace disk {

/* pqrst is there only to squelch -fpermissive
 * GCC warnings. sorry for this. */

template<typename T, size_t DIM>
using pqrst = point<T,DIM>;

template<typename T, size_t DIM>
class quadrature_point
{
    pqrst<T, DIM>   quad_point;
    T               quad_weight;

public:
    quadrature_point()
    {}

    quadrature_point(const pqrst<T, DIM>& qp, const T& qw)
        : quad_point(qp), quad_weight(qw)
    {}

    pqrst<T, DIM> point() const {
        return quad_point;
    }

    T weight() const {
        return quad_weight;
    }
};

template<typename T, size_t DIM>
quadrature_point<T, DIM>
make_qp(const point<T, DIM>& point, const T& weight)
{
    return quadrature_point<T, DIM>(point, weight);
}

//
template<typename T, size_t DIM>
std::ostream&
operator<<(std::ostream& os, const quadrature_point<T,DIM>& qp)
{
    os << qp.point() << " " << qp.weight();
    return os;
}

} // namespace disk