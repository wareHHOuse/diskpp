#pragma once

namespace disk::basis {

namespace priv {
template<typename Basis>
struct dot_evaluator
{
    const Basis& basis;
    double n;

    using point_type = typename Basis::point_type;

    dot_evaluator() = delete;
    dot_evaluator(const dot_evaluator&) = delete;

    dot_evaluator(const Basis& pb, double pn)
        : basis(pb), n(pn)
    {}

    auto operator()(const point_type& pt) const {
        return basis(pt)*n;
    }

    size_t size() const {
        return basis.size();
    }

    size_t degree() const {
        return basis.degree();
    }

    size_t integration_degree() const {
        return basis.degree();
    }
};

template<typename Basis>
struct grad_evaluator
{
    const Basis& basis;

    using mesh_type = typename Basis::mesh_type;
    using point_type = typename Basis::point_type;

    grad_evaluator() = delete;
    grad_evaluator(const grad_evaluator&) = delete;

    grad_evaluator(const Basis& pb)
        : basis(pb)
    {}

    auto operator()(const point_type& pt) const {
        return basis.grad(pt);
    }

    auto dot(const double& n) {
        return dot_evaluator(*this, n);
    }

    size_t size() const {
        return basis.size();
    }

    size_t degree() const {
        return basis.degree();
    }

    size_t integration_degree() const {
        if (basis.degree() > 0)
            return basis.degree()-1;

        return basis.degree();
    }
};

} //namespace priv

template<typename Basis>
auto grad(const Basis& basis) {
    return priv::grad_evaluator(basis);
}

} //namespace disk::basis