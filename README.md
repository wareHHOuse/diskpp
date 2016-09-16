### DiSk++

This is DiSk++, a template library for Discontinuous Skeletal methods (discontinuous Galerkin, HHO, ...) written in C++. The philosophy behind the library is "write once, run on any kind of mesh". This means that the generic nature of the library allows to code a numerical method only one time, and run it on any kind of mesh (1D, 2D, 3D, simplicial, hexahedral, ...)

For example, to compute a mass matrix inside an element:

```C++
scaled_monomial_scalar_basis<mesh_type, cell_type> cb(degree);
quadrature<mesh_type, cell_type> cq(2*degree);
for (auto& cl : msh)
{
    auto quadpoints = cq.integrate(msh, cl);
    for (auto& qp : quadpoints)
    {
        auto pt = qp.point();
        auto phi  = cb.eval_functions(msh, cl, pt);

        for (size_t i = 0; i < phi.size(); i++)
            for (size_t j = 0; j < phi.size(); j++)
                mass(i,j) = qp.weight() * mm_prod(phi[i], phi[j]);
    }
}
```

this code will work independently of the dimension of the problem (1D, 2D, 3D) and the kind of mesh.
