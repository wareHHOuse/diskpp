### DiSk++

This is __DiSk++__, a C++template library for Discontinuous Skeletal methods like discontinuous Galerkin and Hybrid-High Order.

In mathematical literature, discontinuous skeletal methods are always presented and studied in a dimension-indipendent and element-indipendent fashon. That is, the dimension of the space (1D, 2D, 3D) and the shape of the elements in which the space is discretized are just ininfluent details for the mathematical discussion of these methods. The software implementing them, however, does not always take the same approach: it is common to see codes capable to run only on few very specific kinds of mesh, or only in 1D or 2D or 3D.

The philosophy behind DiSk++ is "write once, run on any kind of mesh". This means that the generic nature of the library allows to code a numerical method only one time, and run it on any kind of mesh (1D, 2D, 3D, simplicial, hexahedral, ...)

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

## Installation
The library runs on Unix systems. The main development is made on Mac OS X, but it compiles and runs fine also on Linux and FreeBSD. It is written in C++14 and requires a recent compiler to be compiled (GCC >= 5.0 or Clang >= 3.8). Older compilers may work but they are neither supported nor tested.

### Mac OS X
If you want to just run the examples, in Mac OS X is very easy. You just need to install [Homebrew](http://brew.sh) and then

    brew tap datafl4sh/code
    brew install --HEAD datafl4sh/code/diskpp

You will end up with an executable named `diffusion` and some meshes in `/usr/local/share/meshes`. The problem the demo program solves is the classical diffusion problem

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bcases%7D%0D%0A%5CDelta+u+%3D+f+%26+%5Ctext%7Bin%5C%3B%5C%3B%7D+%5COmega%5C%5C%0D%0Au+%3D+0+%26+%5Ctext%7Bon%5C%3B%5C%3B%7D+%5Cpartial%5COmega%0D%0A%5Cend%7Bcases%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

with the forcing term

![equation2](http://www.sciweavers.org/tex2img.php?eq=f%3D%5Cpi%5E2+sin%28%5Cpi+x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

The solution is of course

![equation3](http://www.sciweavers.org/tex2img.php?eq=sin%28%5Cpi+x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

#### Running a 1D simulation
To run an 1D simulation it is sufficient to call the program without arguments. If you want, you can change the polynomial degree with `-k` and the number of mesh elements with `-n`

A file named `plot.dat` will be generated, to view the result it is sufficient to launch Gnuplot and at the prompt to say

    plot 'plotnew.dat' using 1:2 with linespoints

#### Running a 2D simulation
To run a 2D simulation it is sufficient to pass as last argument the path of a 2D mesh. In the `meshes` directory you will find some of the FVCA5 benchmark meshes, which have extension `.typ1` and a couple of simplicial meshes generated with Netgen, which have extension `.mesh2d`.

    diffusion -k <degree> /usr/local/share/meshes/<meshname>

A file named `plot.dat` will be generated, to view the result it is sufficient to launch Gnuplot and at the prompt to say

    splot 'plotnew.dat' using 1:2:3 with points palette ps 1 pt 7

#### Running a 3D simulation
The same as the 2D case. The 3D meshes have extension `.mesh`.

    diffusion -k <degree> /usr/local/share/meshes/<meshname>

A file named `plot.dat` will be generated, to view the result it is sufficient to launch Gnuplot and at the prompt to say

    splot 'plotnew.dat' using 1:2:3:4 with points palette ps 1 pt 7
