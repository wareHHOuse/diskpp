### DiSk++

This is __DiSk++__, a C++ template library for Discontinuous Skeletal methods like discontinuous Galerkin (dG) and Hybrid-High Order (HHO).

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

        mass += qp.weight() * phi * phi.transpose();
    }
}
```

this code will work independently of the dimension of the problem (1D, 2D, 3D) and the kind of mesh.

## Cloning the repository
The repository has some dependencies that must be fetched. When you clone, make sure to do a recursive clone by using

```git clone --recursive https://github.com/datafl4sh/diskpp```

## Installation
The library runs on Unix systems. The main development is made on Mac OS X, but it compiles and runs fine also on Linux and FreeBSD. It is written in C++14 and requires a recent compiler to be compiled (GCC >= 5.0 or Clang >= 3.8). Older compilers may work but they are neither supported nor tested.

### Mac OS X
If you want to just run the examples, in Mac OS X is very easy. You just need to install [Homebrew](http://brew.sh) and then

    brew tap datafl4sh/code
    brew install --HEAD datafl4sh/code/diskpp

You will end up with an executable named `diffusion` and some meshes in `/usr/local/share/meshes`. The problem the demo program solves is a classical diffusion problem.

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

### Linux
The project requires several packages to be installed on the system:

1. Intel Math Kernel Library (MKL): https://software.intel.com/en-us/mkl
diskpp-master will be able to find MKL after a correct installation. Check cmake/FindMKL.cmake if an error occures.

2. Silo library for reading and writing a wide variety of scientific data to binary, disk files: https://wci.llnl.gov/simulation/computer-codes/silo
If diskpp-master can not find Silo after installation, add LIBRARY and INCLUDE path to cmake/FindSILO.cmake.

3. Lua lightweight embeddable scripting language: https://www.lua.org/download.html
Lua is also available on most Linux platforms. Version 5.2 has been tested, make sure to install both main and development files, e.g.:

```
    sudo apt-get install lua5.2
    sudo apt-get install liblua5.2-dev
```

4. Eigen C++ template library for linear algebra: http://eigen.tuxfamily.org/