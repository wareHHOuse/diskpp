### DiSk++

This is __DiSk++__, a C++ template library for Discontinuous Skeletal methods
like discontinuous Galerkin (dG) and Hybrid-High Order (HHO).

In mathematical literature, discontinuous skeletal methods are always presented
and studied in a dimension-indipendent and element-indipendent fashon. That is,
the dimension of the space (1D, 2D, 3D) and the shape of the elements in which
the space is discretized are just ininfluent details for the mathematical
discussion of these methods. The software implementing them, however, does not
always take the same approach: it is common to see codes capable to run only
on few very specific kinds of mesh, or only in 1D or 2D or 3D.

The philosophy behind DiSk++ is "write once, run on any kind of mesh". This
means that the generic nature of the library allows to code a numerical method
only one time, and run it on any kind of mesh (1D, 2D, 3D, simplicial,
hexahedral, ...)

## Cloning the repository
The repository has some dependencies that must be fetched. When you clone,
make sure to do a recursive clone by using

```
git clone --recursive https://github.com/wareHHOuse/diskpp
```

## Installation
The library runs on Unix systems only. Windows is explicitly NOT supported.
The development is focused on Linux and Mac OS X in particular. DiSk++ is
written in C++17, and C++20 is coming. Be sure to use a recent compiler.

### Linux
DiSk++ requires several packages to be installed in order to compile correctly.

```
apt-get install -y make cmake build-essential git lua5.4 liblua5.4-dev libmumps-seq-5.3 libmumps-seq-dev gmsh libgmsh-dev libsilo-dev libsiloh5-0 libeigen3-dev
```

