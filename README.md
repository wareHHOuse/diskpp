### DiSk++

This is __DiSk++__, a C++ template library for Discontinuous Skeletal methods like discontinuous Galerkin (dG) and Hybrid-High Order (HHO).

In mathematical literature, discontinuous skeletal methods are always presented and studied in a dimension-indipendent and element-indipendent fashon. That is, the dimension of the space (1D, 2D, 3D) and the shape of the elements in which the space is discretized are just ininfluent details for the mathematical discussion of these methods. The software implementing them, however, does not always take the same approach: it is common to see codes capable to run only on few very specific kinds of mesh, or only in 1D or 2D or 3D.

The philosophy behind DiSk++ is "write once, run on any kind of mesh". This means that the generic nature of the library allows to code a numerical method only one time, and run it on any kind of mesh (1D, 2D, 3D, simplicial, hexahedral, ...)

## Cloning the repository
The repository has some dependencies that must be fetched. When you clone, make sure to do a recursive clone by using

```git clone --recursive https://github.com/datafl4sh/diskpp```

After cloning, be sure to checkout the branch `devel` if you want all the latest stuff. The branch master usually lags much behind `devel`.

## Installation
The library runs on Unix systems and the development is focused on Linux and Mac OS X in particular. DiSk++ is written mainly in C++14, however we are using C++17 and c++20 more and more. Be sure to use a recent compiler.

### Linux
The project requires several packages to be installed on the system:

1. Intel Math Kernel Library (MKL): https://software.intel.com/en-us/mkl
The build system expects MKL to be installed in `/opt/intel`, check cmake/FindMKL.cmake if an error occours.

2. Silo library for reading and writing a wide variety of scientific data to binary files: https://wci.llnl.gov/simulation/computer-codes/silo
Install with
```
    sudo apt-get install libsilo-dev libsilo-bin
```
On Linux recent versions of Silo complain that `silo_exports.h` is not present. Please grab it from the LLNL Silo package and copy it in `/usr/include`. It is not clean, but it is the current workaround.

3. Lua lightweight embeddable scripting language: https://www.lua.org/download.html
Lua is also available on most Linux platforms. Version 5.2 has been tested, make sure to install both main and development files, e.g.:

```
    sudo apt-get install lua5.2
    sudo apt-get install liblua5.2-dev
```

4. Eigen C++ template library for linear algebra: http://eigen.tuxfamily.org/

5. HDF5 development files.
```
    apt-get install libhdf5-dev
```
