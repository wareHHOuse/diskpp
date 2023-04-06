#pragma once

template<typename T>
struct computation_info
{
    T           l2_error_e;
    T           l2_error_h;
    T           nrg_error;
    int         mflops;
    size_t      dofs;
    long        nonzeros;
};
