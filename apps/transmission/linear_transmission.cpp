/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *  
 * Matteo Cicuttin (C) 2020,2021
 * matteo.cicuttin@uliege.be
 *
 * University of Liège - Montefiore Institute
 * Applied and Computational Electromagnetics group
 */

/* Implementation of "A mixed Hybrid High-Order formulation for
 * linear interior transmission elliptic problems"
 * by R. Bustinza & J. Munguia La-Cotera
 */

#include <iostream>
#include <iomanip>
#include <regex>

#include <sstream>
#include <string>

#include <unistd.h>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"

#include "sol/sol.hpp"
#include "mumps.hpp"

#include "sgr.hpp"

#define PARDISO_IN_CORE                 0
#define PARDISO_OUT_OF_CORE_IF_NEEDED   1
#define PARDISO_OUT_OF_CORE_ALWAYS      2

template<typename T>
struct pardiso_params
{
    bool    report_factorization_Mflops;
    int     out_of_core; //0: IC, 1: IC, OOC if limits passed, 2: OOC
    int     mflops;

    pardiso_params() :
        report_factorization_Mflops(false), out_of_core(0), mflops(0)
    {}
};

template<typename T>
bool
mkl_pardiso(pardiso_params<T>& params,
            const Eigen::SparseMatrix<T>& A,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
            Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    Eigen::PardisoLU<Eigen::SparseMatrix<T>>  solver;

    if (params.out_of_core >= 0 && params.out_of_core <= 2)
        solver.pardisoParameterArray()[59] = params.out_of_core;

    if (params.report_factorization_Mflops)
        solver.pardisoParameterArray()[18] = -1; //report flops

    solver.analyzePattern(A);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: analyzePattern failed" << std::endl;
       return false;
    }

    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: Could not factorize the matrix" << std::endl;
       std::cerr << "Try to tweak MKL_PARDISO_OOC_MAX_CORE_SIZE" << std::endl;
       return false;
    }

    x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
       std::cerr << "ERROR: Could not solve the linear system" << std::endl;
       return false;
    }

    if (params.report_factorization_Mflops)
    {
        int mflops = solver.pardisoParameterArray()[18];
        params.mflops = mflops;
        std::cout << "[PARDISO] Factorization Mflops: " << mflops << std::endl;
    }

    return true;
}

template<typename Mesh>
std::set<size_t>
subdomain_tags(const Mesh& msh)
{
    std::set<size_t> tags;

    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        tags.insert( di.tag() );
    }

    return tags;
}


struct element_counts
{
    size_t cells_int;
    size_t faces_int;
    size_t cells_ext;
    size_t faces_ext;
    size_t faces_iface;

    element_counts()
        : cells_int(0), faces_int(0), cells_ext(0), faces_ext(0), faces_iface(0)
    {}
};

std::ostream& operator<<(std::ostream& os, const element_counts& ec)
{
    os << "IC: " << ec.cells_int << ", IF: " << ec.faces_int << ", EC: " << ec.cells_ext;
    os << ", EF: " << ec.faces_ext << ", JF: " << ec.faces_iface;
    return os;
}

template<typename Mesh>
element_counts count_mesh_elements(const Mesh& msh, size_t internal_tag)
{
    element_counts ret;

    std::set<size_t> ifaces;
    std::set<size_t> efaces;

    for (auto& cl : msh)
    {
        auto fcs = faces(msh, cl);
        auto di = msh.domain_info(cl);

        if (di.tag() == internal_tag)
            ret.cells_int++;
        else
            ret.cells_ext++;
        for (auto& fc : fcs)
        {
            if (di.tag() == internal_tag)
                ifaces.insert( offset(msh, fc) );
            else
                efaces.insert( offset(msh, fc) );

            auto bi = msh.boundary_info(fc);
            if (bi.is_internal() and di.tag() == internal_tag)
                ret.faces_iface++;
        }
    }

    ret.faces_int = ifaces.size();
    ret.faces_ext = efaces.size();

    return ret;
}

struct dof_bases
{
    size_t cells_int_base;
    size_t cells_ext_base;
    size_t global_constrain_base;
    size_t faces_int_base;
    size_t faces_ext_base;
    size_t multiplier_base;
    size_t total_dofs;
};

template<typename Mesh>
dof_bases compute_dof_bases(const Mesh& msh, const disk::hho_degree_info& hdi, size_t internal_tag)
{
    auto ec  = count_mesh_elements(msh, internal_tag);
    auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

    dof_bases ret;

    ret.cells_int_base          = 0;
    ret.cells_ext_base          = ret.cells_int_base + cbs * ec.cells_int;
    ret.global_constrain_base   = ret.cells_ext_base + cbs * ec.cells_ext;
    ret.faces_int_base          = ret.global_constrain_base + 1;
    ret.faces_ext_base          = ret.faces_int_base + fbs * ec.faces_int;
    ret.multiplier_base         = ret.faces_ext_base + fbs * ec.faces_ext;
    ret.total_dofs              = ret.multiplier_base + fbs * ec.faces_iface;

    //std::cout << ret.cells_int_base << " " << ret.cells_ext_base << " ";
    //std::cout << ret.global_constrain_base << " " << ret.faces_int_base << " ";
    //std::cout << ret.faces_ext_base << " " << ret.multiplier_base << " ";
    //std::cout << ret.total_dofs << std::endl;
    //std::cout << ec << std::endl;

    return ret;
}
struct offsets_g2s
{
    std::map<size_t, size_t>    cells_int_g2s;
    std::map<size_t, size_t>    faces_int_g2s;
    std::map<size_t, size_t>    cells_ext_g2s;
    std::map<size_t, size_t>    faces_ext_g2s;
    std::map<size_t, size_t>    interface_g2s;
};

/* Compute the subdomain-local element numberings */
template<typename Mesh>
void compute_g2s(const Mesh& msh, size_t internal_tag, offsets_g2s& g2s)
{
    size_t ci = 0;
    size_t ce = 0;
    size_t fi = 0;
    size_t fe = 0;
    size_t ii = 0;

    auto count = [](size_t offset, std::map<size_t, size_t>& g2s_l, size_t& cnt) {
        if ( g2s_l.find(offset) != g2s_l.end() )
            return;

        g2s_l[offset] = cnt++;
    };

    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        if (di.tag() == internal_tag)
            g2s.cells_int_g2s[ offset(msh, cl) ] = ci++;
        else
            g2s.cells_ext_g2s[ offset(msh, cl) ] = ce++;

        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            auto fcofs = offset(msh, fc);
            if (di.tag() == internal_tag)
                count( fcofs, g2s.faces_int_g2s, fi);
            else
                count( fcofs, g2s.faces_ext_g2s, fe);

            auto bi = msh.boundary_info(fc);
            if ( bi.is_internal() )
                count( fcofs, g2s.interface_g2s, ii);
        }
    }
}

/***************************************************************************/
template<typename Mesh>
struct problem_data_base
{
    using mesh_type     = Mesh;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, Mesh::dimension, 1>;

    scalar_type     avgcomp1_val;
    scalar_type     avgcomp2_val;

    problem_data_base()
        : avgcomp1_val(0.0), avgcomp2_val(0.0)
    {}

    virtual scalar_type u1(const mesh_type&, const point_type&) = 0;
    virtual scalar_type u2(const mesh_type&, const point_type&) = 0;

    virtual vector_type grad_u1(const mesh_type&, const point_type&) = 0;
    virtual vector_type grad_u2(const mesh_type&, const point_type&) = 0;

    virtual scalar_type f1(const mesh_type&, const point_type&) = 0;
    virtual scalar_type f2(const mesh_type&, const point_type&) = 0;

    scalar_type avgcomp1(void) const {
        return avgcomp1_val;
    }

    void avgcomp1(scalar_type val) {
        avgcomp1_val = val;
    }

    scalar_type avgcomp2(void) const {
        return avgcomp2_val;
    }

    void avgcomp2(scalar_type val) {
        avgcomp2_val = val;
    }

    scalar_type g(const mesh_type& msh, const point_type& pt) {
        return u1(msh, pt) - u2(msh, pt);
    }

    scalar_type g1(const mesh_type& msh, const cell_type& cl, const face_type& fc, const point_type& pt) {
        auto n = normal(msh, cl, fc);
        return grad_u1(msh, pt).dot(n) - grad_u2(msh, pt).dot(n);
    }

    scalar_type g2(const mesh_type& msh, const cell_type& cl, const face_type& fc, const point_type& pt) {
        auto n = normal(msh, cl, fc);
        return grad_u2(msh, pt).dot(n);
    }

    scalar_type xi(const mesh_type& msh, const cell_type& cl, const face_type& fc, const point_type& pt) {
        auto n = normal(msh, cl, fc);
        return grad_u1(msh, pt).dot(n);
    }

    virtual std::string solver() {
        return "pardiso";
    }

    virtual int omega1(void) {
        return 1;
    }

    virtual int omega2(void) {
        return 2;
    }
};

/***************************************************************************/
template<typename T> struct problem_data;

#if 0
template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct problem_data<Mesh<T,2,Storage>> : problem_data_base<Mesh<T,2,Storage>>
{
    using mesh_type     = Mesh<T,2,Storage>;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, 2, 1>;

    scalar_type u1(const mesh_type&, const point_type& pt) {
        return std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y());
    }

    scalar_type u2(const mesh_type& msh, const point_type& pt) {
        return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    }

    vector_type grad_u1(const mesh_type&, const point_type& pt) {
        vector_type ret;
        ret (0) = -M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y());
        ret (1) = -M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    }

    vector_type grad_u2(const mesh_type& msh, const point_type& pt) {
        vector_type ret;
        ret (0) = M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y());
        ret (1) = M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y());
        return ret;
    }

    scalar_type f1(const mesh_type&, const point_type& pt) {
        return 2*M_PI*M_PI*std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y());
    }

    scalar_type f2(const mesh_type& msh, const point_type& pt) {
        return 2*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    }
};
#endif



template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct problem_data<Mesh<T,2,Storage>> : problem_data_base<Mesh<T,2,Storage>>
{
    using mesh_type     = Mesh<T,2,Storage>;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, 2, 1>;

    scalar_type u1(const mesh_type&, const point_type& pt) {
        double aa = (13*std::log(2)-5*std::log(5))/14;
        double aux = std::pow(pt.x(),2)+std::pow(pt.y(),2);
        return pt.x()*pt.y()*std::pow(aux, -1)-aa;
    }

    scalar_type u2(const mesh_type& msh, const point_type& pt) {
        std::complex<double> mycomplex (pt.y()-pt.x(),-pt.x()-pt.y()); 
        double a = std::arg(mycomplex)+0.75*M_PI;
        double r = std::pow(pt.x(), 2)+std::pow(pt.y(), 2);
        return std::pow(r,1./3)*sin((2./3)* a) - 1.08782791064557;
    }

    vector_type grad_u1(const mesh_type&, const point_type& pt) {
        double aux = std::pow(pt.x(),2)+std::pow(pt.y(),2);
        
        vector_type ret;
        ret(0) = pt.y()*(std::pow(pt.y(),2)-std::pow(pt.x(),2))*std::pow(aux,-2),
        ret(1) = pt.x()*(std::pow(pt.x(),2)-std::pow(pt.y(),2))*std::pow(aux,-2);
        return ret;
    }

    vector_type grad_u2(const mesh_type& msh, const point_type& pt) {
        std::complex<double> mycomplex (pt.y()-pt.x(),-pt.x()-pt.y()); 
        double a = std::arg(mycomplex)+0.75*M_PI;
        double r = std::pow(pt.x(), 2)+std::pow(pt.y(), 2);
        vector_type ret;
                     
        ret(0) = -(2./3)*std::pow(r,-1./6) * sin(a/3),
        ret(1) = (2./3)*std::pow(r,-1./6) * cos(a/3);          
        return ret;
    }

    scalar_type f1(const mesh_type&, const point_type& pt) {
        double aux = std::pow(pt.x(),2)+std::pow(pt.y(),2);          
        return 4*pt.x()*pt.y()*std::pow(aux, -2);
    }

    scalar_type f2(const mesh_type& msh, const point_type& pt) {
        return 0.0;
    }
};




template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct problem_data<Mesh<T,3,Storage>> : problem_data_base<Mesh<T,3,Storage>>
{
    using mesh_type     = Mesh<T,3,Storage>;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, 3, 1>;

    scalar_type u1(const mesh_type&, const point_type& pt) {
        return std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::cos(M_PI*pt.z());
    }

    scalar_type u2(const mesh_type& msh, const point_type& pt) {
        return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*(pt.z()-1.5);
    }

    vector_type grad_u1(const mesh_type&, const point_type& pt) {
        vector_type ret;
        ret (0) = -M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::cos(M_PI*pt.z());
        ret (1) = -M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::cos(M_PI*pt.z());
        ret (2) = -M_PI*std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::sin(M_PI*pt.z());
        return ret;
    }

    vector_type grad_u2(const mesh_type& msh, const point_type& pt) {
        vector_type ret;
        ret (0) = M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y())*(pt.z()-1.5);
        ret (1) = M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y())*(pt.z()-1.5);
        ret (2) = std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        return ret;
    }

    scalar_type f1(const mesh_type&, const point_type& pt) {
        return 3*M_PI*M_PI*std::cos(M_PI*pt.x())*std::cos(M_PI*pt.y())*std::cos(M_PI*pt.z());
    }

    scalar_type f2(const mesh_type& msh, const point_type& pt) {
        return 2*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*(pt.z()-1.5);
    }
};

/***************************************************************************/
template<typename Mesh>
struct problem_data_lua_base : public problem_data_base<Mesh>
{
    sol::state lua;

    bool m_state_valid;

    problem_data_lua_base()
        : m_state_valid(false)
    {}

    void
    open_problem_definitions(const std::string& script_filename)
    {
        lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::io);
        auto script = lua.script_file(script_filename);
        if (not script.valid())
        {
            std::cout << "Script not found. Expected filename is '";
            std::cout << script_filename << "'" << std::endl;
            throw std::runtime_error("");
        }

        if (not lua["u1"].valid())
            throw std::runtime_error("Function 'u1' not defined.");

        if (not lua["u2"].valid())
            throw std::runtime_error("Function 'u2' not defined.");

        if (not lua["grad_u1"].valid())
            throw std::runtime_error("Function 'grad_u1' not defined.");

        if (not lua["grad_u2"].valid())
            throw std::runtime_error("Function 'grad_u2' not defined.");

        if (not lua["f1"].valid())
            throw std::runtime_error("Function 'f1' not defined.");

        if (not lua["f2"].valid())
            throw std::runtime_error("Function 'f2' not defined.");
        
        if (not lua["testcase_name"].valid())
            throw std::runtime_error("Test case name not defined");

        std::string tn = lua["testcase_name"];
        std::cout << "Test case: " << sgr::yellowfg << tn << sgr::nofg << std::endl;

        m_state_valid = true;
    }

    void ensure_state_valid() {
        if (!m_state_valid)
            throw std::runtime_error("Problem definition not valid");
    }

    std::string solver(){
        ensure_state_valid();
        auto sname = lua["solver"];
        if (sname.valid())
            return sname;
        return std::string("internal");
    }

    int omega1(void) {
        ensure_state_valid();
        return lua["omega1"].get_or(1);
    }

    int omega2(void) {
        ensure_state_valid();
        return lua["omega2"].get_or(2);
    }

    bool compensate_average1(void) {
        ensure_state_valid();
        return lua["avgcomp1"].get_or(false);
    }

    bool compensate_average2(void) {
        ensure_state_valid();
        return lua["avgcomp2"].get_or(false);
    }
};

/***************************************************************************/
template<typename T> struct problem_data_lua;

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct problem_data_lua<Mesh<T,2,Storage>> : public problem_data_lua_base<Mesh<T,2,Storage>>
{
    using mesh_type     = Mesh<T,2,Storage>;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, 2, 1>;

    problem_data_lua() {}

    scalar_type u1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type u1_val = this->lua["u1"](pt.x(), pt.y());
        return u1_val - this->avgcomp1_val;
    }

    scalar_type u2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type u2_val = this->lua["u2"](pt.x(), pt.y());
        return u2_val - this->avgcomp2_val;
    }

    vector_type grad_u1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type vx, vy;
        sol::tie(vx, vy) = this->lua["grad_u1"](pt.x(), pt.y());
        return vector_type(vx, vy);
    }

    vector_type grad_u2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type vx, vy;
        sol::tie(vx, vy) = this->lua["grad_u2"](pt.x(), pt.y());
        return vector_type(vx, vy);
    }

    scalar_type f1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        return this->lua["f1"](pt.x(), pt.y());
    }

    scalar_type f2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        return this->lua["f2"](pt.x(), pt.y());
    }
};

template<template<typename, size_t, typename> class Mesh, typename T, typename Storage>
struct problem_data_lua<Mesh<T,3,Storage>> : public problem_data_lua_base<Mesh<T,3,Storage>>
{
    using mesh_type     = Mesh<T,3,Storage>;
    using cell_type     = typename mesh_type::cell_type;
    using face_type     = typename mesh_type::face_type;
    using point_type    = typename mesh_type::point_type;
    using scalar_type   = typename mesh_type::coordinate_type;
    using vector_type   = Eigen::Matrix<scalar_type, 3, 1>;

    problem_data_lua() {}

    scalar_type u1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type u1_val = this->lua["u1"](pt.x(), pt.y(), pt.z());
        return u1_val - this->avgcomp1_val;
    }

    scalar_type u2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type u2_val = this->lua["u2"](pt.x(), pt.y(), pt.z());
        return u2_val - this->avgcomp2_val;
    }

    vector_type grad_u1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type vx, vy, vz;
        sol::tie(vx, vy, vz) = this->lua["grad_u1"](pt.x(), pt.y(), pt.z());
        return vector_type(vx, vy, vz);
    }

    vector_type grad_u2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        scalar_type vx, vy, vz;
        sol::tie(vx, vy, vz) = this->lua["grad_u2"](pt.x(), pt.y(), pt.z());
        return vector_type(vx, vy, vz);
    }

    scalar_type f1(const mesh_type&, const point_type& pt) {
        this->ensure_state_valid();
        return this->lua["f1"](pt.x(), pt.y(), pt.z());
    }

    scalar_type f2(const mesh_type& msh, const point_type& pt) {
        this->ensure_state_valid();
        return this->lua["f2"](pt.x(), pt.y(), pt.z());
    }
};

/***************************************************************************/
template<typename Mesh>
void lt_solver(Mesh& msh, size_t degree, const std::string& pbdefs_fn)
{
    using T = typename Mesh::coordinate_type;

    auto subdom_tags = subdomain_tags(msh);
    if (subdom_tags.size() != 2)
    {
        std::cout << "The mesh must have exactly two subdomains" << std::endl;
        return;
    }

    struct problem_data_lua<Mesh> pd;
    pd.open_problem_definitions(pbdefs_fn);

    if ( subdom_tags.find( pd.omega1() ) == subdom_tags.end() )
    {
        std::cout << "Subdomain " << pd.omega1() << " not present in mesh" << std::endl;
        exit(-1);
    }

    if ( subdom_tags.find( pd.omega2() ) == subdom_tags.end() )
    {
        std::cout << "Subdomain " << pd.omega2() << " not present in mesh" << std::endl;
        exit(-1);
    }

    size_t internal_tag = pd.omega1();
    disk::hho_degree_info hdi(degree);

    dof_bases db = compute_dof_bases(msh, hdi, internal_tag);
    offsets_g2s g2s;
    compute_g2s(msh, internal_tag, g2s);

    auto cbs = disk::scalar_basis_size(hdi.cell_degree(), Mesh::dimension);
    auto fbs = disk::scalar_basis_size(hdi.face_degree(), Mesh::dimension-1);

    /* Global system data */
    std::cout << "Total system DoFs: " << db.total_dofs << std::endl;
    Eigen::SparseMatrix<T> LHS(db.total_dofs, db.total_dofs);
    Eigen::Matrix<T, Eigen::Dynamic, 1> RHS = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(db.total_dofs);
    using triplet = Eigen::Triplet<T>;
    std::vector< triplet > triplets;

    /* Lambda to assemble the laplacian contribution */
    auto asm_lapl = [&](const Mesh& msh, const typename Mesh::cell_type& cl,
        disk::dynamic_matrix<T>& lhs, disk::dynamic_vector<T>& rhs, size_t tag) {
        std::vector<size_t> l2g( lhs.rows() );
        size_t clofs = offset(msh, cl);
        size_t clbase;
        if (tag == internal_tag) 
            clbase = db.cells_int_base + cbs*g2s.cells_int_g2s.at(clofs);
        else
            clbase = db.cells_ext_base + cbs*g2s.cells_ext_g2s.at(clofs);

        for (size_t i = 0; i < cbs; i++)
            l2g[i] = clbase+i;

        size_t face_i = 0;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            size_t fcofs = offset(msh, fc);
            size_t fcbase;
            if (tag == internal_tag) 
                fcbase = db.faces_int_base + fbs*g2s.faces_int_g2s.at(fcofs);
            else
                fcbase = db.faces_ext_base + fbs*g2s.faces_ext_g2s.at(fcofs);

            for(size_t i = 0; i < fbs; i++)
                l2g[cbs + face_i*fbs + i] = fcbase+i;

            face_i++;
        }

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            for (size_t j = 0; j < lhs.cols(); j++)
            {
                assert(l2g[i] < LHS.rows());
                assert(l2g[j] < LHS.cols());
                triplets.push_back( triplet( l2g[i], l2g[j], lhs(i,j) ) );
            }
        }
                
        RHS.segment(clbase, cbs) = rhs;
    };

    /* Lambda to assemble integral global constraint */
    auto asm_constraint = [&](const Mesh& msh, const typename Mesh::cell_type& cl,
        disk::dynamic_vector<T>& average, size_t tag) {
        size_t clofs = offset(msh, cl);
        size_t clbase;
        if (tag == internal_tag)
            clbase = db.cells_int_base + cbs*g2s.cells_int_g2s.at(clofs);
        else
            clbase = db.cells_ext_base + cbs*g2s.cells_ext_g2s.at(clofs);

        for (size_t i = 0; i < cbs; i++)
        {
            assert( db.global_constrain_base < LHS.rows() and db.global_constrain_base < LHS.cols() );
            assert( clbase+i < LHS.rows() and clbase+i < LHS.cols() );
            triplets.push_back( triplet(db.global_constrain_base, clbase+i, average[i]) );
            triplets.push_back( triplet(clbase+i, db.global_constrain_base, average[i]) );
        }
    };

    /* Lambda to assemble lagrange multipliers */
    auto asm_lagrange = [&](const Mesh& msh, const typename Mesh::face_type& fc,
        disk::dynamic_matrix<T>& M, disk::dynamic_vector<T>& m, size_t tag) {
        size_t fcofs = offset(msh, fc);
        size_t lfcbase = db.multiplier_base + fbs*g2s.interface_g2s.at(fcofs);
        double sign = (tag == internal_tag) ? -1. : +1.;
        size_t fcbase;
        if (tag == internal_tag)
            fcbase = db.faces_int_base + fbs*g2s.faces_int_g2s.at(fcofs);
        else
            fcbase = db.faces_ext_base + fbs*g2s.faces_ext_g2s.at(fcofs);

        for (size_t i = 0; i < fbs; i++)
        {
            for (size_t j = 0; j < fbs; j++)
            {
                /* M is symmetric */
                assert(fcbase+i < LHS.rows() and fcbase+i < LHS.cols() );
                assert(lfcbase+j < LHS.rows() and lfcbase+j < LHS.cols() );
                triplets.push_back( triplet(fcbase+i, lfcbase+j, sign*M(i,j)) );
                triplets.push_back( triplet(lfcbase+i, fcbase+j, -sign*M(i,j)) );
            }
        }
        if (tag == internal_tag)
            RHS.segment(lfcbase, fbs) = m;
    };

    /* Lambda to assemble interface fluxes */
    auto asm_interface_flux = [&](const Mesh& msh, const typename Mesh::face_type& fc,
        disk::dynamic_vector<T>& g1, size_t tag) {
        size_t fcofs = offset(msh, fc);
        size_t fcbase;
        if (tag != internal_tag)
        {
            fcbase = db.faces_ext_base + fbs*g2s.faces_ext_g2s.at(fcofs);
            RHS.segment(fcbase, fbs) -= g1;
        }
    };

    /* Lambda to assemble boundary fluxes */
    auto asm_boundary_flux = [&](const Mesh& msh, const typename Mesh::face_type& fc,
        disk::dynamic_vector<T>& g2) {
        size_t fcofs = offset(msh, fc);
        size_t fcbase = db.faces_ext_base + fbs*g2s.faces_ext_g2s.at(fcofs);
        RHS.segment(fcbase, fbs) += g2;
    };

    /* Compute averages to check */
    double volume1 = 0.0;
    double volume2 = 0.0;
    double omega1_integral = 0.0;
    double omega2_integral = 0.0;
    for (auto& cl : msh) {
        auto qps = integrate(msh, cl, 10);

        auto di = msh.domain_info(cl);
        if ( di.tag() == internal_tag ) {
            for (auto& qp : qps) {
                volume1 += qp.weight();
                omega1_integral += qp.weight() * pd.u1(msh, qp.point());
            }
        }
        else {
            for (auto& qp : qps) {
                volume2 += qp.weight();
                omega2_integral += qp.weight() * pd.u2(msh, qp.point());
            }
        }
    }

    std::cout << "mean of u₁ on Ω₁: " << std::setprecision(15);
    std::string comp1 = pd.compensate_average1() ? " (compensation ON)" : "";
    std::cout << omega1_integral/volume1 << comp1 << std::endl;
    std::cout << "mean of u₂ on Ω₂: " << std::setprecision(15);
    std::string comp2 = pd.compensate_average2() ? " (compensation ON)" : "";
    std::cout << omega2_integral/volume2 << comp2 << std::defaultfloat << std::endl;

    if ( pd.compensate_average1() )
        pd.avgcomp1( omega1_integral/volume1 );

    if ( pd.compensate_average2() )
        pd.avgcomp2( omega2_integral/volume2 );

    double omega1_integral_aftercomp = 0.0;
    double omega2_integral_aftercomp = 0.0;

    for (auto& cl : msh) {
        auto qps = integrate(msh, cl, 10);

        auto di = msh.domain_info(cl);
        if ( di.tag() == internal_tag ) {
            for (auto& qp : qps) {
                omega1_integral_aftercomp += qp.weight() * pd.u1(msh, qp.point());
            }
        }
        else {
            for (auto& qp : qps) {
                omega2_integral_aftercomp += qp.weight() * pd.u2(msh, qp.point());
            }
        }
    }

    if ( pd.compensate_average1() ) {
        std::cout << "mean of u₁ on Ω₁: " << std::setprecision(15);
        std::cout << omega1_integral_aftercomp/volume1 << " (after compensation)";
        std::cout << std::endl;
    }

    if ( pd.compensate_average2() ) {
        std::cout << "mean of u₂ on Ω₂: " << std::setprecision(15);
        std::cout << omega2_integral_aftercomp/volume2 << " (after compensation)";
        std::cout << std::endl;
    }

    std::cout << "Assembling " << msh.cells_size() << " elements" << std::endl;

    timecounter tc;
    tc.tic();

    /* Assembly loop */
    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);

        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

        /* Standard HHO Laplacian + stabilization */
        auto cell_rhs_fun = [&](const typename Mesh::point_type& pt) -> auto {
            if ( di.tag() == internal_tag )
                return pd.f1(msh, pt);
            else
                return pd.f2(msh, pt);
        };
        disk::dynamic_matrix<T> lhs = gr.second + stab;
        disk::dynamic_vector<T> rhs = make_rhs(msh, cl, cb, cell_rhs_fun, 1);
        asm_lapl(msh, cl, lhs, rhs, di.tag());

        /* Pure Neumann global constraint */
        disk::dynamic_vector<T> average = disk::dynamic_vector<T>::Zero(cb.size());
        auto qps = integrate(msh, cl, hdi.cell_degree());
        for (auto& qp : qps) {
            auto phi = cb.eval_functions(qp.point());
            average += qp.weight() * phi;
        }
        asm_constraint(msh, cl, average, di.tag());

        /* Do the faces */
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            auto bi = msh.boundary_info(fc);

            if (bi.is_internal())
            {
                auto lm_fun = [&](const typename Mesh::point_type& pt) -> auto {
                    return pd.g(msh, pt);
                };

                auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                disk::dynamic_matrix<T> M = make_mass_matrix(msh, fc, fb);
                disk::dynamic_vector<T> m = make_rhs(msh, fc, fb, lm_fun, 1);
                asm_lagrange(msh, fc, M, m, di.tag());

                auto g1_fun = [&](const typename Mesh::point_type& pt) -> auto {
                    return pd.g1(msh, cl, fc, pt);
                };
                disk::dynamic_vector<T> g1 = make_rhs(msh, fc, fb, g1_fun, 1);
                asm_interface_flux(msh, fc, g1, di.tag());
            }
            else if (bi.is_boundary())
            {
                auto fb = make_scalar_monomial_basis(msh, fc, hdi.face_degree());

                auto g2_fun = [&](const typename Mesh::point_type& pt) -> auto {
                    return pd.g2(msh, cl, fc, pt);
                };
                disk::dynamic_vector<T> g2 = make_rhs(msh, fc, fb, g2_fun, 1);
                asm_boundary_flux(msh, fc, g2);
            }
        }
    }

    /* Initialize global matrix */
    LHS.setFromTriplets(triplets.begin(), triplets.end());
    std::cout << "Triplets: " << triplets.size() << std::endl;
    triplets.clear();

    std::cout << "Assembly time: " << tc.toc() << " seconds" << std::endl;

    std::cout << LHS.nonZeros() << " " << LHS.rows() << " " << LHS.cols() << std::endl;

    /* Run solver */
    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(LHS.rows());

    std::string solver_name = pd.solver();

    tc.tic();

    if (solver_name == "pardiso")
    {
        std::cout << "Running pardiso" << std::endl;
        pardiso_params<T> pparams;
        pparams.report_factorization_Mflops = true;
        pparams.out_of_core = PARDISO_OUT_OF_CORE_IF_NEEDED;

        bool success = mkl_pardiso(pparams, LHS, RHS, sol);
        if (!success)
        {
            std::cout << "Pardiso failed" << std::endl;
            return;
        }
    }
    else
    {
        std::cout << "Eigen SparseLU" << std::endl;
        SparseLU<SparseMatrix<T>, COLAMDOrdering<int> >   solver;
        solver.analyzePattern(LHS);
        /*
        if ( solver.info() != Eigen::Success )
        {
            std::cout << "SparseLU failed." << std::endl;
            return; 
        }
        */
        solver.factorize(LHS);
        sol = solver.solve(RHS); 
    }

    std::cout << "Solver time: " << tc.toc() << " seconds" << std::endl;

    tc.tic();

    /* Postprocess */
    T l2_error_sq = 0.0;
    T l2_error_mm_sq = 0.0;
    T l2_gradient_error_sq = 0.0;
    T l2_reconstruction_error_sq = 0.0;
    T l2_reconstruction_error_mm_sq = 0.0;
    T l2_trace_error_mm_sq = 0.0;

    std::vector<T> vec_u, vec_u_ref;
    vec_u.reserve( msh.cells_size() );
    vec_u_ref.reserve( msh.cells_size() );
    for (auto& cl : msh)
    {
        auto di = msh.domain_info(cl);
        size_t cell_ofs;
        if (di.tag() == internal_tag)
            cell_ofs = db.cells_int_base + cbs*g2s.cells_int_g2s.at( offset(msh, cl) );
        else
            cell_ofs = db.cells_ext_base + cbs*g2s.cells_ext_g2s.at( offset(msh, cl) );

        auto cb     = make_scalar_monomial_basis(msh, cl, hdi.cell_degree());
        auto gr     = make_scalar_hho_laplacian(msh, cl, hdi);
        auto stab   = make_scalar_hho_stabilization(msh, cl, gr.first, hdi);

        disk::dynamic_matrix<T> lhs = gr.second + stab;

        disk::dynamic_vector<T> solH = disk::dynamic_vector<T>::Zero(lhs.rows());
        solH.segment(0, cbs) = sol.segment(cell_ofs, cbs);
        vec_u.push_back( solH(0) );

        auto bar = barycenter(msh, cl);
        if (di.tag() == internal_tag)
            vec_u_ref.push_back( pd.u1(msh, bar) );
        else
            vec_u_ref.push_back( pd.u2(msh, bar) );

        size_t face_i = 0;
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            size_t fcofs = offset(msh, fc);
            size_t fcbase;
            if (di.tag() == internal_tag)
            {
                fcbase = db.faces_int_base + fbs*g2s.faces_int_g2s.at(fcofs);
             
                auto xi_fun = [&](const typename Mesh::point_type& pt) -> auto {
                    return pd.xi(msh, cl, fc, pt);
                };
             
                auto bi = msh.boundary_info(fc);
                if ( bi.is_internal() )
                {
                    auto multbase = db.multiplier_base + fbs*g2s.interface_g2s.at(fcofs);
                    disk::dynamic_vector<T> solM = sol.segment(multbase, fbs);
                    auto fb = disk::make_scalar_monomial_basis(msh, fc, hdi.face_degree());
                    disk::dynamic_matrix<T> massM = disk::make_mass_matrix(msh, fc, fb);
                    disk::dynamic_vector<T> rhsM = disk::make_rhs(msh, fc, fb, xi_fun, 1);
                    disk::dynamic_vector<T> projM = massM.llt().solve(rhsM);
                    disk::dynamic_vector<T> diffM = projM - solM;
                    l2_trace_error_mm_sq += diffM.dot(massM*diffM)*diameter(msh, fc);
                }
            }
            else
                fcbase = db.faces_ext_base + fbs*g2s.faces_ext_g2s.at(fcofs);

            solH.segment(cbs + face_i*fbs, fbs) = sol.segment(fcbase, fbs);

            face_i++;
        }
        
        auto rb = make_scalar_monomial_basis(msh, cl, hdi.cell_degree()+1);
        disk::dynamic_vector<T> solR = disk::dynamic_vector<T>::Zero(rb.size());
        disk::dynamic_vector<T> avefix = disk::dynamic_vector<T>::Zero(rb.size());
        solR.segment(1, rb.size()-1) = gr.first*solH;
        avefix.head(cb.size()) = solH.head(cb.size());
        avefix.tail(rb.size()-1) -= solR.tail(rb.size()-1);

        auto qps = integrate(msh, cl, 2*hdi.cell_degree()+3);

        T avg = 0.0;
        for (auto& qp : qps)
        {
            auto rphi = rb.eval_functions(qp.point());
            for (size_t i = 0; i < rb.size(); i++)
                avg += qp.weight()*avefix(i)*rphi(i);
        }
        avg /= measure(msh,cl);


        solR(0) = avg;

        auto sol_fun = [&](const typename Mesh::point_type& pt) -> auto {
            if ( di.tag() == internal_tag )
                return pd.u1(msh, pt);
            else
                return pd.u2(msh, pt);
        };

        for (auto& qp : qps)
        {
            auto phi = cb.eval_functions(qp.point());
            auto uh_val = solH.segment(0, cbs).dot(phi);

            auto rphi = rb.eval_functions(qp.point());
            auto rh_val = solR.dot(rphi);
            
            auto u_val = sol_fun(qp.point());
            l2_error_sq += qp.weight() * std::pow(uh_val - u_val, 2.);
            l2_reconstruction_error_sq += qp.weight() * std::pow(rh_val - u_val, 2.);
        }



        disk::dynamic_matrix<T> massR = disk::make_mass_matrix(msh, cl, rb);
        disk::dynamic_vector<T> projH = disk::project_function(msh, cl, hdi, sol_fun, 1);
        disk::dynamic_vector<T> projR = disk::project_function(msh, cl, hdi.cell_degree()+1, sol_fun, 1);
        
        disk::dynamic_matrix<T> massC = massR.block(0,0,cbs,cbs);

        disk::dynamic_vector<T> diffH = projH - solH;
        disk::dynamic_vector<T> diffC = diffH.segment(0, cbs);
        disk::dynamic_vector<T> diffR = projR - solR;

        l2_gradient_error_sq += diffH.dot(lhs*diffH);
        l2_error_mm_sq += diffC.dot(massC*diffC);
        l2_reconstruction_error_mm_sq += diffR.dot(massR*diffR);
    }

    std::cout << "Postpro time: " << tc.toc() << " seconds" << std::endl;

    std::cout << "h = " << disk::average_diameter(msh) << std::endl;
    std::cout << "L2 error: " << std::sqrt(l2_error_sq) << std::endl;
    std::cout << "L2 error MM: " << std::sqrt(l2_error_mm_sq) << std::endl;
    std::cout << "Reconstruction L2 error: " << std::sqrt(l2_reconstruction_error_sq) << std::endl;
    std::cout << "Reconstruction L2 error MM: " << std::sqrt(l2_reconstruction_error_mm_sq) << std::endl;
    std::cout << "Gradient error: " << std::sqrt(l2_gradient_error_sq) << std::endl;
    std::cout << "Trace error: " << std::setprecision(15) << std::sqrt(l2_trace_error_mm_sq) << std::endl;


    std::cout << "h ; L_2 u ; L_2 R(u) ; Grad ; trace" << std::endl;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << disk::average_diameter(msh) << " ; " << std::sqrt(l2_error_mm_sq);
    std::cout << " ; " << std::sqrt(l2_reconstruction_error_mm_sq);
    std::cout << " ; " << std::sqrt(l2_gradient_error_sq) << " ; ";
    std::cout << std::sqrt(l2_trace_error_mm_sq) << std::endl;

    disk::silo_database silo_db;
    silo_db.create("transmission.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<T> u("u", vec_u);
    silo_db.add_variable("mesh", u);

    disk::silo_zonal_variable<T> u_ref("u_ref", vec_u_ref);
    silo_db.add_variable("mesh", u_ref);
}

static void
usage(const char *progname)
{
    std::cout << "Usage: " << progname << " <options>" << std::endl;
    std::cout << " -k <integer>     Polynomial degree" << std::endl;
    std::cout << " -n               Crash if NaNs are generated" << std::endl;
    std::cout << " -m <string>      Mesh file name" << std::endl;
    std::cout << " -p <string>      Problem definition file name" << std::endl;
    std::cout << " -s <float>       Scale factor (Only FVCA5)" << std::endl;
    std::cout << " -x <float>       X shift (Only FVCA5)" << std::endl;
    std::cout << " -y <float>       Y shift (Only FVCA5)" << std::endl;
}

int main(int argc, char * const *argv)
{
    using T = double;

    int degree = 0;
    const char *mesh_filename = nullptr;
    const char *pbdefs_filename = nullptr;
    T scalefactor = 1.0;
    T displacement_x = 0.0;
    T displacement_y = 0.0;

    int ch;
    while ( (ch = getopt(argc, argv, "k:nm:p:s:x:y:")) != -1 )
    {
        switch(ch)
        {
            case 'k': /* Polynomial degree */
                degree = std::stoi(optarg);
                break;

            case 'n': /* Make program crash if NaNs are generated */
                _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
                break;

            case 'm': /* Mesh filename */
                mesh_filename = optarg;
                break;

            case 'p': /* Problem definitions filename */
                pbdefs_filename = optarg;
                break;

            case 's': /* Scale factor (only for FVCA5) */
                scalefactor = std::stod(optarg);
                break;

            case 'x': /* Displacement in X (only for FVCA5) */
                displacement_x = std::stod(optarg);
                break;

            case 'y': /* Displacement in X (only for FVCA5) */
                displacement_y = std::stod(optarg);
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                usage(argv[0]);
                return -1;
        }
    }

    if (!mesh_filename)
    {
        std::cout << "Please specify mesh filename (-m)" << std::endl;
        usage(argv[0]);
        return -1;
    }

    if (!pbdefs_filename)
    {
        std::cout << "Please specify problem definitions (-p)" << std::endl;
        usage(argv[0]);
        return -1;
    }

#ifdef HAVE_GMSH
    /* GMSH 2D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo2s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 2D simplicials" << std::endl;
        disk::simplicial_mesh<T,2> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,2> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        lt_solver(msh, degree, pbdefs_filename);
    }

    /* GMSH 3D simplicials */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3s$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D simplicials" << std::endl;
        disk::simplicial_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::simplicial_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        lt_solver(msh, degree, pbdefs_filename);
    }

    /* GMSH 3D generic */
    if (std::regex_match(mesh_filename, std::regex(".*\\.geo3g$") ))
    {
        std::cout << "Guessed mesh format: GMSH 3D generic" << std::endl;
        disk::generic_mesh<T,3> msh;
        disk::gmsh_geometry_loader< disk::generic_mesh<T,3> > loader;
        
        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        lt_solver(msh, degree, pbdefs_filename);
    }
#endif /* HAVE_GMSH */

    /* FVCA5 2D */
    if (std::regex_match(mesh_filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;
        disk::generic_mesh<T,2> msh;
        disk::fvca5_mesh_loader<T,2> loader;

        loader.read_mesh(mesh_filename);
        loader.populate_mesh(msh);

        using point_type = typename disk::generic_mesh<T,2>::point_type;
        auto tr = [&](const point_type& pt) -> auto {
            auto npx = pt.x()*scalefactor + displacement_x;
            auto npy = pt.y()*scalefactor + displacement_y;
            return point_type(npx, npy);
        };

        msh.transform(tr);

        lt_solver(msh, degree, pbdefs_filename);
        return 0;
    }

    return 0;
}
