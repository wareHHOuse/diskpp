/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

/* Test for 1D raw monomial bases: interpolate a sine, measure the error
 * on the interpolation and its gradient. */

#include <iostream>
#include <fstream>

#include <silo.h>

#include "diskpp/bases/bases_raw.hpp"
//#include "diskpp/quadratures/quad_raw_gauss.hpp"
#include "diskpp/quadratures/quad_raw_dunavant.hpp"

using namespace disk;
using namespace disk::quadrature;
using namespace disk::basis;

template<typename T>
bool add_mesh(DBfile *db, const std::vector<point<T,2>>& pts, const std::string& name)
{
    std::cout << "Point mesh: " << pts.size() << std::endl;
    std::vector<double> x_coords, y_coords;
    x_coords.reserve(pts.size());
    y_coords.reserve(pts.size());

    for (auto& pt : pts) 
    {
        x_coords.push_back(pt.x());
        y_coords.push_back(pt.y());
    }

    double *coords[] = { x_coords.data(), y_coords.data() };

    int ndims = 2;
    int nels = pts.size();

    auto err = DBPutPointmesh(db, name.c_str(), ndims, coords, nels, DB_DOUBLE, nullptr);
    return err == 0;
}

template<typename T>
bool add_mesh(DBfile *db, const std::vector<point<T,3>>& pts, const std::string& name)
{
    std::cout << "Point mesh: " << pts.size() << std::endl;
    std::vector<double> x_coords, y_coords, z_coords;
    x_coords.reserve(pts.size());
    y_coords.reserve(pts.size());
    z_coords.reserve(pts.size());

    for (auto& pt : pts) 
    {
        x_coords.push_back(pt.x());
        y_coords.push_back(pt.y());
        z_coords.push_back(pt.z());
    }

    double *coords[] = { x_coords.data(), y_coords.data(), z_coords.data() };

    int ndims = 3;
    int nels = pts.size();

    auto err = DBPutPointmesh(db, name.c_str(), ndims, coords, nels, DB_DOUBLE, nullptr);
    return err == 0;
}

bool add_variable(DBfile *db, const std::string& mesh_name, const std::string& var_name,
    const std::vector<double>& var)
{
    int nels = var.size();

    auto err = DBPutPointvar1(db, var_name.c_str(), mesh_name.c_str(),
        var.data(), nels, DB_DOUBLE, nullptr);

    return err == 0;
}

template<typename T, size_t DIM>
DBfile * open_silo(const std::vector<point<T,DIM>>& pts)
{
    std::string output_filename = "monomial_bases_2D_i2.silo";
    if (DIM == 3)
        output_filename = "monomial_bases_2D_i3.silo";
    
    DBfile *db = DBCreate(output_filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
    if (!db) {
        std::cout << "Can't open SILO DB" << std::endl;
        return nullptr;
    }

    add_mesh(db, pts, "mesh");

    return db;
}

template<typename T, size_t DIM>
struct test_data;

template<typename T>
struct test_data<T,2>
{
    auto operator()(const point<T,2>& pt) const {
        return std::sin(pt.x() + pt.y());
    }

    auto grad(const point<T,2>& pt) const {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = std::cos(pt.x()+pt.y());
        ret(1) = std::cos(pt.x()+pt.y());
        return ret;
    }
};

template<typename T>
struct test_data<T,3>
{
    auto operator()(const point<T,3>& pt) const {
        return std::sin(pt.x() + pt.y() + pt.z());
    }

    auto grad(const point<T,3>& pt) const {
        Eigen::Matrix<T,3,1> ret;
        ret(0) = std::cos(pt.x()+pt.y()+pt.z());
        ret(1) = std::cos(pt.x()+pt.y()+pt.z());
        ret(2) = std::cos(pt.x()+pt.y()+pt.z());
        return ret;
    };
};

template<typename T, size_t DIM>
int test_basis_functions(const point<T,DIM>& a, const point<T,DIM>& b, const point<T,DIM>& c)
{
    size_t degree = 8;

    using point_type = point<T,DIM>;

    auto bar = (a+b+c)/3.0;

    scaled_monomial_basis<T, T, DIM, 2> basis(b.to_vector(), c.to_vector(), bar, degree);

    test_data<T,DIM> td;
    auto fun = [&](const point_type& pt) {
        return td(pt);
    };

    auto grad = [&](const point_type& pt) {
        return td.grad(pt);
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mass =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(basis.size(), basis.size());

    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(basis.size());

    auto qps = dunavant(2*degree+2, a, b, c);
    for (auto& qp : qps) {
        auto phi = basis(qp.point());
        mass += qp.weight() * phi * phi.transpose();
        rhs += qp.weight() * fun(qp.point()) * phi;
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> dofs = mass.ldlt().solve(rhs);





    std::vector<point_type>     points;
    std::vector<std::vector<T>> funs, grads_x, grads_y;
    std::vector<T> pf, pgx, pgy;
    funs.resize( basis.size() );
    grads_x.resize( basis.size() );
    grads_y.resize( basis.size() );

    T fun_error = 0.0, fun_norm = 0.0;
    T grad_error = 0.0, grad_norm = 0.0;
    for (auto& qp : qps) {
        points.push_back(qp.point());

        auto phi = basis(qp.point());
        auto fnum = dofs.dot(phi);
        auto fval = fun(qp.point());
        auto diff_fun = fval - fnum;
        pf.push_back(fnum);
        fun_error += qp.weight() * diff_fun * diff_fun;
        fun_norm += qp.weight() * fval * fval;

        auto dphi = basis.grad(qp.point());
        auto gnum = dphi.transpose()*dofs;
        auto gval = grad(qp.point());
        auto diff_grad = gval - gnum;
        std::cout << "******" << std::endl;
        std::cout << gnum.transpose() << std::endl;
        std::cout << gval.transpose() << std::endl;
        std::cout << diff_grad.transpose() << std::endl;
        pgx.push_back(gnum(0));
        pgy.push_back(gnum(1));
        grad_error += qp.weight() * diff_grad.dot(diff_grad);
        grad_norm += qp.weight() * gval.dot(gval);

        for (size_t i = 0; i < basis.size(); i++) {
            funs[i].push_back( phi(i) );
            grads_x[i].push_back( dphi(i,0) );
            grads_y[i].push_back( dphi(i,1) );
        }
    }

    DBfile *db = open_silo(points);
    if (db) {
        add_variable(db, "mesh", "pf", pf);
        add_variable(db, "mesh", "pgx", pgx);
        add_variable(db, "mesh", "pgy", pgy);
        for (size_t i = 0; i < funs.size(); i++) {
            add_variable(db, "mesh", "fun_" + std::to_string(i), funs[i]);

            std::string gradx_name = "grad_x_" + std::to_string(i);
            add_variable(db, "mesh", gradx_name.c_str(), grads_x[i]);
            
            std::string grady_name = "grad_y_" + std::to_string(i);
            add_variable(db, "mesh", grady_name.c_str(), grads_y[i]);
            
            std::string expr_name = "grad_" + std::to_string(i);
            std::string expr_definition = "{" + gradx_name + ", " + grady_name + "}";
            const char *name[] = { expr_name.c_str() };
            const char *def[] = { expr_definition.c_str() };
            auto expr_type = DB_VARTYPE_VECTOR;
            std::string defname = "def_" + std::to_string(i);
            DBPutDefvars(db, defname.c_str(), 1, name, &expr_type, def, NULL);
        }

        DBClose(db);
    }


    auto fe = 100.0*std::sqrt(fun_error/fun_norm);
    auto ge = 100.0*std::sqrt(grad_error/grad_norm);

    const double DBL_FE_MAX_ERR = 2.72e-14;
    const double DBL_GE_MAX_ERR = 1.56e-06;

    if (true or fe > DBL_FE_MAX_ERR or ge > DBL_GE_MAX_ERR) {
        std::cout << __FILE__ << ": Test FAILED: fe = " << fe << "%, ge = " << ge << "%" << std::endl;
        return 1; 
    }

    return 0;
}

int main(void)
{
    using T = double;

    //point<T,2> a2{0, 0};
    //point<T,2> b2{M_PI, 0};
    //point<T,2> c2{0, M_PI};

    point<T,2> a2{0, 0};
    point<T,2> b2{M_PI, M_PI};
    point<T,2> c2{-M_PI, M_PI};

    test_basis_functions(a2, b2, c2);

    //point<T,3> a3{0, 0, 0};
    //point<T,3> b3{M_PI, 0, 1};
    //point<T,3> c3{0, M_PI, 1};

    point<T,3> a3{0, 0, 0};
    point<T,3> b3{M_PI, M_PI, 1};
    point<T,3> c3{-M_PI, M_PI, 1};

    test_basis_functions(a3, b3, c3);
}