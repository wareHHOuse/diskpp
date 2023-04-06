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

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/quadratures/quadratures.hpp"

using namespace disk;
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

template<typename Mesh>
int test_basis_functions(const Mesh& msh)
{
    size_t degree = 8;
    using mesh_type = Mesh;
    using T = typename mesh_type::coordinate_type;
    using point_type = typename mesh_type::point_type;

    auto cl = *msh.cells_begin();

    auto f = [&](const point_type& pt) {
        return std::sin(pt.x() + pt.y());
    };

    auto grad_f = [&](const point_type& pt) {
        Eigen::Matrix<T, 2, 1> ret;
        ret(0) = std::cos(pt.x() + pt.y());
        ret(1) = std::cos(pt.x() + pt.y());
        return ret;
    };

    auto basis = scaled_monomial_basis(msh, cl, degree);

    L2_projector π(msh, cl, basis);
    auto dofs_f = π(f);

    /* Error on the projected function */
    auto errf = [&](const point_type& pt) {
        T s = f(pt) - eval(pt, dofs_f, basis);
        return s*s;
    };
    auto fun_error = integrate(msh, cl, 2*degree+1, errf); 

    /* Norm of the projected function */
    auto normf = [&](const point_type& pt) {
        T s = f(pt);
        return s*s;
    };
    auto fun_norm = integrate(msh, cl, 2*degree+1, normf);

    /* Error on the gradient of the projected function */
    auto grad_error = integrate(msh, cl, 2*degree+1, [&](const point_type& pt) {
        Eigen::Matrix<T,2,1> s = grad_f(pt) - eval(pt, dofs_f, grad(basis));
        return s.dot(s);
    });

    /* Norm of the gradient of the projected function */
    auto grad_norm = integrate(msh, cl, 2*degree+1, [&](const point_type& pt) {
        auto s = grad_f(pt);
        return s.dot(s);
    });

    auto fe = 100.0*std::sqrt(fun_error/fun_norm);
    auto ge = 100.0*std::sqrt(grad_error/grad_norm);

    bool fepass = false, gepass = false;
    if (fe < 0.012834) fepass = true;
    if (ge < 0.13469) gepass = true;

    auto passfail = [](bool pass) {
        return pass ? "[PASS]" : "[FAIL]";
    };

    std::cout << __FILE__ << std::endl;
    std::cout << "  Function relative error = " << fe << "% " << passfail(fepass) << std::endl;
    std::cout << "  Gradient relative error = " << ge << "% " << passfail(gepass) << std::endl;

    auto fcs = faces(msh, cl);
    for (auto& fc : fcs)
    {
        auto face_basis = scaled_monomial_basis(msh, fc, degree);
        L2_projector πF(msh, fc, face_basis);
        auto dofs_f_F = πF(f);
        T fun_error_f = 0.0;
        T fun_norm_f = 0.0;
        auto fqps = integrate(msh, fc, 2*degree+2);
        for (auto& qp : fqps) {
            T diff = f(qp.point()) - eval(qp.point(), dofs_f_F, face_basis);
            fun_error_f += qp.weight() * (diff*diff);
            fun_norm_f += qp.weight() * normf(qp.point());
        }
        std::cout << "    Face Relerr: " << 100*std::sqrt(fun_error_f/fun_norm_f) << std::endl;
    }


    return not (fepass and gepass);
}

int main(void)
{
    using T = double;

    simplicial_mesh<T,2> msh_simp;
    make_single_element_mesh(msh_simp, {0,0}, {1,0}, {0,1});
    test_basis_functions(msh_simp);

    cartesian_mesh<T,2> msh_cart;
    make_single_element_mesh(msh_cart, {0.5, 0.5}, 1.0, M_PI);
    test_basis_functions(msh_cart);

    return 0;
}
