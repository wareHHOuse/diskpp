/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#include <iostream>
#include <regex>

#include "loaders/loader.hpp"
#include "hho/hho.hpp"
#include "solvers/solver.hpp"

using namespace Eigen;


size_t basis_size(size_t k, size_t d)
{
    switch (d)
    {
        case 1:
            return k+1;

        case 2:
            return (k+2)*(k+1)/2;

        case 3:
            return (k+3)*(k+2)*(k+1)/6;

        default:
            throw std::invalid_argument("degree must be 1, 2 or 3");
    }
}

/*
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;
*/


template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 0)
{
    typedef Mesh                                                        mesh_type;
    typedef typename mesh_type::cell                                    cell_type;
    typedef disk::quadrature<mesh_type, cell_type>                      cell_quadrature_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;

    cell_basis_type cb(degree);
    auto cbs = cb.size();

    cell_quadrature_type cq(2*(degree+di));

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);

    auto qps = cq.integrate(msh, cl);

    for (auto& qp : qps)
    {
        auto phi = cb.eval_functions(msh, cl, qp.point());
        ret += qp.weight() * phi * phi.transpose();
    }

    return ret;
}

template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree, size_t di = 0)
{
    typedef Mesh                                                        mesh_type;
    typedef typename mesh_type::face                                    face_type;
    typedef disk::quadrature<mesh_type, face_type>                      face_quadrature_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    face_basis_type fb(degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);

    face_quadrature_type fq(2*(degree+di));

    auto qps = fq.integrate(msh, fc);

    for (auto& qp : qps)
    {
        auto phi = fb.eval_functions(msh, fc, qp.point());
        ret += qp.weight() * phi * phi.transpose();
    }

    return ret;
}

template<typename Mesh, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
         size_t degree, const Function& f, size_t di = 0)
{
    typedef Mesh                                                        mesh_type;
    typedef typename mesh_type::cell                                    cell_type;
    typedef disk::quadrature<mesh_type, cell_type>                      cell_quadrature_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    using T = typename Mesh::coordinate_type;

    cell_basis_type cb(degree);
    auto cbs = cb.size();

    cell_quadrature_type cq(2*(degree+di));

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    auto qps = cq.integrate(msh, cl);

    for (auto& qp : qps)
    {
        auto phi = cb.eval_functions(msh, cl, qp.point());
        ret += qp.weight() * phi * f(qp.point());
    }

    return ret;
}

template<typename Mesh, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const typename Mesh::face_type& fc,
         size_t degree, const Function& f, size_t di = 0)
{
    typedef Mesh                                                        mesh_type;
    typedef typename mesh_type::face                                    face_type;
    typedef disk::quadrature<mesh_type, face_type>                      face_quadrature_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;
    using T = typename Mesh::coordinate_type;

    face_basis_type fb(degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(fbs);

    face_quadrature_type fq(2*(degree+di));
    auto qps = fq.integrate(msh, fc);

    for (auto& qp : qps)
    {
        auto phi = fb.eval_functions(msh, fc, qp.point());
        ret += qp.weight() * phi * f(qp.point());
    }

    return ret;
}

template<typename Mesh, typename Function, typename BQData>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const typename Mesh::cell_type& cl,
                 const BQData& bqd, const Function& f, size_t di = 0)
{
    using T = typename Mesh::coordinate_type;

    auto celdeg = bqd.cell_degree();
    auto facdeg = bqd.face_degree();

    auto cbs = basis_size(celdeg, Mesh::dimension);
    auto fbs = basis_size(facdeg, Mesh::dimension-1);

    auto fcs = faces(msh, cl);
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs+fcs.size()*fbs);

    Matrix<T, Dynamic, Dynamic> cell_mm = make_mass_matrix(msh, cl, celdeg, di);
    Matrix<T, Dynamic, 1> cell_rhs = make_rhs(msh, cl, celdeg, f, di);
    ret.block(0, 0, cbs, 1) = cell_mm.llt().solve(cell_rhs);

    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        Matrix<T, Dynamic, Dynamic> face_mm = make_mass_matrix(msh, fc, facdeg, di);
        Matrix<T, Dynamic, 1> face_rhs = make_rhs(msh, fc, facdeg, f, di);
        ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.llt().solve(face_rhs);
    }

    return ret;
}










template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto itor = std::lower_bound(msh.cells_begin(), msh.cells_end(), cl);
    if ( itor == msh.cells_end() )
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells_begin(), itor);
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto itor = std::lower_bound(msh.faces_begin(), msh.faces_end(), fc);
    if ( itor == msh.faces_end() )
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces_begin(), itor);
}



/* Assembler for the obstacle problem (see "Bubbles enriched quadratic finite
 * element method for the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi",
 * eqn. 5.1 onwards) */


template<typename Mesh, typename BQData>
class obstacle_assembler
{
    using T = typename Mesh::coordinate_type;
    std::vector<size_t>                 face_ct, A_ct, B_ct; //compress tables
    std::vector<size_t>                 face_et, A_et, B_et; //expand tables
    std::vector<bool>                   is_in_set_A;

    std::vector< Triplet<T> >           triplets;

    size_t      num_all_faces, num_dirichlet_faces, num_other_faces;
    size_t      num_all_cells, num_A_cells, num_I_cells;

    const BQData& bqd;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    obstacle_assembler(const Mesh& msh, const std::vector<bool>& in_A, const BQData& p_bqd)
        : bqd(p_bqd), is_in_set_A(in_A)
    {
        auto is_dirichlet = [&](const typename Mesh::face_type& fc) -> bool {
            return msh.is_boundary(fc);
        };

        num_all_faces = msh.faces_size();
        num_dirichlet_faces = std::count_if(msh.faces_begin(), msh.faces_end(), is_dirichlet);
        num_other_faces = num_all_faces - num_dirichlet_faces;

        assert( is_in_set_A.size() == msh.cells_size() );

        num_all_cells = msh.cells_size();
        num_A_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), true);
        num_I_cells = std::count(is_in_set_A.begin(), is_in_set_A.end(), false);

        /* Make A tables: keep the unknowns of cells in set I */
        A_ct.resize( num_all_cells );
        A_et.resize( num_I_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( !is_in_set_A.at(i) )
            {
                A_ct.at(i) = co;
                A_et.at(co) = i;
                co++;
            }
        }

        /* Make face tables */
        face_ct.resize( num_all_faces );
        face_et.resize( num_other_faces );
        for (size_t i = 0, co = 0; i < num_all_faces; i++)
        {
            auto fc = *std::next(msh.faces_begin(), i);
            if ( !is_dirichlet(fc) )
            {
                face_ct.at(i) = co;
                face_et.at(co) = i;
                co++;
            }
        }

        /* Make B tables: keep the unknowns of cells in set A */
        B_ct.resize( num_all_cells );
        B_et.resize( num_A_cells );
        for (size_t i = 0, co = 0; i < num_all_cells; i++)
        {
            auto cl = *std::next(msh.cells_begin(), i);
            if ( is_in_set_A.at(i) )
            {
                B_ct.at(i) = co;
                B_et.at(co) = i;
                co++;
            }
        }

        auto celdeg = bqd.cell_degree();
        auto facdeg = bqd.face_degree();

        auto cbs = basis_size(celdeg, Mesh::dimension);
        auto fbs = basis_size(facdeg, Mesh::dimension-1);

        auto system_size = cbs * (num_I_cells + num_A_cells) + fbs * num_other_faces;

        assert( system_size == cbs * msh.cells_size() + fbs * num_other_faces );

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void dump_tables() const
    {
        std::cout << "Compress table for A: " << std::endl;
        for (size_t i = 0; i < A_ct.size(); i++)
            std::cout << i << " -> " << A_ct.at(i) << std::endl;

        std::cout << "Compress table for faces: " << std::endl;
        for (size_t i = 0; i < face_ct.size(); i++)
            std::cout << i << " -> " << face_ct.at(i) << std::endl;

        std::cout << "Compress table for B: " << std::endl;
        for (size_t i = 0; i < B_ct.size(); i++)
            std::cout << i << " -> " << B_ct.at(i) << std::endl;
    }

    template<typename Function>
    void
    assemble(const Mesh& msh, const typename Mesh::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs,
             const Matrix<T, Dynamic, 1>& gamma, const Function& dirichlet_bf)
    {
        auto celdeg = bqd.cell_degree();
        auto facdeg = bqd.face_degree();

        auto cbs = basis_size(celdeg, Mesh::dimension);
        auto fbs = basis_size(facdeg, Mesh::dimension-1);

        std::vector<assembly_index> asm_map_row, asm_map_col;
        //asm_map.reserve(cbs + 4*fbs);

        /* Cell dofs local to global */
        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = A_ct.at(cell_offset) * cbs;
        bool cell_needs_asm_A   = !is_in_set_A.at(cell_offset);
        bool cell_needs_asm_B   = is_in_set_A.at(cell_offset);
        
        for (size_t i = 0; i < cbs; i++)
        {
            asm_map_row.push_back( assembly_index(cell_offset+i, true) );
            asm_map_col.push_back( assembly_index(cell_LHS_offset+i, cell_needs_asm_A) );
        }

        /* Face dofs local to global */
        auto fcs = faces(msh, cl);
        auto numfcs = fcs.size();
        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(cbs + numfcs*fbs);
        for (size_t face_i = 0; face_i < numfcs; face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = offset(msh, fc);
            auto face_LHS_offset_rows = cbs * num_all_cells + face_ct.at(face_offset)*fbs;
            auto face_LHS_offset_cols = cbs * num_I_cells + face_ct.at(face_offset)*fbs;

            bool dirichlet = msh.is_boundary(fc);

            for (size_t i = 0; i < fbs; i++)
            {
                asm_map_row.push_back( assembly_index(face_LHS_offset_rows+i, !dirichlet) );
                asm_map_col.push_back( assembly_index(face_LHS_offset_cols+i, !dirichlet) );
            }

            if (dirichlet)
            {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, facdeg, dirichlet_bf);
                dirichlet_data.block(cbs+face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
        }

        //assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if ( !asm_map_row[i].assemble() )
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map_col[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map_row[i], asm_map_col[j], lhs(i,j)) );
                else
                {
                    if (j < cbs)
                        RHS(asm_map_row[i]) -= lhs(i,j)*gamma(cell_offset);
                    else
                        RHS(asm_map_row[i]) -= lhs(i,j)*dirichlet_data(j);
                }
            }
        }
        
        
        /* Needed both in case A and I */
        RHS.block(cell_offset, 0, cbs, 1) += rhs.block(0, 0, cbs, 1);

        if (cell_needs_asm_B)
        {
            auto offset_row = cell_offset * cbs;
            auto offset_col = num_I_cells * cbs + num_other_faces * fbs + B_ct.at(cell_offset);
            triplets.push_back( Triplet<T>(offset_row, offset_col, 1.0) );
        }

    } // assemble_A()

    
    template<typename Function>
    void
    expand_solution(const Mesh& msh,
    const Matrix<T, Dynamic, 1>& solution, const Function& dirichlet_bf,
    const Matrix<T, Dynamic, 1>& gamma,
    Matrix<T, Dynamic, 1>& alpha, Matrix<T, Dynamic, 1>& beta)
    {
        auto celdeg = bqd.cell_degree();
        auto facdeg = bqd.face_degree();

        auto cbs = basis_size(celdeg, Mesh::dimension);
        auto fbs = basis_size(facdeg, Mesh::dimension-1);

        alpha.block(0, 0, num_all_cells*cbs, 1) = gamma;
        for (size_t i = 0; i < num_I_cells; i++)
        {
            auto exp_ofs = A_et.at(i);
            alpha.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(i*cbs, 0, cbs, 1);
        }

        beta = Matrix<T, Dynamic, 1>::Zero(msh.cells_size());
        for (size_t i = 0; i < num_A_cells; i++)
        {
            auto exp_ofs = B_et.at(i);
            beta.block(exp_ofs*cbs, 0, cbs, 1) = solution.block(num_I_cells*cbs + num_other_faces*fbs + i*cbs, 0, cbs, 1);
        }

        for (auto itor = msh.faces_begin(); itor != msh.faces_end(); itor++)
        {
            auto fc = *itor;
            auto face_ofs = offset(msh, fc);

            bool dirichlet = msh.is_boundary(fc);

            if (dirichlet)
            {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, facdeg);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, facdeg, dirichlet_bf);
                alpha.block(num_all_cells*cbs + face_ofs*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = cbs * num_I_cells + face_ct.at(face_offset)*fbs;
                alpha.block(num_all_cells*cbs + face_ofs*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }
    }
    
    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }
};

template<typename T, typename Mesh, typename BQData>
Matrix<T, Dynamic, 1>
take_local_data(const Mesh& msh, const typename Mesh::cell_type& cl, const BQData& bqd,
                const Matrix<T, Dynamic, 1>& expanded_solution)
{
    auto celdeg = bqd.cell_degree();
    auto facdeg = bqd.face_degree();

    auto cbs = basis_size(celdeg, Mesh::dimension);
    auto fbs = basis_size(facdeg, Mesh::dimension-1);

    auto cell_offset        = offset(msh, cl);
    auto cell_SOL_offset    = cell_offset * cbs;

    auto fcs = faces(msh, cl);
    auto numfcs = fcs.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + numfcs*fbs);
    ret.block(0, 0, cbs, 1) = expanded_solution.block(cell_SOL_offset, 0, cbs, 1);

    for (size_t face_i = 0; face_i < numfcs; face_i++)
    {
        auto fc = fcs[face_i];
        auto face_offset = offset(msh, fc);
        auto face_SOL_offset = cbs * msh.cells_size() + face_offset*fbs;
        ret.block(cbs+face_i*fbs, 0, fbs, 1) = expanded_solution.block(face_SOL_offset, 0, fbs, 1);
    }

    return ret;
}


template<typename Mesh, typename BQData>
auto make_obstacle_assembler(const Mesh& msh, const std::vector<bool>& in_A, const BQData& bqd)
{
    return obstacle_assembler<Mesh, BQData>(msh, in_A, bqd);
}












template<typename Mesh>
bool
hho_solver(const Mesh& msh)
{
    typedef Mesh                                       mesh_type;
    typedef typename mesh_type::scalar_type            scalar_type;
    typedef typename mesh_type::cell                   cell_type;
    typedef typename mesh_type::face                   face_type;

    typedef disk::quadrature<mesh_type, cell_type>      cell_quadrature_type;
    typedef disk::quadrature<mesh_type, face_type>      face_quadrature_type;

    typedef disk::scaled_monomial_scalar_basis<mesh_type, cell_type>    cell_basis_type;
    typedef disk::scaled_monomial_scalar_basis<mesh_type, face_type>    face_basis_type;

    typedef
    disk::basis_quadrature_data<mesh_type,
                                disk::scaled_monomial_scalar_basis,
                                disk::quadrature> bqdata_type;

    typedef disk::gradient_reconstruction_bq<bqdata_type>               gradrec_type;
    typedef disk::diffusion_like_stabilization_bq<bqdata_type>          stab_type;
    typedef disk::diffusion_like_static_condensation_bq<bqdata_type>    statcond_type;

    size_t          degree      = 1;

    auto bqd = bqdata_type(0, degree);

    auto gradrec    = gradrec_type(bqd);
    auto stab       = stab_type(bqd);
    auto statcond   = statcond_type(bqd);


    auto r0 = 0.7;
/*
    auto rhs_fun = [=](const typename Mesh::point_type& pt) -> scalar_type {
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() );

        if ( r > r0 )
            return -16 *r*r + 8*r0*r0;
        else
            return -8.0*( r0*r0*(r0*r0 + 1) ) + 8*r0*r0*r*r;
    };

    auto sol_fun = [=](const typename Mesh::point_type& pt) -> scalar_type {
        auto r = sqrt(pt.x()*pt.x() + pt.y()*pt.y());
        auto s = r*r - r0*r0;
        auto t = std::max(s, 0.0);
        return t*t;
    };
    
    auto bcs_fun = [&](const typename Mesh::point_type& pt) -> scalar_type {
        return sol_fun(pt);
    };

    auto obstacle_fun = [&](const typename Mesh::point_type& pt) -> scalar_type {
        return 0.0;
    };
*/

    auto rhs_fun = [=](const typename Mesh::point_type& pt) -> scalar_type {
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z() );

        if ( r > r0 )
            return -4 * (2*r*r + 3*(r*r - r0*r0));
        else
            return -8.0*r0*r0*(1 - r*r - r0*r0);
    };

    auto sol_fun = [=](const typename Mesh::point_type& pt) -> scalar_type {
        auto r = sqrt( pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z() );
        auto s = r*r - r0*r0;
        auto t = std::max(s, 0.0);
        return t*t;
    };
    
    auto bcs_fun = [&](const typename Mesh::point_type& pt) -> scalar_type {
        return sol_fun(pt);
    };

    auto obstacle_fun = [&](const typename Mesh::point_type& pt) -> scalar_type {
        return 0.0;
    };



    auto num_cells = msh.cells_size();
    auto num_faces = msh.faces_size();
    auto fbs = basis_size(degree, Mesh::dimension-1);
    auto num_alpha_dofs = num_cells + fbs*num_faces;

    Matrix<scalar_type, Dynamic, 1> alpha = Matrix<scalar_type, Dynamic, 1>::Zero( num_alpha_dofs );
    Matrix<scalar_type, Dynamic, 1> beta  = Matrix<scalar_type, Dynamic, 1>::Ones( num_cells );
    Matrix<scalar_type, Dynamic, 1> gamma = Matrix<scalar_type, Dynamic, 1>::Zero( num_cells );
    scalar_type c = 1.0;

    size_t quadrature_degree_increase = 1;


    size_t i = 0;
    for (auto& cl : msh)
    {
        auto bar = barycenter(msh, cl);
        gamma(i++) = obstacle_fun(bar);
    }

    size_t iter = 0;
    bool has_to_iterate = true;
    while ( has_to_iterate && iter < 50 )
    {
        std::cout << "Iteration " << iter+1 << std::endl;


        /* Compute the beta quantity (see "Bubbles enriched quadratic finite element method for
        *  the 3D-elliptic obstacle problem - S. Gaddam, T. Gudi", eqn. 5.1 onwards) */
        Matrix<scalar_type, Dynamic, 1> diff = beta + c * ( alpha.head(num_cells) - gamma );
        std::vector<bool> in_A;
        in_A.resize(diff.size());
        Matrix<scalar_type, Dynamic, 1> active = Matrix<scalar_type, Dynamic, 1>::Zero(diff.size());

        for (size_t i = 0; i < diff.size(); i++)
        {
            in_A.at(i) = (diff(i) < 0);
            active(i) = (diff(i) < 0);
        }

        auto assembler  = make_obstacle_assembler(msh, in_A, bqd);

        for (auto& cl : msh)
        {
            gradrec.compute(msh, cl);
            stab.compute(msh, cl, gradrec.oper);
            Matrix<scalar_type, Dynamic, 1> f = make_rhs(msh, cl, 0, rhs_fun, quadrature_degree_increase);
            dynamic_matrix<scalar_type> lc = gradrec.data + stab.data;

            assembler.assemble(msh, cl, lc, f, gamma, bcs_fun);
        }

        std::cout << std::endl;


        assembler.finalize();

        size_t systsz = assembler.LHS.rows();
        size_t nnz = assembler.LHS.nonZeros();

        std::cout << "Starting linear solver..." << std::endl;
        std::cout << " * Solving for " << systsz << " unknowns." << std::endl;
        std::cout << " * Matrix fill: " << 100.0*double(nnz)/(systsz*systsz) << "%" << std::endl;

        dynamic_vector<scalar_type> sol = dynamic_vector<scalar_type>::Zero(systsz);



        disk::solvers::pardiso_params<scalar_type> pparams;

        mkl_pardiso(pparams, assembler.LHS, assembler.RHS, sol);


        Matrix<scalar_type, Dynamic, 1> alpha_prev = alpha;
        assembler.expand_solution(msh, sol, sol_fun, gamma, alpha, beta);


        if ( (alpha_prev - alpha).norm() < 1e-7 )
            break;

        iter++;

    }

    scalar_type error = 0.0;
    for (auto& cl : msh)
    {
        gradrec.compute(msh, cl);
        stab.compute(msh, cl, gradrec.oper);
        dynamic_matrix<scalar_type> lc = gradrec.data + stab.data;

        Matrix<scalar_type, Dynamic, 1> local = take_local_data(msh, cl, bqd, alpha);
        Matrix<scalar_type, Dynamic, 1> proj  = project_function(msh, cl, bqd, sol_fun, quadrature_degree_increase);
        
        Matrix<scalar_type, Dynamic, 1> diff = local - proj;

        error += diff.dot(lc*diff);
    }
    
    std::cout << "Error: " << std::sqrt(error) << std::endl;
    
    return true;
}









template<typename MeshType>
void
process_mesh(const MeshType& msh)
{
    std::cout << "Mesh average h: " << mesh_h(msh) << std::endl;
    hho_solver(msh);
}


int main(int argc, char **argv)
{
    using RealType = double;

    char    *filename       = nullptr;
    int ch;


    if (argc == 1)
    {
        std::cout << "Please specify a 2D or 3D mesh" << std::endl;

        return 0;
    }

    filename = argv[1];
/*
    if (std::regex_match(filename, std::regex(".*\\.typ1$") ))
    {
        std::cout << "Guessed mesh format: FVCA5 2D" << std::endl;

        typedef disk::generic_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::fvca5_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.mesh2d$") ))
    {
        std::cout << "Guessed mesh format: Netgen 2D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 2>  mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    if (std::regex_match(filename, std::regex(".*\\.msh$") ))
    {
        std::cout << "Guessed mesh format: FVCA6 3D" << std::endl;

        typedef disk::generic_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::fvca6_mesh_loader<RealType, 3> loader;


        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }

        loader.populate_mesh(msh);

        process_mesh(msh);
    }
*/
    if (std::regex_match(filename, std::regex(".*\\.mesh$") ))
    {
        std::cout << "Guessed mesh format: Netgen 3D" << std::endl;

        typedef disk::simplicial_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::netgen_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }
/*
    if (std::regex_match(filename, std::regex(".*\\.quad$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 2D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 2>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 2> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        auto tr = [](const typename mesh_type::point_type& pt) -> auto {

            auto px = -1 * ( 1-pt.x() ) + 1 * pt.x();
            auto py = -1 * ( 1-pt.y() ) + 1 * pt.y();
            return typename mesh_type::point_type({px, py});
        };

        msh.transform(tr);

        process_mesh(msh);
    }
*/
    if (std::regex_match(filename, std::regex(".*\\.hex$") ))
    {
        std::cout << "Guessed mesh format: Cartesian 3D" << std::endl;

        typedef disk::cartesian_mesh<RealType, 3>   mesh_type;

        mesh_type msh;
        disk::cartesian_mesh_loader<RealType, 3> loader;
        if (!loader.read_mesh(filename))
        {
            std::cout << "Problem loading mesh." << std::endl;
            return 1;
        }
        loader.populate_mesh(msh);

        process_mesh(msh);
    }

    return 0;
}
