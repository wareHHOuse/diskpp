/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#include <functional>
using namespace std::placeholders;

#include "hho.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

#include "common/eigen.hpp"

namespace disk {

template<typename T>
static_matrix<T,2,2>
make_material_tensor(const point<T,2>& pt)
{
    static_matrix<T,2,2> ret = static_matrix<T,2,2>::Identity();
    //return ret;
    //return ret * 6.72071;

    auto c = cos(M_PI * pt.x()/0.02);
    auto s = sin(M_PI * pt.y()/0.02);
    return ret * (1 + 100*c*c*s*s);
}

template<typename T>
size_t
local_face_num(const simplicial_mesh<T, 2>& msh,
               const typename simplicial_mesh<T, 2>::cell& cl,
               const typename simplicial_mesh<T, 2>::face& fc)
{
    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
        if (fcs[i] == fc)
            return i;

    throw std::invalid_argument("Face not part of the specified cell");
}


template<typename Mesh>
class multiscale_local_basis
{
public:
    typedef Mesh                                        outer_mesh_type;
    typedef typename outer_mesh_type::scalar_type       scalar_type;
    typedef typename outer_mesh_type::cell              outer_cell_type;
    typedef typename outer_mesh_type::face              outer_face_type;

    typedef dynamic_matrix<scalar_type>                 matrix_type;
    typedef dynamic_vector<scalar_type>                 vector_type;
    typedef Eigen::SparseMatrix<scalar_type>            sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>                 triplet_type;

    typedef submesher<outer_mesh_type>                  submesher_type;
    typedef typename submesher_type::inner_mesh_type    inner_mesh_type;
    typedef typename inner_mesh_type::cell              inner_cell_type;
    typedef typename inner_mesh_type::face              inner_face_type;
    typedef typename inner_mesh_type::point_type        point_type;

private:

    typedef basis_quadrature_data<inner_mesh_type,
                                  scaled_monomial_scalar_basis,
                                  quadrature>           bqdata_type;

    inner_mesh_type               m_inner_mesh;
    sparse_matrix_type            m_matrix;
    matrix_type                   m_rhs, m_X;

    size_t                        m_degree, m_refinement_levels;

    bqdata_type                   bqd;

    std::vector<matrix_type>    inner_gr_opers, inner_gr_data;

    matrix_type
    take_local_dofs(const inner_cell_type& cl)
    {
        auto num_cell_dofs = bqd.cell_basis.size();
        auto num_face_dofs = bqd.face_basis.size();

        auto cid = find_element_id(m_inner_mesh.cells_begin(), m_inner_mesh.cells_end(), cl);
        if (!cid.first)
            throw std::invalid_argument("This is a bug: cell not found");

        auto cell_offset = cid.second * num_cell_dofs;

        auto fcs = faces(m_inner_mesh, cl);

        auto cell_total_dofs = num_cell_dofs + fcs.size() * num_face_dofs;

        matrix_type locdata = matrix_type::Zero(cell_total_dofs, m_X.cols());

        locdata.block(0,0,num_cell_dofs, m_X.cols()) = m_X.block(cell_offset, 0, num_cell_dofs, m_X.cols());

        size_t fn = 0;
        for (auto& fc : fcs)
        {
            auto fid = find_element_id(m_inner_mesh.faces_begin(), m_inner_mesh.faces_end(), fc);
            if (!fid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_offset = fid.second * num_face_dofs;

            locdata.block(num_cell_dofs + fn*num_face_dofs, 0, num_face_dofs, m_X.cols()) =
                m_X.block(num_cell_dofs * m_inner_mesh.cells_size() + face_offset, 0, num_face_dofs, m_X.cols());

            fn++;
        }

        return locdata;
    }

public:
    template<typename OuterCellBasis, typename OuterFaceBasis>
    multiscale_local_basis(const outer_mesh_type& outer_msh,
                           const outer_cell_type& outer_cl,
                           const OuterCellBasis& outer_cell_basis,
                           const OuterFaceBasis& outer_face_basis,
                           size_t degree, size_t refinement_levels)
        : m_degree(degree), m_refinement_levels(refinement_levels)
    {
        submesher_type                                  submesher;
        bqd = bqdata_type(m_degree, m_degree);
        gradient_reconstruction_bq<bqdata_type>         gradrec(bqd);
        diffusion_like_stabilization_bq<bqdata_type>    stab(bqd);

        m_inner_mesh = submesher.generate_mesh(outer_msh, outer_cl, m_refinement_levels);
        //std::cout << "Inner cells: " << m_inner_mesh.cells_size() << std::endl;
        //std::cout << "Inner h: " << average_diameter(m_inner_mesh) << std::endl;

        auto quadpair = make_quadrature(outer_msh, 2*outer_cell_basis.degree()+6,
                                                    2*outer_face_basis.degree()+6); //XXX

        auto num_cell_faces = howmany_faces(outer_msh, outer_cl);

        auto num_cell_dofs = howmany_dofs(bqd.cell_basis);
        auto num_face_dofs = howmany_dofs(bqd.face_basis);

        auto num_cell_funcs = howmany_dofs(outer_cell_basis);
        auto num_face_funcs = howmany_dofs(outer_face_basis);

        auto matrix_face_offset = num_cell_dofs * m_inner_mesh.cells_size();
        auto matrix_mult_offset = matrix_face_offset + num_face_dofs * m_inner_mesh.faces_size();
        auto system_size = num_cell_dofs * m_inner_mesh.cells_size() +  // HHO cell dofs
                           num_face_dofs * m_inner_mesh.faces_size() +  // HHO face dofs
                           num_face_funcs * num_cell_faces;             // Lagrange multipliers

        m_matrix = sparse_matrix_type(system_size, system_size);

        auto rhs_size = num_cell_funcs + num_face_funcs * num_cell_faces;
        m_rhs = matrix_type::Zero(system_size, rhs_size);

        std::vector<triplet_type> triplets;
        triplets.reserve(100*m_inner_mesh.cells_size());

        /* list of the boundary ids on the fine face. there must be as many
         * boundary ids as the faces in the coarse cell */
        auto bid_list = m_inner_mesh.boundary_id_list();
        assert( bid_list.size() == howmany_faces(outer_msh, outer_cl) );

        /* Assemble standard HHO part */
        //std::cout << "Inner assembly" << std::endl;
        size_t cell_idx = 0;
        for (auto& cl : m_inner_mesh)
        {
            std::vector<size_t> l2g(num_cell_dofs + num_cell_faces * num_face_dofs);

            /* Build DOF offset table: cell */
            for (size_t i = 0; i < num_cell_dofs; i++)
                l2g[i] = cell_idx * num_cell_dofs + i;

            /* Build DOF offset table: faces */
            auto fcs = faces(m_inner_mesh, cl);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto fc = fcs[i];
                auto eid = find_element_id(m_inner_mesh.faces_begin(), m_inner_mesh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                auto face_id = eid.second;
                /* global offset of current face */
                auto face_offset = matrix_face_offset + face_id * num_face_dofs;

                /* offset in the DOF table */
                auto dt_ofs = num_cell_dofs + i * num_face_dofs;

                for (size_t j = 0; j < num_face_dofs; j++)
                    l2g[dt_ofs+j] = face_offset+j;
            }

            auto tf = std::bind(make_material_tensor<scalar_type>, _1);

            /* Compute HHO element contribution */
            gradrec.compute(m_inner_mesh, cl, tf);
            stab.compute(m_inner_mesh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            assert(loc.rows() == l2g.size());
            assert(loc.cols() == l2g.size());

            inner_gr_opers.push_back(gradrec.oper);
            inner_gr_data.push_back(gradrec.data);

            /* Assemble into the matrix */
            for (size_t i = 0; i < l2g.size(); i++)
                for (size_t j = 0; j < l2g.size(); j++)
                    triplets.push_back( triplet_type(l2g[i], l2g[j], loc(i,j)) );

            /* Now compute the multiscale-specific stuff */
            matrix_type cell_rhs;
            cell_rhs = matrix_type::Zero(num_cell_dofs, num_cell_funcs);

            auto cell_quadpoints = bqd.cell_quadrature.integrate(m_inner_mesh, cl);
            //auto cell_quadpoints = quadpair.first.integrate(m_inner_mesh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto phi = bqd.cell_basis.eval_functions(m_inner_mesh, cl, qp.point());
                auto phi_coarse = outer_cell_basis.eval_functions(outer_msh, outer_cl, qp.point());
                cell_rhs += qp.weight() * phi * phi_coarse.transpose();
            }

            auto cell_offset = cell_idx * num_cell_dofs;
            m_rhs.block(cell_offset, 0, num_cell_dofs, num_cell_funcs) += cell_rhs;
            cell_idx++;
        } // for (auto& cl : m_inner_mesh)

        //std::cout << "Boundaries" << std::endl;
        for (auto itor = m_inner_mesh.boundary_faces_begin();
                  itor != m_inner_mesh.boundary_faces_end();
                  itor++)
        {
            auto fc = *itor;
            //std::cout << fc << std::endl;
            auto eid = find_element_id(m_inner_mesh.faces_begin(), m_inner_mesh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            /* fine face number */
            auto fine_face_id = eid.second;

            auto face_offset = matrix_face_offset + fine_face_id * num_face_dofs;

            /* lookup the boundary number of the current face */
            size_t bnd_id = m_inner_mesh.boundary_id(fc);

            /* since the boundary id on the fine mesh is inherited from the face
             * id on the coarse mesh, we use it to recover the coarse face */
            auto coarse_fc = *std::next(outer_msh.faces_begin(), bnd_id);

            size_t coarse_cell_face_num = local_face_num(outer_msh, outer_cl, coarse_fc);
            assert (coarse_cell_face_num < howmany_faces(outer_msh, outer_cl));
            //std::cout << coarse_fc << " " << fc << " (->) " << coarse_cell_face_num << std::endl;

            matrix_type face_matrix, ff;
            face_matrix = matrix_type::Zero(num_face_dofs, num_face_funcs);
            ff = matrix_type::Zero(num_face_funcs, num_face_funcs);

            //std::cout << " INT " << std::endl;
            auto face_quadpoints = bqd.face_quadrature.integrate(m_inner_mesh, fc);
            //auto face_quadpoints = quadpair.second.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                //std::cout << "Quadrature point: " << qp << std::endl;
                auto phi = bqd.face_basis.eval_functions(m_inner_mesh, fc, qp.point());
                auto phi_coarse = outer_face_basis.eval_functions(outer_msh, coarse_fc, qp.point());
                face_matrix += qp.weight() * phi * phi_coarse.transpose();
                //std::cout << "Fine face basis:   " << phi.transpose() << std::endl;
                //std::cout << "Coarse face basis: " << phi_coarse.transpose() << std::endl;
            }

            //std::cout << "Face matrix:" << std::endl;
            //std::cout << face_matrix << std::endl;

            //auto coarse_face_quadpoints = bqd.face_quadrature.integrate(outer_msh, coarse_fc);
            auto coarse_face_quadpoints = quadpair.second.integrate(outer_msh, coarse_fc);

            for (auto& qp : coarse_face_quadpoints)
            {
                //std::cout << "Quadrature point: " << qp << std::endl;
                auto phi_coarse = outer_face_basis.eval_functions(outer_msh, coarse_fc, qp.point());
                ff += qp.weight() * phi_coarse * phi_coarse.transpose();
                //std::cout << "Coarse face basis: " << phi_coarse.transpose() << std::endl;
            }

            auto mult_offset = matrix_mult_offset + coarse_cell_face_num * num_face_funcs;
            auto rhs_offset = num_cell_funcs + coarse_cell_face_num * num_face_funcs;

            for (size_t j = 0; j < face_matrix.rows(); j++)
            {
                for (size_t k = 0; k < face_matrix.cols(); k++)
                {
                    size_t row = face_offset + j;
                    size_t col = mult_offset + k;

                    triplets.push_back( triplet_type(row, col, -face_matrix(j,k)) );
                    triplets.push_back( triplet_type(col, row, face_matrix(j,k)) );
                }
            }

            m_rhs.block(mult_offset, rhs_offset, ff.rows(), ff.cols()) += ff;

        }

        m_matrix.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();

        //std::cout << "solve" << std::endl;


#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

        solver.analyzePattern(m_matrix);
        solver.factorize(m_matrix);
        m_X = solver.solve(m_rhs);

        //std::cout << "solver done" << std::endl;
    }

    inner_mesh_type
    inner_mesh(void)
    {
        return m_inner_mesh;
    }

#define EVAL_WITH_RECONSTRUCTION

    vector_type
    eval_functions(const inner_cell_type& cl, const point_type& pt)
    {
        matrix_type local_dofs = take_local_dofs(cl);

        //gradient_reconstruction_bq<bqdata_type> gradrec(bqd);
        //gradrec.compute(m_inner_mesh, cl);

        vector_type ret = vector_type::Zero( local_dofs.cols() );

        auto tf = std::bind(make_material_tensor<scalar_type>, _1);
        //gradient_reconstruction_bq<bqdata_type> gradrec(bqd);
        //gradrec.compute(m_inner_mesh, cl, tf);

        auto cell_id = m_inner_mesh.lookup(cl);
        auto gr_oper = inner_gr_opers.at(cell_id);
        for (size_t i = 0; i < local_dofs.cols(); i++)
        {
#ifdef EVAL_WITH_RECONSTRUCTION
            vector_type func_dofs = local_dofs.block(0, i, local_dofs.rows(), 1);

            vector_type R = gr_oper * func_dofs;

            vector_type phi = bqd.cell_basis.eval_functions(m_inner_mesh, cl, pt, 0, m_degree+1);

            assert(R.size() == phi.rows()-1);

            for (size_t j = 1; j < phi.rows(); j++)
                ret(i) += R(j-1)*phi(j);

            ret(i) += func_dofs(0)*phi(0);
#else
            auto phi = bqd.cell_basis.eval_functions(m_inner_mesh, cl, pt, 0, m_degree);
            vector_type func_dofs = local_dofs.block(0, i, phi.size(), 1);

            ret(i) = func_dofs.dot(phi);

#endif
        }

        return ret;
    }

    matrix_type
    eval_gradients(const inner_cell_type& cl, const point_type& pt)
    {
        matrix_type local_dofs = take_local_dofs(cl);

        matrix_type ret = matrix_type::Zero( local_dofs.cols(), 2 );

        auto tf = std::bind(make_material_tensor<scalar_type>, _1);
        //gradient_reconstruction_bq<bqdata_type> gradrec(bqd);
        //gradrec.compute(m_inner_mesh, cl, tf);

        auto cell_id = m_inner_mesh.lookup(cl);
        auto gr_oper = inner_gr_opers.at(cell_id);
        for (size_t i = 0; i < local_dofs.cols(); i++)
        {

#ifdef EVAL_WITH_RECONSTRUCTION
            vector_type func_dofs = local_dofs.block(0, i, local_dofs.rows(), 1);

            vector_type R = gr_oper * func_dofs;

            matrix_type dphi = bqd.cell_basis.eval_gradients(m_inner_mesh, cl, pt, 0, m_degree+1);

            assert(R.size() == dphi.rows()-1);

            for (size_t j = 1; j < dphi.rows(); j++)
                ret.block(i,0,1,2) += R(j-1)*dphi.block(j,0,1,2);

            ret.block(i,0,1,2) += func_dofs(0)*dphi.block(0,0,1,2);
#else
            matrix_type dphi = bqd.cell_basis.eval_gradients(m_inner_mesh, cl, pt, 0, m_degree);
            vector_type func_dofs = local_dofs.block(0, i, local_dofs.rows(), 1);
            
            for (size_t j = 0; j < dphi.rows(); j++)
                ret.block(i,0,1,2) += func_dofs(j)*dphi.block(j,0,1,2);
#endif

        }

        return ret;
    }

    matrix_type //vector_type
    get_lagrange_multipliers()
    {
        
        auto num_cell_dofs = howmany_dofs(bqd.cell_basis);
        auto num_face_dofs = howmany_dofs(bqd.face_basis);
        auto begin = num_cell_dofs * m_inner_mesh.cells_size() + num_face_dofs * m_inner_mesh.faces_size();
        auto end = m_X.rows();
        //std::cout << m_X.rows() << " " << m_X.cols() << std::endl;
        //std::cout << begin << " " << end << std::endl;
        matrix_type ret = m_X.block(begin, 0, end-begin, m_X.cols());
        return ret;
    }

    size_t size() const
    {
        return m_X.cols();
    }

    size_t degree() const
    {
        return m_degree;
    }

};

template<typename Mesh, typename CellBasis, typename FaceBasis, typename Function>
dynamic_vector<typename Mesh::scalar_type>
make_rhs(const Mesh& msh, const typename Mesh::cell& cl,
         const CellBasis& cb, const FaceBasis& fb, const Function& f)
{
    typedef dynamic_vector<typename Mesh::scalar_type> vector_type;
    auto fcs = faces(msh, cl);
    //std::sort(fcs.begin(), fcs.end());
    auto num_faces = fcs.size();
    auto rhs_size = cb.size() + num_faces * fb.size();
    vector_type rhs = vector_type::Zero(rhs_size);

    quadrature<Mesh, typename Mesh::cell> cq(2*cb.degree());
    auto cell_quadpoints = cq.integrate(msh, cl);
    auto cbs = cb.size();
    for (auto& qp : cell_quadpoints)
    {
        auto phi    = cb.eval_functions(msh, cl, qp.point());
        auto f_val  = f(qp.point());
        rhs.block(0, 0, cbs, 1) += qp.weight() * f_val * phi;
    }

    auto fbs = fb.size();
    quadrature<Mesh, typename Mesh::face> fq(2*fb.degree());
    for (size_t face_i = 0; face_i < num_faces; face_i++)
    {
        auto fc = fcs[face_i];
        auto face_quadpoints = fq.integrate(msh, fc);
        auto fbs = fb.size();
        for (auto& qp : face_quadpoints)
        {
            auto offset = cbs + face_i*fbs;
            auto phi    = fb.eval_functions(msh, fc, qp.point());
            auto f_val  = f(qp.point());
            rhs.block(offset, 0, fbs, 1) += qp.weight() * f_val * phi;
        }
    }

    return rhs;
}


template<typename Mesh, typename CellBasis, typename FaceBasis>
class projector_new
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;

    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;

    matrix_type                                 mass_matrix;
    Eigen::LLT<matrix_type>                     dec_mass_matrix;

    const mesh_type&                            mesh;
    const cell_type&                            cell;
    const CellBasis&                            cell_basis;
    const FaceBasis&                            face_basis;

public:
    projector_new(const mesh_type& msh, const cell_type& cl,
                  const CellBasis& cb, const FaceBasis& fb)
        : mesh(msh), cell(cl), cell_basis(cb), face_basis(fb)
    {
        auto fcs = faces(msh, cl);
        //std::sort(fcs.begin(), fcs.end());
        auto num_faces = fcs.size();
        auto mass_matrix_size = cell_basis.size() + num_faces * face_basis.size();
        mass_matrix = matrix_type::Zero(mass_matrix_size, mass_matrix_size);

        quadrature<mesh_type, cell_type> cq(2*cell_basis.degree());
        auto cell_quadpoints = cq.integrate(msh, cl);
        auto cbs = cell_basis.size();
        for (auto& qp : cell_quadpoints)
        {
            auto phi = cell_basis.eval_functions(msh, cl, qp.point());
            mass_matrix.block(0, 0, cbs, cbs) +=
                qp.weight() * phi * phi.transpose();
        }

        auto fbs = face_basis.size();
        quadrature<mesh_type, face_type> fq(2*face_basis.degree());
        for (size_t face_i = 0; face_i < num_faces; face_i++)
        {
            auto fc = fcs[face_i];
            auto face_quadpoints = fq.integrate(msh, fc);
            auto fbs = face_basis.size();
            for (auto& qp : face_quadpoints)
            {
                auto offset = cbs + face_i*fbs;
                auto phi = face_basis.eval_functions(msh, fc, qp.point());
                mass_matrix.block(offset, offset, fbs, fbs) +=
                    qp.weight() * phi * phi.transpose();
            }
        }

        dec_mass_matrix = mass_matrix.llt();
    }

    template<typename Function>
    vector_type project(const Function& f)
    {
        auto rhs = make_rhs(mesh, cell, cell_basis, face_basis, f);
        return dec_mass_matrix.solve(rhs);
    }
};



template<typename Mesh, typename CellBasis, typename FaceBasis>
projector_new<Mesh, CellBasis, FaceBasis>
make_projector(const Mesh& msh, const typename Mesh::cell& cl,
               const CellBasis& cb, const FaceBasis& fb)
{
    return projector_new<Mesh, CellBasis, FaceBasis>(msh, cl, cb, fb);
}



template<typename Mesh, typename CellBasis>
typename Mesh::scalar_type
evaluate_dofs(const Mesh& msh, const typename Mesh::cell& cl,
              const CellBasis& cb,
              const dynamic_vector<typename Mesh::scalar_type>& dofs,
              const typename Mesh::point_type& pt)
{
    auto phi = cb.eval_functions(msh, cl, pt);
    
    typename Mesh::scalar_type ret{};

    for (size_t i = 0; i < cb.size(); i++)
        ret += phi(i) * dofs(i);

    return ret;
}

//#define USE_22_FOR_GR

template<typename Mesh, typename OuterCellBasis, typename OuterFaceBasis>
class gradient_reconstruction_multiscale
{
    typedef Mesh                                                mesh_type;
    typedef typename mesh_type::scalar_type                     scalar_type;
    typedef typename mesh_type::cell                            cell_type;
    typedef typename mesh_type::face                            face_type;
    typedef OuterCellBasis                                      cell_basis_type;
    typedef OuterFaceBasis                                      face_basis_type;
    typedef quadrature<mesh_type, cell_type>                    cell_quad_type;
    typedef quadrature<mesh_type, face_type>                    face_quad_type;
    typedef dynamic_matrix<scalar_type>                         matrix_type;
    typedef dynamic_vector<scalar_type>                         vector_type;

    typedef multiscale_local_basis<mesh_type>                   multiscale_basis_type;

    typedef typename multiscale_basis_type::inner_mesh_type     ms_inner_mesh_type;
    typedef typename multiscale_basis_type::inner_cell_type     ms_inner_cell_type;

    multiscale_basis_type&                  multiscale_basis;
    cell_basis_type&                        cell_basis;
    face_basis_type&                        face_basis;

public:
    matrix_type     oper;
    matrix_type     data;
    matrix_type     stiff_mat;

    gradient_reconstruction_multiscale(const mesh_type& outer_mesh,
                                       const typename mesh_type::cell& outer_cl,
                                       multiscale_basis_type& mb,
                                       cell_basis_type& cb,
                                       face_basis_type& fb)
        : multiscale_basis(mb), cell_basis(cb), face_basis(fb)
    {
        auto inner_mesh = multiscale_basis.inner_mesh();
        auto ms_basis_size = multiscale_basis.size();

        /*matrix_type*/ stiff_mat = matrix_type::Zero(ms_basis_size+1, ms_basis_size+1);

        quadrature<ms_inner_mesh_type, ms_inner_cell_type> ms_cell_quad(2*multiscale_basis.degree()+6);
        auto quadpair = make_quadrature(inner_mesh, 2*cb.degree()+6,
                                             2*fb.degree()+6); //XXX

        matrix_type rhs;
        auto rhs_cols_size = cb.size() + howmany_faces(outer_mesh, outer_cl) * fb.size();
        rhs = matrix_type::Zero(ms_basis_size+1, rhs_cols_size);

        auto bid_list = inner_mesh.boundary_id_list();
        assert( bid_list.size() == howmany_faces(outer_mesh, outer_cl) );

#ifdef USE_22_FOR_GR
        matrix_type lmlm = multiscale_basis.get_lagrange_multipliers().transpose();
        //std::cout << "LM" << std::endl;
        //std::cout << lmlm << std::endl;
#endif
        for (auto& icl : inner_mesh)
        {   
            //cell_info(inner_mesh, icl);

            /* stiff_mat will be the lhs */
            auto icl_quadpoints = ms_cell_quad.integrate(inner_mesh, icl);
            for (auto& qp : icl_quadpoints)
            {
                vector_type phi = multiscale_basis.eval_functions(icl, qp.point());
                matrix_type vt_phi = cb.eval_functions(outer_mesh, outer_cl, qp.point());
                matrix_type dphi = multiscale_basis.eval_gradients(icl, qp.point());
                stiff_mat.block(0,0,ms_basis_size,ms_basis_size) +=
                    qp.weight() * dphi * make_material_tensor(qp.point()) * dphi.transpose();

                // Lagrange multiplier to fix the average of the reconstruction
                stiff_mat.block(0, ms_basis_size, ms_basis_size, 1) += qp.weight() * phi;
                stiff_mat.block(ms_basis_size, 0, 1, ms_basis_size) += qp.weight() * phi.transpose();
                rhs.block(ms_basis_size, 0, 1, vt_phi.size()) += qp.weight() * vt_phi.transpose();
            }

            auto rhs_qps = quadpair.first.integrate(inner_mesh, icl);

            for (auto& qp : rhs_qps)
            {
#ifdef USE_22_FOR_GR
                vector_type vt_phi_outer = cb.eval_functions(outer_mesh, outer_cl, qp.point());
                rhs.block(0, 0, cb.size(), cb.size()) +=
                    qp.weight() * vt_phi_outer * vt_phi_outer.transpose();
#else
                matrix_type ms_dphi = multiscale_basis.eval_gradients(icl, qp.point());
                matrix_type vt_dphi = cb.eval_gradients(outer_mesh, outer_cl, qp.point());                
                rhs.block(0, 0, ms_basis_size, cb.size()) +=
                    qp.weight() * (ms_dphi * make_material_tensor(qp.point()).transpose()) * vt_dphi.transpose();
#endif
            }

            auto fcs = faces(inner_mesh, icl);
            auto num_faces = fcs.size();

#ifdef USE_22_FOR_GR
            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto fc = fcs[face_i];
                if ( !inner_mesh.is_boundary(fc) )
                    continue;

                size_t bnd_id = inner_mesh.boundary_id(fc);

                /* since the boundary id on the fine mesh is inherited from the face
                 * id on the coarse mesh, we use it to recover the coarse face */
                auto outer_fc = *std::next(outer_mesh.faces_begin(), bnd_id);
                size_t outer_face_num = local_face_num(outer_mesh, outer_cl, outer_fc);
                assert (outer_face_num < howmany_faces(outer_mesh, outer_cl));
                
                matrix_type face_contrib = matrix_type::Zero(ms_basis_size, fb.size());
                auto face_quadpoints = quadpair.second.integrate(inner_mesh, fc);
                for (auto& qp : face_quadpoints)
                {
                    auto begin = face_i * fb.size();
                    auto lm = lmlm.block(0, begin, lmlm.rows(), fb.size());
                    vector_type vf_phi = fb.eval_functions(outer_mesh, outer_fc, qp.point());
                    face_contrib += (lm*vf_phi) * vf_phi.transpose();
                }
                auto rhs_ofs = cb.size() + face_i*fb.size();
                rhs.block(0, rhs_ofs, ms_basis_size, fb.size()) += face_contrib;
            }
#else
            for (size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto fc = fcs[face_i];
                if ( !inner_mesh.is_boundary(fc) )
                    continue;

                size_t bnd_id = inner_mesh.boundary_id(fc);

                /* since the boundary id on the fine mesh is inherited from the face
                 * id on the coarse mesh, we use it to recover the coarse face */
                auto outer_fc = *std::next(outer_mesh.faces_begin(), bnd_id);
                size_t outer_face_num = local_face_num(outer_mesh, outer_cl, outer_fc);
                assert (outer_face_num < howmany_faces(outer_mesh, outer_cl));
                
                //std::cout << outer_fc << " " << fc << " -> " << outer_face_num << std::endl;

                auto n = normal(outer_mesh, outer_cl, outer_fc);
                auto col_ofs = cb.size() + outer_face_num*fb.size();
                auto face_quadpoints = quadpair.second.integrate(inner_mesh, fc);
                for (auto& qp : face_quadpoints)
                {
                    matrix_type ms_dphi = multiscale_basis.eval_gradients(icl, qp.point());
                    vector_type vf_phi = fb.eval_functions(outer_mesh, outer_fc, qp.point());
                    vector_type vt_phi = cb.eval_functions(outer_mesh, outer_cl, qp.point()); 
                    rhs.block(0, col_ofs, ms_basis_size, fb.size()) += 
                        qp.weight() * ((ms_dphi * make_material_tensor(qp.point()).transpose()) * n) * vf_phi.transpose();
                    rhs.block(0, 0, ms_basis_size, cb.size()) -=
                        qp.weight() * ((ms_dphi * make_material_tensor(qp.point()).transpose()) * n) * vt_phi.transpose();

                }
            }
#endif
        }
        oper = stiff_mat.lu().solve(rhs);
        data = rhs.transpose() * oper;
    }
};


















template<typename Mesh>
class multiscale_local_problem
{
    typedef Mesh                                mesh_type;
    typedef typename mesh_type::scalar_type     scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef typename mesh_type::face            face_type;
    typedef dynamic_matrix<scalar_type>         matrix_type;
    typedef dynamic_vector<scalar_type>         vector_type;
    typedef Eigen::SparseMatrix<scalar_type>    sparse_matrix_type;
    typedef Eigen::Triplet<scalar_type>         triplet_type;

    typedef basis_quadrature_data<mesh_type,
                                  scaled_monomial_scalar_basis,
                                  quadrature>   bqdata_type;

    size_t                        m_degree;
    size_t                        m_refinement_levels;


    vector_type take_local_dofs(const mesh_type& msh, const cell_type& cl,
                                const vector_type& solution, size_t num_cell_dofs,
                                size_t num_face_dofs)
    {
        auto cid = find_element_id(msh.cells_begin(), msh.cells_end(), cl);
        if (!cid.first)
            throw std::invalid_argument("This is a bug: cell not found");

        auto cell_offset = cid.second * num_cell_dofs;

        auto fcs = faces(msh, cl);
        vector_type locdata = vector_type::Zero(num_cell_dofs + fcs.size() * num_face_dofs);

        locdata.block(0, 0, num_cell_dofs, 1) = solution.block(cell_offset, 0, num_cell_dofs, 1);

        size_t fn = 0;
        for (auto& fc : fcs)
        {
            auto fid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!fid.first)
                throw std::invalid_argument("This is a bug: face not found");

            auto face_offset = fid.second * num_face_dofs;

            locdata.block(num_cell_dofs + fn*num_face_dofs, 0, num_face_dofs, 1) =
                solution.block(num_cell_dofs * msh.cells_size() + face_offset, 0, num_face_dofs, 1);

            fn++;
        }

        return locdata;
    }


public:
    sparse_matrix_type                          matrix;
    matrix_type                                 rhs;

    multiscale_local_problem()
    {
        m_degree                = 1;
        m_refinement_levels     = 4;
    }

    multiscale_local_problem(size_t degree, size_t refinement_levels)
    {
        m_degree                = degree;
        m_refinement_levels     = refinement_levels;
    }



    template<typename OuterCellBasis, typename OuterFaceBasis>
    void assemble(const mesh_type& coarse_msh,
                  const cell_type& coarse_cl,
                  const OuterCellBasis& coarse_cell_basis,
                  const OuterFaceBasis& coarse_face_basis)
    {
        submesher<Mesh>                                 submesher;
        bqdata_type                                     bqd(m_degree, m_degree);
        gradient_reconstruction_bq<bqdata_type>         gradrec(bqd);
        diffusion_like_stabilization_bq<bqdata_type>    stab(bqd);

        auto msh = submesher.generate_mesh(coarse_msh, coarse_cl, m_refinement_levels);

        auto quadpair = make_quadrature(msh, m_degree+coarse_cell_basis.degree(),
                                             m_degree+coarse_face_basis.degree());

        auto num_cell_faces = howmany_faces(coarse_msh, coarse_cl);

        auto num_cell_dofs = howmany_dofs(bqd.cell_basis);
        auto num_face_dofs = howmany_dofs(bqd.face_basis);

        auto num_cell_funcs = howmany_dofs(coarse_cell_basis);
        auto num_face_funcs = howmany_dofs(coarse_face_basis);

        auto matrix_face_offset = num_cell_dofs * msh.cells_size();
        auto matrix_mult_offset = matrix_face_offset + num_face_dofs * msh.faces_size();
        auto system_size = num_cell_dofs * msh.cells_size() +
                           num_face_dofs * msh.faces_size() +
                           num_face_funcs * num_cell_faces;

        matrix = sparse_matrix_type(system_size, system_size);

        auto rhs_size = num_cell_funcs + num_face_funcs * num_cell_faces;
        rhs = matrix_type::Zero(system_size, rhs_size);

        std::vector<triplet_type> triplets;

        /* list of the boundary ids on the fine face. there must be as many
         * boundary ids as the faces in the coarse cell */
        auto bid_list = msh.boundary_id_list();
        assert( bid_list.size() == howmany_faces(coarse_msh, coarse_cl) );

        /* Assemble standard HHO part */

        size_t cell_idx = 0;
        for (auto& cl : msh)
        {
            std::vector<size_t> l2g(num_cell_dofs + num_cell_faces * num_face_dofs);

            /* Build DOF offset table: cell */
            for (size_t i = 0; i < num_cell_dofs; i++)
                l2g[i] = cell_idx * num_cell_dofs + i;

            /* Build DOF offset table: faces */
            auto fcs = faces(msh, cl);
            for (size_t i = 0; i < fcs.size(); i++)
            {
                auto fc = fcs[i];
                auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
                if (!eid.first)
                    throw std::invalid_argument("This is a bug: face not found");

                auto face_id = eid.second;
                /* global offset of current face */
                auto face_offset = matrix_face_offset + face_id * num_face_dofs;

                /* offset in the DOF table */
                auto dt_ofs = num_cell_dofs + i * num_face_dofs;

                for (size_t j = 0; j < num_face_dofs; j++)
                    l2g[dt_ofs+j] = face_offset+j;
            }

            /* Compute HHO element contribution */
            gradrec.compute(msh, cl);
            stab.compute(msh, cl, gradrec.oper);
            dynamic_matrix<scalar_type> loc = gradrec.data + stab.data;
            assert(loc.rows() == l2g.size());
            assert(loc.cols() == l2g.size());

            /* Assemble into the matrix */
            for (size_t i = 0; i < l2g.size(); i++)
                for (size_t j = 0; j < l2g.size(); j++)
                    triplets.push_back( triplet_type(l2g[i], l2g[j], loc(i,j)) );

            /* Now compute the multiscale-specific stuff */
            matrix_type cell_rhs;
            cell_rhs = matrix_type::Zero(num_cell_dofs, num_cell_funcs);

            //auto cell_quadpoints = bqd.cell_quadrature.integrate(msh, cl);
            auto cell_quadpoints = quadpair.first.integrate(msh, cl);
            for (auto& qp : cell_quadpoints)
            {
                auto phi = bqd.cell_basis.eval_functions(msh, cl, qp.point());
                auto phi_coarse = coarse_cell_basis.eval_functions(coarse_msh, coarse_cl, qp.point());
                cell_rhs += qp.weight() * phi * phi_coarse.transpose();
            }

            auto cell_offset = cell_idx * num_cell_dofs;
            rhs.block(cell_offset, 0, num_cell_dofs, num_cell_funcs) += cell_rhs;

            cell_idx++;
        }



        for (auto itor = msh.boundary_faces_begin(); itor != msh.boundary_faces_end(); itor++)
        {
            auto fc = *itor;
            //std::cout << fc << std::endl;
            auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
            if (!eid.first)
                throw std::invalid_argument("This is a bug: face not found");

            /* fine face number */
            auto fine_face_id = eid.second;

            auto face_offset = matrix_face_offset + fine_face_id * num_face_dofs;

            /* lookup the boundary number of the current face */
            size_t bnd_id = msh.boundary_id(fc);

            size_t coarse_cell_face_num;
            for (size_t i = 0; i < bid_list.size(); i++)
                if (bid_list[i] == bnd_id)
                {
                    coarse_cell_face_num = i;
                    break;
                }

            assert (coarse_cell_face_num < howmany_faces(coarse_msh, coarse_cl));

            /* since the boundary id on the fine mesh is inherited from the face
             * id on the coarse mesh, we use it to recover the coarse face */
            auto coarse_fc = *(coarse_msh.faces_begin() + bnd_id);

            matrix_type face_matrix, ff;
            face_matrix = matrix_type::Zero(num_face_dofs, num_face_funcs);
            ff = matrix_type::Zero(num_face_funcs, num_face_funcs);

            //std::cout << " INT " << std::endl;
            //auto face_quadpoints = bqd.face_quadrature.integrate(msh, fc);
            auto face_quadpoints = quadpair.second.integrate(msh, fc);
            for (auto& qp : face_quadpoints)
            {
                //std::cout << "Quadrature point: " << qp << std::endl;
                auto phi = bqd.face_basis.eval_functions(msh, fc, qp.point());
                auto phi_coarse = coarse_face_basis.eval_functions(coarse_msh, coarse_fc, qp.point());
                face_matrix += qp.weight() * phi * phi_coarse.transpose();
                //std::cout << "Fine face basis:   " << phi.transpose() << std::endl;
                //std::cout << "Coarse face basis: " << phi_coarse.transpose() << std::endl;
            }

            //std::cout << "Face matrix:" << std::endl;
            //std::cout << face_matrix << std::endl;

            auto coarse_face_quadpoints = bqd.face_quadrature.integrate(coarse_msh, coarse_fc);
            //auto coarse_face_quadpoints = quadpair.second.integrate(coarse_msh, coarse_fc);
            for (auto& qp : coarse_face_quadpoints)
            {
                //std::cout << "Quadrature point: " << qp << std::endl;
                auto phi_coarse = coarse_face_basis.eval_functions(coarse_msh, coarse_fc, qp.point());
                ff += qp.weight() * phi_coarse * phi_coarse.transpose();
                //std::cout << "Coarse face basis: " << phi_coarse.transpose() << std::endl;
            }

            auto mult_offset = matrix_mult_offset + coarse_cell_face_num * num_face_funcs;
            auto rhs_offset = num_cell_funcs + coarse_cell_face_num * num_face_funcs;

            for (size_t j = 0; j < face_matrix.rows(); j++)
            {
                for (size_t k = 0; k < face_matrix.cols(); k++)
                {
                    size_t row = face_offset + j;
                    size_t col = mult_offset + k;

                    triplets.push_back( triplet_type(row, col, -face_matrix(j,k)) );
                    triplets.push_back( triplet_type(col, row, face_matrix(j,k)) );
                }
            }

            rhs.block(mult_offset, rhs_offset, ff.rows(), ff.cols()) += ff;

        }


        std::ofstream ofsr("rhs.dat");

        for (size_t i = 0; i < rhs.rows(); i++)
        {
            for (size_t j = 0; j < rhs.cols(); j++)
                ofsr << rhs(i,j) << " ";
            ofsr << std::endl;
        }

        ofsr.close();


        matrix.setFromTriplets(triplets.begin(), triplets.end());

/*
        std::stringstream ssw;
        ssw << "matrix.mat";
        std::ofstream ofs(ssw.str());

        for (int k=0; k<matrix.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(matrix,k); it; ++it)
                ofs << it.row() << " " << it.col() << " " << it.value() << std::endl;

        ofs.close();
*/
        triplets.clear();



#ifdef HAVE_INTEL_MKL
        Eigen::PardisoLU<Eigen::SparseMatrix<scalar_type>>  solver;
#else
        Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>>   solver;
#endif

        solver.analyzePattern(matrix);
        solver.factorize(matrix);
        matrix_type X = solver.solve(rhs);

        auto eid = find_element_id(coarse_msh.cells_begin(), coarse_msh.cells_end(), coarse_cl);
        if (!eid.first)
            throw std::invalid_argument("This is a bug: cell not found");

        auto cell_id = eid.second;

        for (size_t funcnum = 0; funcnum < num_cell_funcs + 3*num_face_funcs; funcnum++)
        {
            std::stringstream ss;
            ss << "mshho_ko_" << coarse_cell_basis.degree()+1 << "_ki_" << m_degree << "_rl_" << m_refinement_levels << "_fn_" << funcnum << ".dat";
            //ss << "mshho_fn_" << funcnum << ".dat";
            std::ofstream ofs_sol(ss.str());
            cell_idx = 0;
            for (auto& cl : msh)
            {
                gradrec.compute(msh, cl);
                auto cell_sol = X.block(num_cell_dofs * cell_idx, funcnum, num_cell_dofs, 1);

                vector_type locsol = X.block(0, funcnum, X.rows(), 1);
                vector_type locdata = take_local_dofs(msh, cl, locsol, num_cell_dofs, num_face_dofs );

                vector_type R = gradrec.oper * locdata;

                //auto qps = bqd.cell_quadrature.integrate(msh, cl);
                auto qps = make_test_points(msh, cl, 15);
                for (auto& qp : qps)
                {
                    auto phi = bqd.cell_basis.eval_functions(msh, cl, qp, 0, m_degree+1);

                    scalar_type pot = 0.0;
                    for (size_t i = 1; i < phi.size(); i++)
                        pot += phi[i] * R(i-1);

                    pot += locdata(0);

                    auto tp = qp;
                    for (size_t i = 0; i < mesh_type::dimension; i++)
                        ofs_sol << tp[i] << " ";
                    ofs_sol << pot << std::endl;
                }

                cell_idx++;
            }

            ofs_sol.close();
        }
    }
};


} //namespace disk
