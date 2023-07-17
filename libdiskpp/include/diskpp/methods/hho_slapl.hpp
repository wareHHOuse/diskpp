/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#pragma once

#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/bases/bases_traits.hpp"
#include "hho_assemblers.hpp"

namespace disk {
namespace hho {
namespace slapl {

template<typename Mesh, typename ScalT = typename Mesh::coordinate_type>
struct hho_space
{
    using mesh_type = Mesh;
    using scalar_type = ScalT;
    using cell_type = typename Mesh::cell_type;
    using face_type = typename Mesh::face_type;
    using cell_basis_type = disk::basis::scalar_monomial<mesh_type, cell_type, scalar_type>;
    using face_basis_type = disk::basis::scalar_monomial<mesh_type, face_type, scalar_type>;
    using reco_basis_type = disk::basis::scalar_monomial<mesh_type, cell_type, scalar_type>;
};

struct degree_info {
    size_t cell;
    size_t face;
    size_t reco;

    degree_info()
        : cell(0), face(0), reco(1)
    {}

    degree_info(size_t deg)
        : cell(deg), face(deg), reco(deg+1)
    {}

    degree_info(size_t cdeg, size_t fdeg)
        : cell(cdeg), face(fdeg), reco(fdeg+1)
    {}
};

template<typename Space>
auto
space_dimensions(const degree_info& di)
{
    auto szT = Space::cell_basis_type::size_of_degree(di.cell);
    auto szF = Space::face_basis_type::size_of_degree(di.face);
    auto szR = Space::cell_basis_type::size_of_degree(di.reco);

    return std::tuple(szT, szF, szR);
}

template<typename Mesh, typename Space = hho_space<Mesh>>
auto
local_operator(const Mesh& msh, const typename Mesh::cell_type& cl,
    const degree_info& di)
{
    using namespace disk::basis;
    using T = typename Space::scalar_type;

    auto phiR = typename Space::cell_basis_type(msh, cl, di.reco);
    auto phiT = typename Space::cell_basis_type(msh, cl, di.cell);

    auto [szT, szF, szR] = space_dimensions<Space>(di);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();
    assert(szR > 0);
    auto rows = szR - 1;
    auto cols = szT + num_faces*szF;

    dynamic_matrix<T> stiffness = integrate(msh, cl, grad(phiR), grad(phiR));
    dynamic_matrix<T> lhs = stiffness.bottomRightCorner(szR-1, szR-1);
    dynamic_matrix<T> rhs = dynamic_matrix<T>::Zero(rows, cols);
    rhs.block(0,0,szR-1,szT) = stiffness.bottomLeftCorner(szR-1, szT);

    size_t offset = szT;
    for (const auto& fc : fcs)
    {
        auto n = normal(msh, cl, fc);
        auto phiF = typename Space::face_basis_type(msh, fc, di.face);

        rhs.block(0,offset,szR-1,szF) += 
            integrate(msh, fc, phiF, grad(phiR).dot(n)).block(1,0,szR-1,szF);

        rhs.block(0,0,szR-1,szT) -=
            integrate(msh, fc, phiT, grad(phiR).dot(n)).block(1,0,szR-1,szT);

        offset += szF;
    }

    dynamic_matrix<T> R = lhs.ldlt().solve(rhs);
    dynamic_matrix<T> A = rhs.transpose()*R;
    return std::pair(R, A);
}


template<typename Mesh, typename Space = hho_space<Mesh>>
auto
local_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl,
    const degree_info& di, const dynamic_matrix<typename Space::scalar_type>& R)
{
    using namespace disk::basis;
    using T = typename Space::scalar_type;

    auto [szT, szF, szR] = space_dimensions<Space>(di);
    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    auto sz_total = szT + num_faces*szF;
    dynamic_matrix<T> S = dynamic_matrix<T>::Zero(sz_total, sz_total);

    auto phiT = typename Space::cell_basis_type(msh, cl, di.cell);
    auto phiR = typename Space::reco_basis_type(msh, cl, di.reco);
    
    auto scale = 1./diameter(msh, cl);

    
    dynamic_matrix<T> MT = integrate(msh, cl, phiT, phiT);
    dynamic_matrix<T> R2T =
        integrate(msh, cl, phiR, phiT).block(0,1,phiT.size(), phiR.size()-1);
    dynamic_matrix<T> P1 = -MT.ldlt().solve(R2T*R);
    P1.block(0, 0, phiT.size(), phiT.size()) +=
        dynamic_matrix<T>::Identity(phiT.size(), phiT.size());

    size_t offset = szT;
    for (const auto& fc : fcs)
    {
        auto phiF = typename Space::face_basis_type(msh, fc, di.face);
        dynamic_matrix<T> MF = integrate(msh, fc, phiF, phiF);
        Eigen::LDLT<dynamic_matrix<T>> MF_llt(MF);
        dynamic_matrix<T> T2F = integrate(msh, fc, phiT, phiF);
        dynamic_matrix<T> P2 = MF_llt.solve(T2F*P1);

        dynamic_matrix<T> R2F =
            integrate(msh, fc, phiR, phiF).block(0, 1, phiF.size(), phiR.size()-1);
        dynamic_matrix<T> P3 = MF_llt.solve(R2F*R);
        P3.block(0,offset,phiF.size(), phiF.size()) -=
            dynamic_matrix<T>::Identity(phiF.size(), phiF.size());

        dynamic_matrix<T> rhs = P2+P3;

        S += scale * rhs.transpose() * MF * rhs;
        offset += szF;
    }

    return S;
}

template<typename Mesh, typename Space = hho_space<Mesh>, typename Function>
auto
local_reduction(const Mesh& msh, const typename Mesh::cell_type& cl,
    const degree_info& di, const Function& f)
{
    using namespace disk::basis;
    using T = typename Space::scalar_type;
    auto [szT, szF, szR] = space_dimensions<Space>(di);
    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    auto sz_total = szT + num_faces*szF;
    dynamic_vector<T> ret = dynamic_vector<T>::Zero(sz_total);
    auto phiT = typename Space::cell_basis_type(msh, cl, di.cell);

    ret.head(szT) = L2_project(msh, cl, f, phiT);

    size_t offset = szT;
    for (const auto& fc : fcs)
    {
        auto phiF = typename Space::face_basis_type(msh, fc, di.face);
        ret.segment(offset, szF) = L2_project(msh, fc, f, phiF);
        offset += szF;
    }

    return ret;
}

template<typename Mesh, typename Space = hho_space<Mesh>>
auto
make_assembler(const Mesh& msh, const degree_info& di)
{
    using fbt = typename Space::face_basis_type;
    return basic_condensed_assembler<Mesh, fbt>(msh, di.face);
}

} // namespace slapl
} // namespace hho
} // namespace disk
