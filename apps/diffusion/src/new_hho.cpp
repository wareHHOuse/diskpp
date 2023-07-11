#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"

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
    int cell;
    int face;
    int reco;
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
    
    auto scale = diameter(msh, cl);

    size_t offset = szT;
    for (const auto& fc : fcs)
    {
        auto phiF = typename Space::face_basis_type(msh, fc, di.face);
        dynamic_matrix<T> MF = integrate(msh, fc, phiF, phiF);
        dynamic_matrix<T> MTF = integrate(msh, fc, phiT, phiF);
        dynamic_matrix<T> rhs = dynamic_matrix<T>::Zero(szF, sz_total);
        rhs.block(0,0,szF,szT) = -MF.ldlt().solve(MTF);
        rhs.block(0,offset,szF, szF) = dynamic_matrix<T>::Identity(szF, szF);
        S += scale * rhs.transpose() * MF * rhs;
        offset += szF;
    }

    return S;
}

#if 0
template<cell_context Ctx, typename Function>
auto
local_reduction(const Ctx& ctxT, const Function& f)
{
    using namespace disk::basis;
    using T = typename Ctx::coordinate_type;// XXX: fix ScalT
    auto [szT, szF, szR] = space_dimension(ctxT);
    auto num_faces = ctxT.num_faces();

    auto sz_total = szT + num_faces*szF;
    dynamic_vector<T> ret = dynamic_vector<T>::Zero(sz_total);
    auto phiT = scaled_monomial_basis(ctxT, ctxT.degree_);

    ret.head(szT) = L2_project(f, phiT, ctxT);

    size_t offset = szT;
    for (const auto& ctxF : ctxT)
    {
        auto phiF = scaled_monomial_basis(ctxF, ctxT.degree_);
        ret.segment(offset, szF) = L2_project(f, phiF, ctxF);
        offset += szF;
    }

    return ret;
}

#endif

} // namespace slapl
} // namespace hho
} // namespace disk

int main(void)
{
    using T = double;

    disk::simplicial_mesh<T,2> msh;
    auto mesher = disk::make_simple_mesher(msh);
    //mesher.refine();
    //mesher.refine();

    disk::hho::slapl::degree_info di;
    di.cell = 0;
    di.face = 0;
    di.reco = 1;

    for (auto& cl : msh)
    {
        auto [R, A] = disk::hho::slapl::local_operator(msh, cl, di);
        std::cout << R << std::endl;
        auto S = disk::hho::slapl::local_stabilization(msh, cl, di, R);
        std::cout << S << std::endl;
    }

    return 0;
}