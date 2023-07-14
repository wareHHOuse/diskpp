/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include "diskpp/quadratures/quadratures.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/mesh/mesh.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/solvers/solver.hpp"
#include "diskpp/output/silo.hpp"

#include "diskpp/methods/hho"

#include "mumps.hpp"

namespace disk {

//template<typename Mesh>
struct source
{
    double operator()(const point<double,2>& pt) const {
        return 2.0*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    }

    double operator()(const point<double,3>& pt) const {
        return 3.0*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
    }
};
/*
template<mesh_2D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
    }
};

template<mesh_3D Mesh>
struct source<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        return M_PI*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y())*std::sin(M_PI*pt.z());
    }
};
*/
}

template<typename Mesh>
void
diffusion_solver(const Mesh& msh)
{
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;

    degree_info di(2);

    disk::hho_degree_info hdi;
    hdi.cell_degree(di.cell);
    hdi.face_degree(di.face);
    hdi.reconstruction_degree(di.reco);

    disk::source f;

    auto assm = make_assembler(msh, di);

    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);

        auto [OR, OA] = disk::make_scalar_hho_laplacian(msh, cl, hdi);
        auto OS = disk::make_scalar_hdg_stabilization(msh, cl, hdi);        

        disk::dynamic_matrix<T> lhs = A+S;
    
        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);

        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);

        assm.assemble(msh, cl, lhsc, rhsc);
    }

    assm.finalize();

    std::cout << assm.LHS.rows() << " " << assm.LHS.nonZeros() << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.dofs());

    disk::solvers::conjugated_gradient_params<T> cgp;
    cgp.verbose = true;
    cgp.max_iter = 0;
    cgp.rr_tol = 1e-10;
    //disk::solvers::conjugated_gradient(cgp, assm.LHS, assm.RHS, sol);

    sol = mumps_lu(assm.LHS, assm.RHS);
    std::vector<T> p;

    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);

        auto [OR, OA] = disk::make_scalar_hho_laplacian(msh, cl, hdi);
        auto OS = disk::make_scalar_hdg_stabilization(msh, cl, hdi);        

        disk::dynamic_matrix<T> lhs = A+S;

        auto phiT = typename hho_space<mesh_type>::cell_basis_type(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
    
        auto locsolF = assm.take_local_solution(msh, cl, sol);

        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        p.push_back(locsol(0));
        //p.push_back( f(barycenter(msh, cl)) );
    }


    disk::silo_database silo_db;
    silo_db.create("diffusion.silo");
    silo_db.add_mesh(msh, "mesh");

    disk::silo_zonal_variable<double> ens("en", p);
    silo_db.add_variable("mesh", ens);

    silo_db.close();
}

int main(void)
{
    using T = double;

    using mesh_type = disk::simplicial_mesh<T,2>;
    mesh_type msh;
    using point_type = typename mesh_type::point_type;
    auto mesher = disk::make_simple_mesher(msh);
    mesher.refine();
    mesher.refine();
    mesher.refine();
    mesher.refine();
    //mesher.refine();
    
    /*
    using mesh_type = disk::generic_mesh<T,2>;
    mesh_type msh;
    auto mesher = disk::make_fvca5_hex_mesher(msh);
    mesher.make_level(5);
    */

    /*
    const char *mesh_filename = ".././../../refactor_old_diskpp_code/meshes/3D_tetras/netgen/cube3.mesh";
    disk::simplicial_mesh<T, 3> msh;
    disk::load_mesh_netgen<T>(mesh_filename, msh);
    */
    diffusion_solver(msh);

    return 0;
}