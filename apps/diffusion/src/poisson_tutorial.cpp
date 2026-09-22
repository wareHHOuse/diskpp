#include "diskpp/mesh/mesh.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/bases/bases.hpp"
#include "diskpp/bases/bases_new.hpp"
#include "diskpp/methods/hho_assemblers.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/dg"
#include "diskpp/solvers/direct_solvers.hpp"
#include "diskpp/solvers/iterative_solvers.hpp"
#include "diskpp/output/silo.hpp"
#include "diskpp/common/timecounter.hpp"
#include "diskpp/loaders/loader.hpp"

template<typename Mesh>
struct source_functor;

template<disk::mesh_1D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        return M_PI*M_PI*sx;
    }
};

template<disk::mesh_2D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return 2.0*M_PI*M_PI*sx*sy;
    }
};

template<disk::mesh_3D Mesh>
struct source_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return 3.0*M_PI*M_PI*sx*sy*sz;
    }
};


template<typename Mesh>
struct solution_functor;

template<disk::mesh_1D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        return sx;
    }
};

template<disk::mesh_2D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        return sx*sy;
    }
};

template<disk::mesh_3D Mesh>
struct solution_functor<Mesh> {
    using point_type = typename Mesh::point_type;
    auto operator()(const point_type& pt) const {
        auto sx = std::sin(M_PI*pt.x());
        auto sy = std::sin(M_PI*pt.y());
        auto sz = std::sin(M_PI*pt.z());
        return sx*sy*sz;
    }
};

template<typename Mesh>
auto make_rhs_function(const Mesh& msh)
{
    return source_functor<Mesh>();
}

template<typename Mesh>
auto make_solution_function(const Mesh& msh)
{
    return solution_functor<Mesh>();
}

template<typename Mesh>
auto
dg_diffusion_solver(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type Cpen,
    disk::silo_database& silo)
{   
    std::cout << "DG solver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef disk::dynamic_matrix<T>     matrix_type;
    typedef disk::dynamic_vector<T>     vector_type;

    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;
    using point_type = typename Mesh::point_type;

    auto basis_rescaling = disk::basis::rescaling_strategy::inertial;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);

    auto alpha = [](const mesh_type&, const cell_type&) {
        return 1.0;
    };
    
    timecounter tc;
    tc.tic();
    for (int cell_i = 0; cell_i < msh.cells_size(); cell_i++)
    {
        const auto& tcl = msh.cell_at(cell_i);
        auto tbasis = disk::basis::scaled_monomial_basis(msh, tcl, degree, basis_rescaling);

        auto my_h = diameter(msh, tcl);
        auto my_p = degree;
        auto my_alpha = alpha(msh, tcl);

        matrix_type M = integrate(msh, tcl, tbasis, tbasis);
        matrix_type K = integrate(msh, tcl, grad(tbasis), grad(tbasis));
        vector_type loc_rhs = integrate(msh, tcl, f, tbasis);

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {   
            auto n     = normal(msh, tcl, fc);
            auto hF    = measure(msh, fc);
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            if (nv) {
                auto ncl = nv.value();
                auto nbasis = disk::basis::scaled_monomial_basis(msh, ncl, degree, basis_rescaling);
                assert(tbasis.size() == nbasis.size());

                auto other_h = diameter(msh, ncl);
                auto other_p = degree;
                auto other_alpha = alpha(msh, ncl);
                
                auto my_omega = other_alpha / (my_alpha + other_alpha);
                auto other_omega = my_alpha / (my_alpha + other_alpha);
                auto gamma = Cpen * 2.0*(my_alpha*other_alpha) / (my_alpha + other_alpha);
                
                auto eta_l = Cpen * std::max(
                    my_alpha*my_p*(my_p+1)/my_h,
                    other_alpha*other_p*(other_p+1)/other_h
                );

                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());

                Att += + (gamma/hF) * integrate(msh, fc, tbasis, tbasis);
                Att += - (my_alpha * my_omega) * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - (my_alpha * my_omega) * integrate(msh, fc, tbasis, grad(tbasis).dot(n));

                Atn += - (gamma/hF) * integrate(msh, fc, nbasis, tbasis);
                Atn += - (other_alpha * other_omega) * integrate(msh, fc, grad(nbasis).dot(n), tbasis);
                Atn += + (other_alpha * other_omega) * integrate(msh, fc, nbasis, grad(tbasis).dot(n));

                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                auto eta_l = Cpen * my_alpha*my_p*(my_p+1)/my_h;
                
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                Att += + Cpen * (my_alpha/hF) * integrate(msh, fc, tbasis, tbasis);
                Att += - my_alpha * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - my_alpha * integrate(msh, fc, tbasis, grad(tbasis).dot(n));
                assm.assemble(msh, tcl, tcl, Att);
                
                //loc_rhs += + Cpen * (my_alpha/hF) * integrate(msh, fc, g, tbasis);
                //loc_rhs += - my_alpha * integrate(msh, fc, g, grad(tbasis).dot(n));
            }   
        }

        assm.assemble(msh, tcl, K, loc_rhs);
    }

    auto syssz = cbs * msh.cells_size();

    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;

    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    tc.tic();
    disk::solvers::sparse_lu(assm.LHS, assm.RHS, sol);
    std::cout << " Solver time: " << tc.toc() << std::endl;

    auto sol_fun = make_solution_function(msh);

    std::vector<double> u;
    u.reserve(msh.cells_size());

    T err = 0.0;
    size_t cell_i = 0;
    T errint = 0.0;
    tc.tic();

    for (auto& cl : msh)
    {
        auto cb = disk::basis::scaled_monomial_basis(msh, cl, degree, basis_rescaling);
        auto MMe = integrate(msh, cl, cb, cb);
        auto arhs = integrate(msh, cl, sol_fun, cb);

        vector_type asol = MMe.llt().solve(arhs);
        vector_type lsol = sol.segment(cell_i*cb.size(), cb.size());
        vector_type diff = lsol - asol;
        err += diff.dot(MMe*diff);
        u.push_back( lsol(0) );

        auto qps = integrate(msh, cl, 2*degree);
        for (auto& qp : qps) {
            auto val = lsol.dot(cb(qp.point())) - sol_fun(qp.point());
            errint += qp.weight() * val * val;
        }

        cell_i++;
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error (mass): " << std::sqrt(err) << std::endl;;
    std::cout << " L2-norm error (int) : " << std::sqrt(errint) << std::endl;;
    silo.add_variable("mesh", "u_dg", u, disk::zonal_variable_t);
    
    return sol;
}

template<typename Mesh>
void
hho_diffusion_solver(const Mesh& msh, size_t degree, disk::silo_database& silo)
{
    std::cout << "HHO solver" << std::endl;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;

    degree_info di(degree);

    auto f = make_rhs_function(msh);

    auto assm = make_assembler(msh, di);

    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    std::vector<std::pair<MT, VT>> lcs;

    timecounter tc;
    tc.tic();
    for (auto& cl : msh)
    {
        auto [R, A] = local_operator(msh, cl, di);
        auto S = local_stabilization(msh, cl, di, R);
        disk::dynamic_matrix<T> lhs = A+S;
        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        disk::dynamic_vector<T> rhs = integrate(msh, cl, f, phiT);
        lcs.push_back({lhs, rhs});
        auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
        assm.assemble(msh, cl, lhsc, rhsc);
    }
    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;
    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();
    disk::dynamic_vector<T> sol;
    disk::solvers::sparse_lu(assm.LHS, assm.RHS, sol);    
    std::cout << " Solver time: " << tc.toc() << std::endl;

    std::vector<T> u_data;
    
    T error = 0.0;
    T L2error = 0.0;
    auto u_sol = make_solution_function(msh);
    tc.tic();
    size_t cell_i = 0;
    for (auto& cl : msh)
    {
        auto phiT = hho_space<mesh_type>::cell_basis(msh, cl, di.cell);
        auto MMe = integrate(msh, cl, phiT, phiT);
        const auto& [lhs, rhs] = lcs[cell_i++];
        disk::dynamic_vector<T> sol_ana = local_reduction(msh, cl, di, u_sol);
        auto locsolF = assm.take_local_solution(msh, cl, sol);
        disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
        u_data.push_back(locsol(0));
        disk::dynamic_vector<T> diff = locsol - sol_ana;
        error += diff.dot(lhs*diff);
        disk::dynamic_vector<T> diffT = diff.head(phiT.size());
        L2error += diffT.transpose() * (MMe*diffT);
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(L2error) << ", ";
    std::cout << "A-norm error: " << std::sqrt(error) << std::endl;

    silo.add_variable("mesh", "u_hho", u_data, disk::zonal_variable_t);
}

int main(int argc, char **argv)
{
    using T = double;

    size_t      levels = 0;
    size_t      degree = 0;
    char *      mesh_filename = nullptr;

    int ch;
    while ( (ch = getopt(argc, argv, "r:k:m:")) != -1 )
    {
        switch(ch)
        {
            case 'r':
                levels = std::stoull(optarg);
                break;

            case 'k':
                degree = std::stoull(optarg);
                break;

            case 'm':
                mesh_filename = optarg;
                break;

            case '?':
            default:
                std::cout << "Invalid option" << std::endl;
                return 1;
        }
    }

    using mesh_type = disk::cartesian_mesh<T,3>;

    //using mesh_type = disk::generic_mesh<T,2>;

    mesh_type msh;
    auto mesher = disk::make_simple_mesher(msh);

    //disk::gmsh_geometry_loader<mesh_type> loader;
    //loader.read_mesh("triquad.geo2g");
    //loader.populate_mesh(msh);

    std::vector<double> dids;
    size_t i = 0;
    for (i = 0; i < levels; i++) {
        mesher.refine();
    }
        //std::cout << "Diameter: " << disk::average_diameter(msh) << std::endl;
        msh.statistics();
        std::string silo_fn = "poisson_level_" + std::to_string(i) + ".silo";
        disk::silo_database db;
        db.create(silo_fn);
        db.add_mesh(msh, "mesh");
        for (auto& cl : msh) {
            auto di = msh.domain_info(cl);
            dids.push_back( di.tag() );
        }
        db.add_variable("mesh", "domain_ids", dids, disk::zonal_variable_t);
        //dg_diffusion_solver(msh, degree+1, 10.0, db);
        hho_diffusion_solver(msh, degree, db);
}