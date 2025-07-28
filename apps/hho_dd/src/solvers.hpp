
namespace disk {

template<typename Mesh>
struct source_functor;

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
    disk::dynamic_vector<T> sol = mumps_lu(assm.LHS, assm.RHS);
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

    silo.add_variable("srcmesh", "u_hho", u_data, disk::zonal_variable_t);
}

template<typename Mesh>
auto
dg_diffusion_solver(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type eta,
    disk::silo_database& silo)
{   
    std::cout << "DG solver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    auto basis_rescaling = disk::basis::rescaling_strategy::inertial;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);
    
    timecounter tc;
    tc.tic();
    for (auto& tcl : msh)
    {
        auto tbasis = disk::basis::scaled_monomial_basis(msh, tcl, degree, basis_rescaling);
        
        matrix_type K = integrate(msh, tcl, grad(tbasis), grad(tbasis));
        vector_type loc_rhs = integrate(msh, tcl, f, tbasis);

        auto fcs = faces(msh, tcl);
        for (auto& fc : fcs)
        {   
            auto n     = normal(msh, tcl, fc);
            auto eta_l = eta / diameter(msh, fc);
            
            auto nv = cvf.neighbour_via(msh, tcl, fc);
            if (nv) {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                matrix_type Atn = matrix_type::Zero(tbasis.size(), tbasis.size());
                
                auto ncl = nv.value();
                auto nbasis = disk::basis::scaled_monomial_basis(msh, ncl, degree, basis_rescaling);
                assert(tbasis.size() == nbasis.size());

                Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                Att += - 0.5 * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - 0.5 * integrate(msh, fc, tbasis, grad(tbasis).dot(n));

                Atn += - eta_l * integrate(msh, fc, nbasis, tbasis);
                Atn += - 0.5 * integrate(msh, fc, grad(nbasis).dot(n), tbasis);
                Atn += + 0.5 * integrate(msh, fc, nbasis, grad(tbasis).dot(n));

                assm.assemble(msh, tcl, tcl, Att);
                assm.assemble(msh, tcl, ncl, Atn);
            }
            else {
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                Att += + eta_l * integrate(msh, fc, tbasis, tbasis);
                Att += - integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - integrate(msh, fc, tbasis, grad(tbasis).dot(n));
                assm.assemble(msh, tcl, tcl, Att);
            }   
        }
        
        assm.assemble(msh, tcl, K, loc_rhs);
    }

    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;

    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    tc.tic();
    sol = mumps_lu(assm.LHS, assm.RHS);
    std::cout << " Solver time: " << tc.toc() << std::endl;

    auto sol_fun = make_solution_function(msh);

    std::vector<double> u;
    u.reserve(msh.cells_size());

    T err = 0.0; size_t cell_i = 0;
    tc.tic();
    for (auto& cl : msh)
    {
        auto cb = disk::basis::scaled_monomial_basis(msh, cl, degree, basis_rescaling);
        auto MMe = integrate(msh, cl, cb, cb);
        auto arhs = integrate(msh, cl, sol_fun, cb);

        Matrix<T, Dynamic, 1> asol = MMe.llt().solve(arhs);
        Matrix<T, Dynamic, 1> lsol = sol.segment(cell_i*cb.size(), cb.size());
        Matrix<T, Dynamic, 1> diff = lsol - asol;
        err += diff.dot(MMe*diff);
        u.push_back( lsol(0) );
        cell_i++;
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(err) << std::endl;;
    silo.add_variable("dstmesh", "u_dg", u, disk::zonal_variable_t);

    return sol;
}

}
