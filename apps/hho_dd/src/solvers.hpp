#include "diskpp/methods/hho_assemblers.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/solvers/direct_solvers.hpp"
#include "diskpp/solvers/iterative_solvers.hpp"
#include "diskpp/solvers/eigensolvers.hpp"
namespace disk {

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

template<typename T>
T cond(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, int ex = 0)
{
    using mtype = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    Eigen::JacobiSVD<mtype> svd(A);
    return svd.singularValues()(0) 
        / svd.singularValues()(svd.singularValues().size()-(1+ex));
}


template<typename Mesh>
class hho_diffusion_solver
{   
    const Mesh&     msh_;
    disk::hho::slapl::degree_info     di_;
    
    using space_type = disk::hho::slapl::hho_space<Mesh>;

    using mesh_type = Mesh;
    using T = typename space_type::scalar_type;
    using face_basis_type = typename space_type::face_basis_type;
    
    using MT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using VT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    
    std::vector<std::pair<MT, VT>> lcs;

    disk::hho::basic_condensed_assembler<Mesh, face_basis_type> assm;

public:
    hho_diffusion_solver(const Mesh& msh, size_t degree = 1)
        : msh_(msh), di_( disk::hho::slapl::degree_info(degree) ),
          assm(msh, degree)
    {
    }

    template<typename RHSFun>
    dynamic_vector<T> solve(RHSFun& f) {
        
        //timecounter tc;
        //tc.tic();
        for (auto& cl : msh_)
        {
            auto [R, A] = local_operator(msh_, cl, di_);
            auto S = local_stabilization(msh_, cl, di_, R);
            disk::dynamic_matrix<T> lhs = A+S;
            auto phiT = space_type::cell_basis(msh_, cl, di_.cell);
            disk::dynamic_vector<T> rhs = integrate(msh_, cl, f, phiT);
            lcs.push_back({lhs, rhs});
            auto [lhsc, rhsc] = disk::hho::schur(lhs, rhs, phiT);
            assm.assemble(msh_, cl, lhsc, rhsc);;
        }
        assm.finalize();
        //std::cout << " Assembly time: " << tc.toc() << std::endl;
        //std::cout << " Unknowns: " << assm.LHS.rows() << " ";
        //std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
        //tc.tic();
        disk::dynamic_vector<T> sol;
        disk::solvers::sparse_lu(assm.LHS, assm.RHS, sol);
        return sol;
    }

    template<typename RefSol>
    std::pair<T,T> compute_errors(const dynamic_vector<T>& sol, RefSol u_sol)
    {
        T Aerror = 0.0;
        T L2error = 0.0;
    
        size_t cell_i = 0;
        for (auto& cl : msh_)
        {
            auto phiT = space_type::cell_basis(msh_, cl, di_.cell);
            auto MMe = integrate(msh_, cl, phiT, phiT);
            const auto& [lhs, rhs] = lcs[cell_i++];
            disk::dynamic_vector<T> sol_ana = local_reduction(msh_, cl, di_, u_sol);
            auto locsolF = assm.take_local_solution(msh_, cl, sol);
            disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            disk::dynamic_vector<T> diff = locsol - sol_ana;
            Aerror += diff.dot(lhs*diff);
            disk::dynamic_vector<T> diffT = diff.head(phiT.size());
            L2error += diffT.transpose() * (MMe*diffT);
        }

        return {L2error, Aerror};
    }

    dynamic_vector<T>
    compute_nodal_values(const dynamic_vector<T>& sol)
    {
        std::vector<std::pair<T, int>> nodedata(msh_.points_size(), {0.0, 0});
        dynamic_vector<T> vals = dynamic_vector<T>::Zero(msh_.points_size());

        for (size_t cell_i = 0; cell_i < msh_.cells_size(); cell_i++)
        {
            const auto& cl =  msh_.cell_at(cell_i);
            auto phiT = space_type::cell_basis(msh_, cl, di_.cell);
            const auto& [lhs, rhs] = lcs[cell_i];
            auto locsolF = assm.take_local_solution(msh_, cl, sol);
            disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);

            auto ptids = cl.point_ids();
            for(auto& ptid : ptids) {
                auto pt = msh_.point_at(ptid);
                auto val = locsol.head(phiT.size()).dot(phiT(pt));
                nodedata[ptid].first += val;
                nodedata[ptid].second++;
            }
        }

        for (size_t node_i = 0; node_i < msh_.points_size(); node_i++) {
            auto& [val, nodes] = nodedata[node_i];
            vals(node_i) = val/nodes;
        }

        return vals;
    }

    dynamic_vector<T>
    compute_zonal_values(const dynamic_vector<T>& sol)
    {
        dynamic_vector<T> vals = dynamic_vector<T>::Zero(msh_.cells_size());

        for (size_t cell_i = 0; cell_i < msh_.cells_size(); cell_i++)
        {
            const auto& cl =  msh_.cell_at(cell_i);
            auto phiT = space_type::cell_basis(msh_, cl, di_.cell);
            const auto& [lhs, rhs] = lcs[cell_i];
            auto locsolF = assm.take_local_solution(msh_, cl, sol);
            disk::dynamic_vector<T> locsol = disk::hho::deschur(lhs, rhs, locsolF, phiT);
            vals(cell_i) = locsol;
        }

        return vals;
    }
};

#if 0
template<typename Mesh, typename RHSFun>
dynamic_vector<typename Mesh::scalar_type>
xx(const Mesh& msh, size_t degree, RHSFun f)
{
    std::cout << "HHO solver" << std::endl;
    using namespace disk::basis;
    using namespace disk::hho::slapl;

    using mesh_type = Mesh;
    using T = typename hho_space<Mesh>::scalar_type;

    degree_info di(degree);

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
        assm.assemble(msh, cl, lhsc, rhsc);;
    }
    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;
    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;
    tc.tic();
    disk::dynamic_vector<T> sol = mumps_lu(assm.LHS, assm.RHS);
    
    //disk::solvers::conjugated_gradient_params<T> cgp;
    //cgp.verbose = true;
    //cgp.max_iter = 3*assm.RHS.rows();
    //disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.RHS.rows());
    //disk::solvers::conjugated_gradient(cgp, assm.LHS, assm.RHS, sol);
    
    std::cout << " Solver time: " << tc.toc() << std::endl;

    std::vector<std::pair<T, int>> nodedata(msh.points_size(), {0.0, 0});
    std::vector<T> u_data;
    std::vector<T> u_data_nodal;

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
        vcondm.push_back(cond(MMe));

        auto ptids = cl.point_ids();
        for(auto& ptid : ptids) {
            auto pt = msh.point_at(ptid);
            auto val = locsol.head(phiT.size()).dot(phiT(pt));
            nodedata[ptid].first += val;
            nodedata[ptid].second++;
        }
    }

    for (auto& [val, nodes] : nodedata) {
        u_data_nodal.push_back(val/nodes);
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error: " << std::sqrt(L2error) << ", ";
    std::cout << "A-norm error: " << std::sqrt(error) << std::endl;


    silo.add_variable(mname, mname+"_u_hho", u_data_nodal, disk::nodal_variable_t);
    //silo.add_variable("srcmesh", "cond_hho", vcond, disk::zonal_variable_t);
    //silo.add_variable("srcmesh", "cond_mass", vcondm, disk::zonal_variable_t);
}
#endif

template<disk::mesh_2D Mesh>
struct rfun {
    using point_type = typename Mesh::point_type;
    using scalar_type = typename Mesh::coordinate_type;
    
    scalar_type A_, B_, a_, kin_, kout_;

    rfun(scalar_type a, scalar_type kin, scalar_type kout)
        : A_(2*kout/(kin + kout)),
          B_(a*a*(kout - kin)/(kin + kout)),
          a_(a), kin_(kin), kout_(kout)
    {}

    auto kappa(const point_type& pt) const {
        auto r = std::hypot(pt.x(), pt.y());
        if (r < a_)
            return kin_;
        return kout_;
    }
    
    auto operator()(const point_type& pt) const {
        auto r = std::hypot(pt.x(), pt.y());
        if (r < a_) {
            return A_ * pt.x();
        }
        
        return pt.x() * (1.0 + B_/(r*r));
    }
};


template<typename Mesh>
auto
dg_diffusion_solver(Mesh& msh, size_t degree,
    const typename Mesh::coordinate_type Cpen,
    disk::silo_database& silo)
{   
    std::cout << "DG solver" << std::endl;
    auto cvf = connectivity_via_faces(msh);
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    using mesh_type = Mesh;
    using cell_type = typename Mesh::cell_type;

    auto basis_rescaling = disk::basis::rescaling_strategy::inertial;

    auto f = make_rhs_function(msh);

    auto cbs = disk::scalar_basis_size(degree, Mesh::dimension);
    auto assm = make_discontinuous_galerkin_assembler(msh, cbs);
    
    using tripv_t = std::vector<Eigen::Triplet<T>>;
    tripv_t precond_triplets;
    tripv_t mass_triplets;

    auto assm_blockdiag = [&cbs](const matrix_type& m, size_t ofs, tripv_t& triplets) {
        auto c_ofs = ofs*cbs;
        for (int i = 0; i < m.rows(); i++)
            for (int j = 0; j < m.cols(); j++)
                triplets.push_back( {int(c_ofs+i), int(c_ofs+j), m(i,j)} );
    };

    rfun<Mesh> g(0.5, 1e6, 1.0);

    timecounter tc;
    tc.tic();
    for (int cell_i = 0; cell_i < msh.cells_size(); cell_i++)
    {
        const auto& tcl = msh.cell_at(cell_i);
        auto tbasis = disk::basis::scaled_monomial_basis(msh, tcl, degree, basis_rescaling);
        
        auto alpha = [&](const mesh_type& msh, const cell_type& cl) {
            auto bar = barycenter(msh, cl);
            //if (bar.x() > 0.5) {
            //    return 10.0;
            //}
            //return 1.0;
            return g.kappa(bar);
        };

        auto my_h = diameter(msh, tcl);
        auto my_p = degree;
        auto my_alpha = alpha(msh, tcl);

        matrix_type M = integrate(msh, tcl, tbasis, tbasis);
        matrix_type K = my_alpha * integrate(msh, tcl, grad(tbasis), grad(tbasis));
        vector_type loc_rhs = vector_type::Zero(tbasis.size());// integrate(msh, tcl, f, tbasis);
        assm_blockdiag(K, cell_i, precond_triplets);

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

                assm_blockdiag(Att, cell_i, precond_triplets);
            }
            else {
                auto eta_l = Cpen * my_alpha*my_p*(my_p+1)/my_h;
                
                matrix_type Att = matrix_type::Zero(tbasis.size(), tbasis.size());
                Att += + Cpen * (my_alpha/hF) * integrate(msh, fc, tbasis, tbasis);
                Att += - my_alpha * integrate(msh, fc, grad(tbasis).dot(n), tbasis);
                Att += - my_alpha * integrate(msh, fc, tbasis, grad(tbasis).dot(n));
                assm.assemble(msh, tcl, tcl, Att);
                
                loc_rhs += + Cpen * (my_alpha/hF) * integrate(msh, fc, g, tbasis);
                loc_rhs += - my_alpha * integrate(msh, fc, g, grad(tbasis).dot(n));

                assm_blockdiag(Att, cell_i, precond_triplets);
            }   
        }

        assm_blockdiag(M, cell_i, mass_triplets);
        assm.assemble(msh, tcl, K, loc_rhs);
    }

    auto syssz = cbs * msh.cells_size();
    Eigen::SparseMatrix<T> precond = Eigen::SparseMatrix<T>(syssz, syssz);
    precond.setFromTriplets( precond_triplets.begin(), precond_triplets.end() );

    Eigen::SparseMatrix<T> mass = Eigen::SparseMatrix<T>(syssz, syssz);
    mass.setFromTriplets( mass_triplets.begin(), mass_triplets.end() );

    assm.finalize();
    std::cout << " Assembly time: " << tc.toc() << std::endl;

    std::cout << " Unknowns: " << assm.LHS.rows() << " ";
    std::cout << " Nonzeros: " << assm.LHS.nonZeros() << std::endl;

    disk::dynamic_vector<T> sol = disk::dynamic_vector<T>::Zero(assm.syssz);

    tc.tic();
    //sol = mumps_lu(assm.LHS, assm.RHS);
    
    Eigen::SparseLU<Eigen::SparseMatrix<T>> lu(precond);
    auto precondOper = [&](const dynamic_vector<T>& v) {
        return lu.solve(v);
    };

    disk::solvers::iterative_solver_params cgp;
    cgp.verbose = true;
    cgp.max_residual = 1000;
    cgp.max_iter = 3*assm.RHS.rows();
    disk::solvers::cg(cgp, assm.LHS, assm.RHS, sol, precondOper);

    std::cout << " Solver time: " << tc.toc() << std::endl;

    auto [evmax, lambdamax] = disk::solvers::powiter(assm.LHS, 1e-5).value();
    std::cout << "lambda max = " << lambdamax << std::endl;

    auto [evmin, lambdamin] = disk::solvers::inv_powiter(assm.LHS, 0.0, 1e-5).value();
    std::cout << "lambda min = " << lambdamin << std::endl;

    std::cout << "Estimated condition number = " << lambdamax/lambdamin << std::endl;

    auto sol_fun = make_solution_function(msh);

    std::vector<double> u, evmax_data, evmin_data;
    u.reserve(msh.cells_size());
    evmax_data.reserve(msh.cells_size());
    evmin_data.reserve(msh.cells_size());

    T err = 0.0;
    size_t cell_i = 0;
    T errint = 0.0;
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

        auto qps = integrate(msh, cl, 2*degree);
        for (auto& qp : qps) {
            auto val = lsol.dot(cb(qp.point())) - g(qp.point());
            errint += qp.weight() * val * val;
        }

        Matrix<T, Dynamic, 1> levmax = evmax.segment(cell_i*cb.size(), cb.size());
        evmax_data.push_back(levmax(0));

        Matrix<T, Dynamic, 1> levmin = evmin.segment(cell_i*cb.size(), cb.size());
        evmin_data.push_back(levmin(0));
        cell_i++;
    }
    std::cout << " Postpro time: " << tc.toc() << std::endl;
    std::cout << " L2-norm error (mass): " << std::sqrt(err) << std::endl;;
    std::cout << " L2-norm error (int) : " << std::sqrt(errint) << std::endl;;
    silo.add_variable("dstmesh", "u_dg", u, disk::zonal_variable_t);
    silo.add_variable("dstmesh", "evmax_dg", evmax_data, disk::zonal_variable_t);
    silo.add_variable("dstmesh", "evmin_dg", evmin_data, disk::zonal_variable_t);

    /*
    disk::solvers::feast_eigensolver_params<T> fep;
    fep.verbose = true;
    fep.tolerance = 10;
    fep.min_eigval = 0;
    fep.max_eigval = 200;
    fep.max_iter = 50;
    fep.subspace_size = 20;
    fep.fis = disk::solvers::feast_inner_solver::mumps;

    dynamic_matrix<T> eigvecs;
    dynamic_vector<T> eigvals;
    auto status = disk::solvers::feast(fep, assm.LHS, mass, eigvecs, eigvals);
    std::cout << status << std::endl;
    std::cout << eigvals.transpose() << std::endl;

    int i = 0;
    dynamic_vector<T> eigvals_ref_tmp = dynamic_vector<T>::Zero(81);
    for (size_t m = 1; m < 10; m++) {
        for (size_t n = 1; n < 10; n++) {
            eigvals_ref_tmp(i++) = M_PI*M_PI*(m*m + n*n);
        }
    }

    std::sort(eigvals_ref_tmp.begin(), eigvals_ref_tmp.end());

    dynamic_vector<T> eigvals_ref = eigvals_ref_tmp.segment(0, eigvals.size());
    dynamic_vector<T> eigvals_diff = (eigvals - eigvals_ref).cwiseAbs();
    dynamic_vector<T> eigvals_relerr =
        (100 * eigvals_diff.cwiseAbs()).array() / eigvals_ref.cwiseAbs().array();

    std::cout << eigvals_ref.transpose() << std::endl;
    std::cout << eigvals_relerr.transpose() << std::endl;
    */
    
    return sol;
}

}
