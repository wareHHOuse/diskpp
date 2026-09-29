/*
 *
 *      (1/kappa) p_tt  -  div(grad p) = 0     in (0,1)^2 x (0,T],
 *                                  p = 0        on the boundary,
 *
*/

#include <iostream>
#include <regex>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <numbers>

#include "diskpp/loaders/loader.hpp"
#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/solvers/direct_solvers.hpp" 
#include "diskpp/common/timecounter.hpp"

#include "diskpp/bases/bases_new.hpp"
#include "diskpp/bases/bases_operations.hpp"
#include "diskpp/methods/hho_slapl.hpp"
#include "diskpp/methods/hho_assemblers.hpp"

using namespace disk;

struct params {
    size_t degree       = 1;      // HHO polynomial degree k
    double kappa        = 1.0;    // bulk modulus (wave speed c = sqrt(kappa))
    double final_time   = 1.0;
    size_t nsteps       = 20;     // time steps at the COARSEST level
    size_t levels       = 4;      // refinement levels per mesh family
    double beta         = 0.4;   // Newmark (1/4, 1/2): average acceleration,
    double gamma        = 0.7;    // 2nd order, unconditionally stable, no dissip.

    bool fixed_dt = false;
    size_t steps_growth = 5;
    
};

struct exact_solution 
{
    static constexpr double omega = std::numbers::sqrt2 * M_PI;
    static constexpr double amp   = 1.0 / omega;

    static double p(double x, double y, double t)
    {
        return amp * std::sin(omega * t) * std::sin(M_PI * x) * std::sin(M_PI * y);
    }

    static double dpdt(double x, double y, double t)
    {
        return amp * omega * std::cos(omega * t) * std::sin(M_PI * x) * std::sin(M_PI * y); 
    }
};


template<typename Mesh>
struct spatial_operator {
    using scalar = typename Mesh::coordinate_type;
    using matrix = disk::dynamic_matrix<scalar>;

    const Mesh&                   msh;      // for condensed_solver + assembler
    disk::hho::slapl::degree_info di;       // for make_assembler
    scalar                        inv_kappa;
    std::vector<matrix>           K;        // b_h cellwise (consistency + stab)
    std::vector<matrix>           M;        // cell mass *1/kappa
    std::vector<size_t>           ncell;    // nT cellwise

    spatial_operator(const Mesh& msh, const params& prm)
        : msh(msh), di(prm.degree), inv_kappa(scalar(1)/prm.kappa)
    {
        using namespace disk::hho::slapl;
        const size_t N = msh.cells_size();
        K.reserve(N);  M.reserve(N);  ncell.reserve(N);

        for (auto& cl : msh)
        {
            auto [R, A] = local_operator(msh, cl, di);          // rec + cons
            auto stab   = local_stabilization(msh, cl, di, R);  // stab TODO: possibility of change between HDG and HHO
            K.push_back(A + stab);                              // b_h local

            auto phi = typename hho_space<Mesh>::cell_basis_type(msh, cl, di.cell);
            M.push_back(inv_kappa * integrate(msh, cl, phi, phi)); //mass for timestep
            ncell.push_back(phi.size());
        }
    }
};


/*
* Condensed solver
*/

template<typename Mesh>
struct condensed_solver {

    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;
    using matrix = disk::dynamic_matrix<scalar>;

    #ifdef HAVE_PARDISO
        using global_solver = Eigen::PardisoLDLT<Eigen::SparseMatrix<scalar>>;
    #else
        using global_solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<scalar>>;
    #endif
        using assembler_type =
            decltype(disk::hho::slapl::make_assembler(std::declval<const Mesh&>(),
                                                  std::declval<const disk::hho::slapl::degree_info&>()));

    const spatial_operator<Mesh>& op; 

    mutable assembler_type            assembler;   
    std::vector<Eigen::LDLT<matrix>>  STT_solver;   // S_TT^{-1} aplicable, cellwise
    std::vector<matrix>               S_FT;         // cellwise S_FT
    std::vector<matrix>               STTinv_STF;   // S_TT^{-1} S_TF, cellwise (deschur)
    global_solver                     fact; 

    condensed_solver(const spatial_operator<Mesh>& sp_op) 
    : op(sp_op)
    , assembler(disk::hho::slapl::make_assembler(sp_op.msh, sp_op.di))
    { }
    
    /*
    * ==================================================
    * FACTORIZATION 
    * ==================================================
    */
    void factor(const std::vector<matrix>& S)
    {
        const size_t N = op.K.size();

        STT_solver.resize(N);
        S_FT.resize(N);
        STTinv_STF.resize(N);

        size_t c = 0;
        for (auto& cl : op.msh)
        {
            const size_t nT = op.ncell[c];
            const size_t nF = S[c].rows() - nT; // faces

            // ======= local variables ========

            matrix S_TT = S[c].block(0,0,nT,nT);  // upper-left
            matrix S_TF = S[c].block(0,nT,nT,nF); // upper-right
            matrix S_FF = S[c].block(nT,nT,nF,nF); // bottom-right

            // ======= export variables =======

            S_FT[c] = S[c].block(nT,0,nF,nT); //bottom-left
            STT_solver[c] = S_TT.ldlt();  // ready to computa S_TT^{-1}
            STTinv_STF[c] = STT_solver[c].solve(S_TF); // S_TT^{-1}*S_TF

            // ======== schur ========
            matrix local_schur = S_FF - S_FT[c]*STTinv_STF[c];
            vector zero_rhs = vector::Zero(nF); //we need only the LHS so we put zero RHS
            assembler.assemble(op.msh, cl, local_schur, zero_rhs);
            ++c;
        }
        // ====== global factorization ======
        assembler.finalize();
        fact.compute(assembler.LHS);
        if (fact.info() != Eigen::Success)
        throw std::runtime_error("factor failed");

    }   

    /*
    * ==================================================
    * SOLVE TIMESTEP 
    * ==================================================
    */
    std::vector<vector>  solve (const std::vector<vector>& R) const
    {
        const size_t N = op.K.size();

        std::vector<vector> STTinv_RT(N);

        assembler.RHS.setZero();

        size_t c = 0;
        for (auto& cl : op.msh)
        {
            const size_t nT = op.ncell[c];
            const size_t nF = S_FT[c].rows();

            STTinv_RT[c] = STT_solver[c].solve(R[c].head(nT)); // S_TT^{-1} R_T
            vector condensed_rhs = R[c].tail(nF) - S_FT[c]*STTinv_RT[c]; // R_F - S_FT*S_TT^{-1}*R_T
            matrix zero_lhs = matrix::Zero(nF,nF);
            assembler.assemble(op.msh, cl, zero_lhs, condensed_rhs);
            ++c;
        }
        assembler.finalize();

        vector a_F_global = fact.solve(assembler.RHS);

        std::vector<vector> a_next(N);
        c = 0;
        for (auto& cl : op.msh)
        {
            const size_t nT = op.ncell[c];
            vector a_F = assembler.take_local_solution(op.msh,cl,a_F_global);
            const size_t nF = a_F.size(); // faces

            a_next[c].resize(nT + nF);
            a_next[c].head(nT) = STTinv_RT[c] - STTinv_STF[c]*a_F;
            a_next[c].tail(nF) = a_F;
            ++c;
        }

        return a_next;

    }
}; // void condensed_solver

/*
* One time step
*/

template<typename Mesh>
struct newmark_step { 
    
    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;
    using matrix = disk::dynamic_matrix<scalar>;

    const spatial_operator<Mesh>& op; // for M and K = b_h

    scalar dt, beta, gamma;

    condensed_solver<Mesh> solver;

    newmark_step(const spatial_operator<Mesh>& sp_op, const params& prm)
        : op(sp_op), dt(prm.final_time/prm.nsteps), beta(prm.beta), gamma(prm.gamma), solver(sp_op)
    {   
        const size_t N = op.K.size();
        std::vector<matrix> S(N);
        
        for (size_t c = 0; c < N; ++c)
        {
            const size_t nT = op.ncell[c];    // cell
            S[c] = beta * dt * dt * op.K[c];  // beta*dt^2*K
            S[c].block(0,0,nT,nT) += op.M[c]; // adding M only in the cell block
        }
        solver.factor(S);
    }

    // mute p, v, a
    void step(std::vector<vector>& p, std::vector<vector>& v, std::vector<vector>& a) const 
    {
        const size_t N = op.K.size();
        std::vector<vector> rhs(N), pstar(N), vstar(N);
        
        // predictor
        for (size_t c = 0; c < N; ++c)
        {
            const size_t nT = op.ncell[c];
            pstar[c] = p[c] + dt * v[c] + 0.5 * dt * dt * (1.0 - 2.0*beta) * a[c];
            vstar[c] = v[c] + dt * (1.0 - gamma) * a[c];
            rhs[c] = -op.K[c] * pstar[c];
        }

        // solve
        auto a_next = solver.solve(rhs);

        // corrector
        for (size_t c = 0; c < N; ++c)
        {
            p[c] = pstar[c] + beta * dt * dt * a_next[c];
            v[c] = vstar[c] + gamma * dt * a_next[c];
            a[c] = a_next[c];

        }
    }

};

template<typename Mesh>
std::pair<double, double>
errors(const spatial_operator<Mesh>& op,
       const std::vector<disk::dynamic_vector<typename Mesh::coordinate_type>>& p,
       double t)
{
    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;
    using matrix = disk::dynamic_matrix<scalar>;
    using cell_basis = typename disk::hho::slapl::hho_space<Mesh>::cell_basis_type;

    scalar l2_num = 0, l2_den = 0, h1_num = 0, h1_den = 0;
    std::function<scalar(const typename Mesh::point_type&)> pex =
        [&](const auto& pt) { return exact_solution::p(pt.x(), pt.y(), t); };

    size_t c = 0;
    for (auto& cl : op.msh) {                       
        const size_t nT = op.ncell[c];
        vector Ip = disk::hho::slapl::local_reduction(op.msh, cl, op.di, pex);
        vector e  = p[c] - Ip;

        cell_basis phi(op.msh, cl, op.di.cell);
        matrix mass = integrate(op.msh, cl, phi, phi);
        vector eT = e.head(nT), IpT = Ip.head(nT);
        l2_num += eT.dot(mass * eT);
        l2_den += IpT.dot(mass * IpT);

        h1_num += e.dot(op.K[c] * e);
        h1_den += Ip.dot(op.K[c] * Ip);
        ++c;
    }
    return { std::sqrt(l2_num / std::max(l2_den, scalar(1e-30))),
             std::sqrt(h1_num / std::max(h1_den, scalar(1e-30))) };
}

template<typename Mesh>
std::pair<double,double> run_acoustic(const Mesh& msh, const params& prm)
{
    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;

    spatial_operator<Mesh> op(msh, prm);
    newmark_step<Mesh> stepper(op, prm);

    const size_t N = msh.cells_size();
    std::vector<vector> p(N), v(N), a(N);

    std::function<scalar(const typename Mesh::point_type&)> p0 =
        [](const auto& pt) { return exact_solution::p(pt.x(), pt.y(), 0.0); };
    std::function<scalar(const typename Mesh::point_type&)> v0 =
        [](const auto& pt) { return exact_solution::dpdt(pt.x(), pt.y(), 0.0); };

    size_t c = 0;
    for (auto& cl : msh) {
        p[c] = disk::hho::slapl::local_reduction(msh, cl, op.di, p0);
        v[c] = disk::hho::slapl::local_reduction(msh, cl, op.di, v0);
        ++c;
    }

    for (size_t c = 0; c < N; ++c)
    {
        a[c] = vector::Zero(p[c].size());
    }

    for (size_t n = 0; n < prm.nsteps; ++n)
    {
        stepper.step(p, v, a);
    }

    return errors(op, p, prm.final_time);
}

template<class MakeMesh>
void energy_trace(const std::string& name, MakeMesh make_mesh,
                  params prm, size_t L)
{
    using Mesh   = decltype(make_mesh(L));
    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;

    auto msh = make_mesh(L);                       // fixed mesh

    spatial_operator<Mesh> op(msh, prm);
    newmark_step<Mesh>     stepper(op, prm);

    const size_t N = msh.cells_size();
    std::vector<vector> p(N), v(N), a(N);

    std::function<scalar(const typename Mesh::point_type&)> p0 =
        [](const auto& pt) { return exact_solution::p(pt.x(), pt.y(), 0.0); };
    std::function<scalar(const typename Mesh::point_type&)> v0 =
        [](const auto& pt) { return exact_solution::dpdt(pt.x(), pt.y(), 0.0); };

    size_t c = 0;
    for (auto& cl : msh) {
        p[c] = disk::hho::slapl::local_reduction(msh, cl, op.di, p0);
        v[c] = disk::hho::slapl::local_reduction(msh, cl, op.di, v0);
        a[c] = vector::Zero(p[c].size());
        ++c;
    }

    const double dt = prm.final_time / prm.nsteps;
    const double E0 = energy(op, p, v);

    std::cout << "\n=== energy: " << name << "  (k=" << prm.degree
              << ", L=" << L << ", nsteps=" << prm.nsteps
              << ", beta=" << prm.beta << ", gamma=" << prm.gamma << ") ===\n"
              << std::setw(6) << "step" << std::setw(14) << "t"
              << std::setw(18) << "E" << std::setw(16) << "(E-E0)/E0" << "\n";

    std::cout << std::setw(6) << 0 << std::setw(14) << std::fixed << std::setprecision(4) << 0.0
              << std::setw(18) << std::scientific << std::setprecision(8) << E0
              << std::setw(16) << 0.0 << "\n";

    for (size_t n = 0; n < prm.nsteps; ++n) {
        stepper.step(p, v, a);
        const double E = energy(op, p, v);
        std::cout << std::setw(6) << n+1
                  << std::setw(14) << std::fixed << std::setprecision(4) << (n+1)*dt
                  << std::setw(18) << std::scientific << std::setprecision(8) << E
                  << std::setw(16) << (E - E0)/E0 << "\n";
    }
}

template<class MakeMesh>
void convergence_table(const std::string& name, MakeMesh make_mesh,
                       params prm, size_t levels)
{
    std::cout << "\n=== " << name << "  (k = " << prm.degree << ") ===\n"
          << std::setw(6)  << "level" << std::setw(14) << "h"
          << std::setw(16) << "L2 error" << std::setw(8) << "rate"
          << std::setw(16) << "H1 error" << std::setw(8) << "rate" << "\n";

    const size_t base_nsteps = prm.nsteps;
    double h_prev = 0, l2_prev = 0, h1_prev = 0;

    for (size_t L = 0; L < levels; ++L) {
        auto   msh = make_mesh(L);
        double h   = disk::average_diameter(msh);

        prm.nsteps = base_nsteps;

        if (!prm.fixed_dt)
            for (size_t i = 0; i < L; ++i) prm.nsteps *= prm.steps_growth;

        auto [eL2, eH1] = run_acoustic(msh, prm);

        std::cout << std::setw(6) << L
                  << std::setw(14) << std::scientific << std::setprecision(4) << h
                  << std::setw(16) << eL2;
        if (L > 0) std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                             << std::log(l2_prev/eL2) / std::log(h_prev/h);
        else       std::cout << std::setw(8) << "--";
        std::cout << std::setw(16) << std::scientific << std::setprecision(4) << eH1;
        if (L > 0) std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                             << std::log(h1_prev/eH1) / std::log(h_prev/h);
        else       std::cout << std::setw(8) << "--";
        std::cout << "\n";

        h_prev = h; l2_prev = eL2; h1_prev = eH1;
    }
}

template<typename Mesh>
double energy(const spatial_operator<Mesh>& op,
              const std::vector<disk::dynamic_vector<typename Mesh::coordinate_type>>& p,
              const std::vector<disk::dynamic_vector<typename Mesh::coordinate_type>>& v)
{
    using scalar = typename Mesh::coordinate_type;
    using vector = disk::dynamic_vector<scalar>;

    scalar kin = 0, pot = 0;
    const size_t N = op.K.size();
    for (size_t c = 0; c < N; ++c) {
        const size_t nT = op.ncell[c];
        vector vT = v[c].head(nT);                 // dtp_T
        kin += vT.dot(op.M[c] * vT);               // || dtp_T ||_L2(1/kappa;Omega)^2
        pot += p[c].dot(op.K[c] * p[c]);           // b(ph, ph)
    }
    return 0.5 * (kin + pot);
}

int main(int argc, char** argv)
{
    using T = double;

    params prm;

    // triangular (simplicial) meshes of the unit square
    auto make_tri = [](size_t L) {
        disk::simplicial_mesh<T, 2> msh;
        auto mesher = disk::make_simple_mesher(msh);
        for (size_t i = 0; i <= L; ++i) mesher.refine();
        return msh;
    };
    // quadrilateral (cartesian) meshes
    auto make_quad = [](size_t L) {
        disk::cartesian_mesh<T, 2> msh;
        auto mesher = disk::make_simple_mesher(msh);
        for (size_t i = 0; i <= L; ++i) mesher.refine();
        return msh;
    };
    // polygonal (FVCA5 hexagonal) meshes
    auto make_hex = [](size_t L) {
        disk::generic_mesh<T, 2> msh;
        auto mesher = disk::make_fvca5_hex_mesher(msh);
        mesher.make_level(L);
        return msh;
    };

    convergence_table("triangular",        make_tri,  prm, prm.levels);
    convergence_table("quadrilateral",     make_quad, prm, prm.levels);
    convergence_table("polygonal (FVCA5)", make_hex,  prm, prm.levels);

    params eprm = prm;
    eprm.nsteps = 200;            // fijo; resolución temporal para la traza
    const size_t Lfix = 2;        // malla fija

    energy_trace("triangular",        make_tri,  eprm, Lfix);
    energy_trace("quadrilateral",     make_quad, eprm, Lfix);
    energy_trace("polygonal (FVCA5)", make_hex,  eprm, Lfix);

    return 0;
}
