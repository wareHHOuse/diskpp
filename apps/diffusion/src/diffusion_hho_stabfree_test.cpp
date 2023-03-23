/*
 * DISK++, a template library for DIscontinuous SKeletal methods.
 *
 * Matteo Cicuttin (C) 2023
 * matteo.cicuttin@polito.it
 *
 * Politecnico di Torino - DISMA
 * Dipartimento di Matematica
 */

#include <iostream>
#include <vector>
#include <regex>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <map>

//#include <matplot/matplot.h>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "mumps.hpp"
#include "diffusion_hho_common.hpp"

#include "sgr.hpp"


template<typename T>
struct convergence_history {
    std::vector<T>  hs;
    std::vector<T>  L2errs;
    std::vector<T>  H1errs;
    std::vector<T>  Aerrs;
    bool            has_failed;

    constexpr static double L2MAX = 1e8;
    constexpr static double H1MAX = 1e8;
    constexpr static double AMAX = 1e8;

    using error_info_type = error_info<T>;

    convergence_history() : has_failed(false)
    {}

    bool failed(void) const {
        return has_failed;
    }

    void failed(bool hf) {
        has_failed = hf;
    }

    bool add(const error_info<T>& ei) {
        if (ei.L2err > L2MAX or ei.H1err > H1MAX or ei.Aerr > AMAX) {
            failed(true);
            return false;
        }
        hs.push_back( ei.h );
        L2errs.push_back( ei.L2err );
        H1errs.push_back( ei.H1err );
        Aerrs.push_back( ei.Aerr );
        return true;
    }
};

template<typename T>
class convergence_database {
    using variant_t = std::string;
    using hdi_t = disk::hho_degree_info;
    using convhist_t = convergence_history<T>;
    using errinfo_t = typename convhist_t::error_info_type;
    
    using varhist_t = std::map<variant_t, convhist_t>;

    std::map<size_t, varhist_t> convergence_histories;

public:
    convergence_database() {}

    void add(size_t deg, const variant_t& variant, const errinfo_t& ch) {
        convergence_histories[deg][variant].add(ch);
    }

    void mark_failed(size_t deg, const variant_t& variant) {
        convergence_histories[deg][variant].failed(true);
    }

    auto begin() const {
        return convergence_histories.begin();
    }

    auto end() const {
        return convergence_histories.end();
    }
};

#if 0
template<typename T>
void make_images(const convergence_database<T>& cdb)
{
    using namespace matplot;

    auto style = [] (const hho_degree_info& hdi) {
        auto cd = hdi.cell_degree();
        auto fd = hdi.face_degree();
        auto rd = hdi.reconstruction_degree();

        if (cd == fd and fd+1 == rd)
            return '.'; /* equal order */
        
        if (cd+1 == fd and fd+1 == rd)
            return '-'; /* mixed order */
        
        return 'x';
    };

    for (const auto& [degree, variants] : cdb) {
        figure();
        /* L2 errors */
        hold(on);
        for (const auto& [variant, history] : variants) {
            if ( history.failed() )
                continue;
            std::cout << "Plot: " << variant << std::endl;
            auto p = loglog(history.hs, history.L2errs);
            p->display_name( variant );
        }
        auto l = legend();
        l->location(legend::general_alignment::bottomright);
        hold(off);

        {
            std::stringstream title_ss;
            title_ss << "Degree " << degree << " L2 error";
            sgtitle(title_ss.str());

            std::stringstream ss;
            ss << "img/convergence_" << degree << "_L2.eps";
            save(ss.str());
        }

        /* H1 errors */
        figure();
        hold(on);
        for (const auto& [variant, history] : variants) {
            if ( history.failed() )
                continue;

            auto p = loglog(history.hs, history.H1errs);
            p->display_name( variant );
        }
        l = legend();
        l->location(legend::general_alignment::bottomright);
        hold(off);

        {
            std::stringstream title_ss;
            title_ss << "Degree " << degree << " energy error";
            sgtitle(title_ss.str());

            std::stringstream ss;
            ss << "img/convergence_" << degree << "_energy.eps";
            save(ss.str());
        }

        /* Operator errors */
        figure();
        hold(on);
        for (const auto& [variant, history] : variants) {
            if ( history.failed() )
                continue;

            auto p = loglog(history.hs, history.Aerrs);
            p->display_name( variant );
        }
        l = legend();
        l->location(legend::general_alignment::bottomright);
        hold(off);

        {
            std::stringstream title_ss;
            title_ss << "Degree " << degree << " operator error";
            sgtitle(title_ss.str());

            std::stringstream ss;
            ss << "img/convergence_" << degree << "_operator.eps";
            save(ss.str());
        }
    }
}
#endif

template<typename T>
auto compute_orders(const std::vector<T>& hs, const std::vector<T>& errs)
{
    assert(hs.size() == errs.size());
    std::vector<T> ret;
    for (size_t i = 1; i < hs.size(); i++) {
        auto num = std::log( errs.at(i-1)/errs.at(i) );
        auto den = std::log( hs.at(i-1)/hs.at(i) );
        ret.push_back(num/den);
    }

    return ret;
}

namespace mpriv {
template<typename T>
struct range_check
{
    T value;
    T lower_limit, upper_limit;
    
    range_check(const T& val, const T& tol)
        : value(val), lower_limit(-tol),upper_limit(tol)
    {}

    range_check(const T& val, const std::array<T,2>& minmax)
        : value(val), lower_limit(minmax[0]), upper_limit(minmax[1])
    {}
};
}//namespace mpriv

template<typename T>
std::ostream&
operator<<(std::ostream& os, const mpriv::range_check<T>& rc) 
{
    if (rc.value < rc.lower_limit)
        os << sgr::redfg << rc.value << sgr::nofg;
    else if (rc.value > rc.upper_limit)
        os << sgr::yellowfg << rc.value << sgr::nofg;
    else
        os << sgr::greenfg << rc.value << sgr::nofg;
    return os;
}

template<typename T>
inline mpriv::range_check<T>
range_check(const T& val, const T& tol)
{
    return mpriv::range_check<T>(val, tol);
}

template<typename T>
inline mpriv::range_check<T>
range_check(const T& val, const std::array<T,2>& minmax)
{
    return mpriv::range_check<T>(val, minmax);
}

template<typename T>
void make_reports(const convergence_database<T>& cdb)
{
    std::ofstream ofs("report.txt");

    for (const auto& [degree, variants] : cdb)
    {
        for (const auto& [variant, history] : variants)
        {
            std::cout << sgr::Bon << " * * * " << variant << " * * *" << sgr::Boff << std::endl;
            if ( history.failed() ) {
                std::cout << "This variant failed to run" << std::endl;
                continue;
            }
            std::cout << " Expected order for L2: " << degree+2 << std::endl << '\t';
            auto L2orders = compute_orders(history.hs, history.L2errs);
            for (const auto& order : L2orders)
                std::cout << range_check(order, {degree+2-0.3, degree+2+0.3}) << " ";
            std::cout << std::endl;
            std::cout << " Expected order for H1: " << degree+1 << std::endl << '\t';
            auto H1orders = compute_orders(history.hs, history.H1errs);
            for (const auto& order : H1orders)
                std::cout << range_check(order, {degree+1-0.3, degree+1+0.3}) << " ";
            std::cout << std::endl;
        }
    }
}

struct test_configuration
{
    size_t          k_min;
    size_t          k_max;
    size_t          rincrmin;
    size_t          rincrmax;
    bool            mixed_order;
    bool            use_projection;
    double          k00;
    double          k11;
    std::string     variant_name;

    test_configuration() : k_min(0), k_max(3), rincrmin(1), rincrmax(4),
        mixed_order(false), use_projection(false), k00(1.0), k11(1.0)
    {}

    test_configuration(const test_configuration&) = default;

    test_configuration(const std::string& vn) : k_min(0), k_max(3), rincrmin(1),
        rincrmax(4), mixed_order(false), use_projection(false), k00(1.0), k11(1.0),
        variant_name(vn)
    {}

    std::string name(void) const { 
        return variant_name;
    }

    void name(const std::string& name) {
        variant_name = name;
    }
};

template<typename Mesh>
void
test_stabfree_hho(Mesh& msh, convergence_database<typename Mesh::coordinate_type>& cdb,
    const test_configuration& tc)
{
    using T = typename Mesh::coordinate_type;

    disk::diffusion_tensor<Mesh> diff_tens = disk::diffusion_tensor<Mesh>::Zero();
    diff_tens(0,0) = tc.k00;
    diff_tens(1,1) = tc.k11;

    auto make_display_name = [] (const std::string& variant, const hho_degree_info& hdi) {
        auto cd = hdi.cell_degree();
        auto fd = hdi.face_degree();
        auto rd = hdi.reconstruction_degree();
        std::stringstream ss;
        ss << variant << '(' << cd << ',' << fd << ',' << rd << ')';
        return ss.str();
    };

    for (size_t k = tc.k_min; k <= tc.k_max; k++)
    {
        for (size_t r = k+tc.rincrmin; r <= k+tc.rincrmax; r++)
        {
            hho_degree_info hdi;
            if (tc.mixed_order)
                hdi.cell_degree(k+1);
            else
                hdi.cell_degree(k);

            hdi.face_degree(k);
            hdi.reconstruction_degree(r);

            try {
                auto error = run_hho_diffusion_solver(msh, hdi, true, true, diff_tens, tc.use_projection);
                cdb.add(hdi.face_degree(), make_display_name(tc.variant_name, hdi), error);
            }
            catch (...) {
                cdb.mark_failed(hdi.face_degree(), make_display_name(tc.variant_name, hdi));
            }
        }
    }
}

int main(int argc, char **argv)
{
    using T = double;
    rusage_monitor rm;

    std::vector<char *> meshes;
    bool meshes_from_file = false;
    size_t max_refinements = 6;
    
    std::string shape = "tri";

    test_configuration default_test_config;
    default_test_config.k00 = 1.0;
    default_test_config.k11 = 1.0;

    int ch;
    while ( (ch = getopt(argc, argv, "m:x:y:h:k:K:r:R:s:")) != -1 )
    {
        switch(ch)
        {
            case 'm':
                meshes_from_file = true;
                meshes.push_back(optarg);
                break;

            case 'x':
                default_test_config.k00 = std::stod(optarg);
                break;
            
            case 'y':
                default_test_config.k11 = std::stod(optarg);
                break;

            case 'h':
                max_refinements = std::stoi(optarg);
                break;

            case 'k':
                default_test_config.k_min = std::stoi(optarg);
                break;

            case 'K':
                default_test_config.k_max = std::stoi(optarg);
                break;

            case 'r':
                default_test_config.rincrmin = std::stoi(optarg);
                break;

            case 'R':
                default_test_config.rincrmax = std::stoi(optarg);
                break;
                
            case 's':
                shape = optarg;
                break;
        }
    }

    test_configuration plain_hho(default_test_config);
    plain_hho.variant_name = "HHO-E";

    test_configuration mixed_hho(default_test_config);
    mixed_hho.variant_name = "HHO-M";
    mixed_hho.mixed_order = true;

    test_configuration mixed_proj_hho(default_test_config);
    mixed_proj_hho.variant_name = "HHO-MP";
    mixed_proj_hho.mixed_order = true;
    mixed_proj_hho.use_projection = true;

    convergence_database<T> cdb;

    if (shape == "tri") {
        disk::triangular_mesh<T> msh;
        auto mesher = make_simple_mesher(msh);
        mesher.refine();

        for (size_t i = 0; i < max_refinements; i++)
        {
            test_stabfree_hho(msh, cdb, plain_hho);
            test_stabfree_hho(msh, cdb, mixed_hho);
            test_stabfree_hho(msh, cdb, mixed_proj_hho);
        }
    }

    if (shape == "hex") {
        disk::generic_mesh<T,2> msh;
        auto mesher = make_fvca5_hex_mesher(msh);

        for (size_t i = 0; i < max_refinements; i++)
        {
            const auto offset = 1;
            mesher.make_level(i+offset);
            test_stabfree_hho(msh, cdb, plain_hho);
            test_stabfree_hho(msh, cdb, mixed_hho);
            test_stabfree_hho(msh, cdb, mixed_proj_hho);
        }
    }

    //make_images(cdb);
    make_reports(cdb);

    return 0;
}

