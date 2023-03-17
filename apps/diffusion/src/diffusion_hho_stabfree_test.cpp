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

#include <matplot/matplot.h>

#include "diskpp/mesh/meshgen.hpp"
#include "diskpp/loaders/loader.hpp"
#include "diskpp/methods/hho"
#include "mumps.hpp"
#include "diffusion_hho_common.hpp"

#if 0

#define CSV_SEPARATOR ';'

template<typename T>
struct error_table {
    std::vector<hho_degree_info>    hdis;
    std::vector<T>                  hs;
    std::vector<std::vector<T>>     errors;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const error_table<T>& et)
{
    assert(et.hs.size() == et.errors.size());

    os << CSV_SEPARATOR;
    for (auto& hdi : et.hdis)
    {
        os << "HHO(" << hdi.cell_degree() << ',' << hdi.face_degree() << ',';
        os << hdi.reconstruction_degree() << ")" << CSV_SEPARATOR << CSV_SEPARATOR;
    }
    os << std::endl;

    for (size_t row = 0; row < et.hs.size(); row++)
    {
        os << et.hs.at(row) << CSV_SEPARATOR;

        if (row == 0) {
            const auto& errs = et.errors.at(row); 
            for (size_t col = 0; col < errs.size(); col++)
                os << errs.at(col) << CSV_SEPARATOR << CSV_SEPARATOR;
            os << std::endl;
        }
        else {
            const auto& errs_prev = et.errors.at(row-1);
            const auto& errs_curr = et.errors.at(row);
            assert(errs_prev.size() == errs_curr.size());

            const auto& h_prev = et.hs.at(row-1);
            const auto& h_curr = et.hs.at(row);
            for (size_t col = 0; col < errs_curr.size(); col++) {
                os << errs_curr.at(col) << CSV_SEPARATOR;
                auto rate = std::log(errs_prev.at(col)/errs_curr.at(col)) / std::log(h_prev/h_curr);
                os << rate << CSV_SEPARATOR;
            }
            os << std::endl;
        }
    }

    return os;
}

#endif

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

    std::map<hdi_t, varhist_t> convergence_histories_byhdi;

public:
    convergence_database() {}

    void add(const hdi_t& hdi, const variant_t& variant, const errinfo_t& ch) {
        convergence_histories_byhdi[hdi][variant].add(ch);
    }

    void mark_failed(const hdi_t& hdi, const variant_t& variant) {
        convergence_histories_byhdi[hdi][variant].failed(true);
    }

    auto begin() const {
        return convergence_histories_byhdi.begin();
    }

    auto end() const {
        return convergence_histories_byhdi.end();
    }
};

template<typename T>
void make_images(const convergence_database<T>& cdb)
{
    using namespace matplot;

    for (const auto& [hdi, variant_histories] : cdb) {
        figure();

        /* L2 errors */
        subplot(3,1,0);
        hold(on);
        for (const auto& [variant, history] : variant_histories) {
            if ( history.failed() )
                continue;

            auto p = loglog(history.hs, history.L2errs);
            p->display_name(variant);
        }
        auto l = legend();
        l->location(legend::general_alignment::bottomleft);
        hold(off);

        /* H1 errors */
        subplot(3,1,1);
        hold(on);
        for (const auto& [variant, history] : variant_histories) {
            if ( history.failed() )
                continue;

            auto p = loglog(history.hs, history.H1errs);
            p->display_name(variant);
        }
        l = legend();
        l->location(legend::general_alignment::bottomleft);
        hold(off);

        /* Operator errors */
        subplot(3,1,2);
        hold(on);
        for (const auto& [variant, history] : variant_histories) {
            if ( history.failed() )
                continue;

            auto p = loglog(history.hs, history.Aerrs);
            p->display_name(variant);
        }
        l = legend();
        l->location(legend::general_alignment::bottomright);
        hold(off);


        std::stringstream title_ss;
        title_ss << "HHO(" << hdi.cell_degree() << ", " << hdi.face_degree() << ", ";
        title_ss << hdi.reconstruction_degree() << ")";
        sgtitle(title_ss.str());

        std::stringstream ss;
        ss << "img/convergence_" << hdi.cell_degree() << "_" << hdi.face_degree() << "_";
        ss << hdi.reconstruction_degree() << ".png";
        save(ss.str());
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

    test_configuration(const std::string& vn) : k_min(0), k_max(3), rincrmin(1),
        rincrmax(4), mixed_order(false), use_projection(false), k00(1.0), k11(1.0),
        variant_name(vn)
    {}
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
                cdb.add(hdi, tc.variant_name, error);
            }
            catch (...) {
                cdb.mark_failed(hdi, tc.variant_name);
            }
        }
    }
}

int main(int argc, char **argv)
{
    rusage_monitor rm;

    std::vector<char *> meshes;

    bool meshes_from_file = false;
    bool use_proj = false;
    int ch;

    double k00 = 0.008;
    double k11 = 1.0;

    while ( (ch = getopt(argc, argv, "m:Px:y:")) != -1 )
    {
        switch(ch)
        {
            case 'm':
                meshes_from_file = true;
                meshes.push_back(optarg);
                break;

            case 'P':
                use_proj = true;
                break;

            case 'x':
                k00 = std::stod(optarg);
                break;
            
            case 'y':
                k11 = std::stod(optarg);
                break;
        }
    }

    test_configuration plain_hho("HHO-E");

    test_configuration mixed_hho("HHO-M");
    mixed_hho.mixed_order = true;

    test_configuration mixed_proj_hho("HHO-MP");
    mixed_proj_hho.mixed_order = true;
    mixed_proj_hho.use_projection = true;

    convergence_database<double> cdb;

    if (meshes_from_file)
    {
        for (const auto& mesh : meshes)
        {
            disk::dispatch_all_meshes(mesh,
                [](auto& ...args) { test_stabfree_hho(args...); },
                cdb, plain_hho );

            disk::dispatch_all_meshes(mesh,
                [](auto& ...args) { test_stabfree_hho(args...); },
                cdb, mixed_hho );

            disk::dispatch_all_meshes(mesh,
                [](auto& ...args) { test_stabfree_hho(args...); },
                cdb, mixed_proj_hho );
        }
    }
    else
    {
        using T = double;
        disk::triangular_mesh<T> msh;
        auto mesher = make_simple_mesher(msh);

        for (size_t i = 0; i < 7; i++)
        {
            mesher.refine();
            test_stabfree_hho(msh, cdb, plain_hho);
            test_stabfree_hho(msh, cdb, mixed_hho);
            test_stabfree_hho(msh, cdb, mixed_proj_hho);
        }
    }

    make_images(cdb);

    return 0;
}

