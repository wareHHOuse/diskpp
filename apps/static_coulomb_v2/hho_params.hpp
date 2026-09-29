#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace hho_contact {

/* Which Nitsche version to use on the contact boundary.
 *
 *   face : the displacement in P_n(u) = sigma_n(u) - gamma u_n is the FACE
 *          unknown u_F.  Contact faces keep their unknowns. Is the default.
 *
 *   cell : the displacement is the trace of the CELL unknown u_T|_F, and the
 *          face unknowns on the contact boundary are removed.
 *
 * Both run on the same mixed-order discretization (cell k+1, face k), so the
 * only difference between them is the trace source plus that DOF removal.  */
enum class trace_variant { face, cell };

struct solver_opts
{
    double tol_newton     = 1e-10;   // absolute, on |du_faces|
    double tol_newton_rel = 1e-8;    // relative, on |du_faces| / |u_faces|
    double tol_picard     = 1e-8;    // on |u - u_picard| over all dofs
    size_t max_newton     = 40;
    size_t max_picard     = 60;
    bool   verbose        = true;
};

struct Params
{
    // discretization
    size_t k       = 1;
    trace_variant variant = trace_variant::face;   // 'variant = cell' in the .dat

    /* Nitsche penalty.  gamma_0 is DIMENSIONLESS
     *    the code turns it into gamma_n = gamma_n_0 * 2 mu / h_F,   gamma_t = gamma_t_0 * 2 mu / h_F
     * (see disk::nitsche_gamma in contact_terms.hpp)
     *
     * gamma_0 sets both; gamma_n_0 / gamma_t_0 override the normal and the
     * tangential one individually, in either order in the .dat file.  A
     * value <= 0 means "not set, take gamma_0".                            */
    double gamma_0   = 10.0;
    double gamma_n_0 = 0.0;          // 0 -> gamma_0
    double gamma_t_0 = 0.0;          // 0 -> gamma_0
    double gamma_s_0 = 10.0; 
    double theta     = -1.0;

    double gn0() const { return gamma_n_0 > 0.0 ? gamma_n_0 : gamma_0; }
    double gt0() const { return gamma_t_0 > 0.0 ? gamma_t_0 : gamma_0; }
    double gs0() const { return gamma_s_0 > 0.0 ? gamma_s_0 : gamma_0; }

    // material: (E, nu) or (mu, lam) directly
    double E   = 0.0;
    double nu  = 0.0;
    double mu  = 0.0;
    double lam = 0.0;

    // friction datum: Tresca s or Coulomb F, experiment decides
    double friction = 0.0;

    // mesh
    std::string mesh_path;
    size_t      mesh_level = 3;      // for built-in meshers

    solver_opts opts;

    void derive_lame()
    {
        if (mu == 0.0 && lam == 0.0 && E > 0.0)
        {
            mu  = E / (2.0*(1.0 + nu));
            lam = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu));   // plane strain 
        }
    }

    static Params load(const std::string& path)
    {
        Params p;
        std::ifstream ifs(path);
        if (!ifs)
        {
            std::cout << "Params: cannot open '" << path
                      << "', using defaults\n";
            return p;
        }

        auto trim = [](std::string s) {
            auto notspace = [](unsigned char c){ return !std::isspace(c); };
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
            s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
            return s;
        };

        std::string line;
        while (std::getline(ifs, line))
        {
            const auto hash = line.find('#');
            if (hash != std::string::npos) line = line.substr(0, hash);
            const auto eq = line.find('=');
            if (eq == std::string::npos) continue;

            const std::string key = trim(line.substr(0, eq));
            const std::string val = trim(line.substr(eq + 1));
            if (key.empty() || val.empty()) continue;

            if      (key == "k")              p.k          = std::stoul(val);
            else if (key == "variant")
            {
                if      (val == "cell") p.variant = trace_variant::cell;
                else if (val == "face") p.variant = trace_variant::face;
                else std::cout << "Params: variant must be 'face' or 'cell', "
                                  "got '" << val << "', keeping face\n";
            }
            else if (key == "gamma_0")        p.gamma_0    = std::stod(val);
            else if (key == "gamma_n_0")      p.gamma_n_0  = std::stod(val);
            else if (key == "gamma_t_0")      p.gamma_t_0  = std::stod(val);
            else if (key == "gamma_s_0")      p.gamma_s_0  = std::stod(val); // for symmetry cases only
            else if (key == "theta")          p.theta      = std::stod(val);
            else if (key == "E")              p.E          = std::stod(val);
            else if (key == "nu")             p.nu         = std::stod(val);
            else if (key == "mu")             p.mu         = std::stod(val);
            else if (key == "lambda")         p.lam        = std::stod(val);
            else if (key == "friction")       p.friction   = std::stod(val);
            else if (key == "mesh")           p.mesh_path  = val;
            else if (key == "level")          p.mesh_level = std::stoul(val);
            else if (key == "tol_newton")     p.opts.tol_newton     = std::stod(val);
            else if (key == "tol_newton_rel") p.opts.tol_newton_rel = std::stod(val);
            else if (key == "tol_picard")     p.opts.tol_picard     = std::stod(val);
            else if (key == "max_newton")     p.opts.max_newton     = std::stoul(val);
            else if (key == "max_picard")     p.opts.max_picard     = std::stoul(val);
            else if (key == "verbose")        p.opts.verbose        = (val == "1" || val == "true");
            else
                std::cout << "Params: unknown key '" << key << "' ignored\n";
        }
        return p;
    }
};

} // namespace hho_contact
