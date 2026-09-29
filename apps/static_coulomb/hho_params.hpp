#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace hho_contact {

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
    double gamma_0 = 10.0;           // gamma_n = gamma_t = gamma_0 / h_F
    double theta   = -1.0;

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
            lam = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu));   // plane strain / 3D
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
            else if (key == "gamma_0")        p.gamma_0    = std::stod(val);
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
