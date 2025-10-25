
//
//  preprocessor.hpp
//  acoustics
//
//  Created by Omar Durán on 4/10/20.
//

#pragma once
#ifndef preprocessor_hpp
#define preprocessor_hpp

#include <getopt.h>
#include "../../../../contrib/sol2/include/sol/sol.hpp"


class simulation_data {
  
public:
  
    size_t m_k_degree;
    size_t m_n_divs;
    bool m_hdg_stabilization_Q;
    bool m_scaled_stabilization_Q;
    bool m_sc_Q;
    size_t m_nt_divs;
    bool m_render_silo_files_Q;
    bool m_report_energy_Q;
    bool m_iterative_solver_Q;
    bool m_polygonal_mesh_Q;
  
    simulation_data() : m_k_degree(2), m_n_divs(2), m_hdg_stabilization_Q(false), m_scaled_stabilization_Q(false), m_sc_Q(false), m_nt_divs(0), m_render_silo_files_Q(true), m_report_energy_Q(false), m_iterative_solver_Q(false),  m_polygonal_mesh_Q(false) {}
  
    simulation_data(const simulation_data & other) {
        m_k_degree               = other.m_k_degree;
        m_n_divs                 = other.m_n_divs;
        m_hdg_stabilization_Q    = other.m_hdg_stabilization_Q;
        m_scaled_stabilization_Q = other.m_scaled_stabilization_Q;
        m_sc_Q                   = other.m_sc_Q;
        m_nt_divs                = other.m_nt_divs;
        m_render_silo_files_Q    = other.m_render_silo_files_Q;
        m_report_energy_Q        = other.m_report_energy_Q;
        m_iterative_solver_Q     = other.m_iterative_solver_Q;
        m_polygonal_mesh_Q       = other.m_polygonal_mesh_Q;
    }
  
    const simulation_data & operator=(const simulation_data & other) {
    
        if (&other == this) {
            return *this;
        }
        
        m_k_degree               = other.m_k_degree;
        m_n_divs                 = other.m_n_divs;
        m_hdg_stabilization_Q    = other.m_hdg_stabilization_Q;
        m_scaled_stabilization_Q = other.m_scaled_stabilization_Q;
        m_sc_Q                   = other.m_sc_Q;
        m_nt_divs                = other.m_nt_divs;
        m_render_silo_files_Q    = other.m_render_silo_files_Q;
        m_report_energy_Q        = other.m_report_energy_Q;
        m_iterative_solver_Q     = other.m_iterative_solver_Q;
        m_polygonal_mesh_Q       = other.m_polygonal_mesh_Q;
        
        return *this;
    
    }

    ~simulation_data(){}
    
    void write_simulation_data(std::ostream & simulation_log = std::cout) {
        simulation_log << "  Input parameters:" << std::endl;
        simulation_log << "    Iterative solver: " << m_iterative_solver_Q << std::endl;
        simulation_log << "    HHO setting: " << std::endl;
        simulation_log << "      - Polynomial degree    -k: " << m_k_degree << "  (Face unknowns)" << std::endl;
        simulation_log << "      - Stabilization type   -s: " << m_hdg_stabilization_Q << "  (Equal-order=0, Mixed-order=1)" << std::endl;
        simulation_log << "      - Scaled stabilization -r: " << m_scaled_stabilization_Q << "  (O(1)=0, O(1/h)=1)" << std::endl;
        simulation_log << "      - Static condensation  -c: " << m_sc_Q << std::endl;
        simulation_log << "    Mesh:" << "- Refinement level     -l : " << m_n_divs << std::endl;
        simulation_log << "    Time refinement level -n : " << m_nt_divs << std::endl;
    };

    void print_simulation_data() {
        std::cout << "\n   " << bold << red << "HHO SETTING: " << reset;
        std::cout << bold << cyan
                  << "\n      " << bold << "Discretization: ";
                  if (m_hdg_stabilization_Q) {
                      std::cout << bold << yellow << "Mixed-order" << "     -s"
                                << "\n                      Cell degree = " << m_k_degree+1 << " -k"
                                << "\n                      Face degree = " << m_k_degree << std::endl;
                  }
                  else {
                      std::cout << bold << yellow << "Equal-order" << "     -s"
                                << "\n                      Cell degree = " << m_k_degree << " -k"
                                << "\n                      Face degree = " << m_k_degree << std::endl;
                  }
                  std::cout << bold << cyan << "      Stabilization scaling: ";
                  if (m_scaled_stabilization_Q) {
                      std::cout << "O(1/h)   -r" << std::endl;
                  }
                  else {
                      std::cout << "O(1)     -r " << std::endl;
                  }
                  std::cout << cyan << "      Static condensation = ";
                  if (m_sc_Q) {
                      std::cout << "YES" << "       -c" << reset << std::endl;
                  }
                  else {
                      std::cout << "NO" << "        -c" << reset << std::endl;
                  }
                  
                  std::cout << "\n   " << bold << red << "SIMULATION PARAMETERS: "                   
                  << "\n      " << bold << cyan << "Mesh:" << yellow << bold << " Cartesian mesh: ";
                  if (m_polygonal_mesh_Q) {
                      std::cout << "NO" << "        -p";
                  }
                  else {
                      std::cout << "YES" << "       -p";
                  }
                  std::cout << "\n            " << "Mesh refinement level: " << m_n_divs << "  -l";
                  if (m_nt_divs > 12) {
                      std::cout << "\n      " << bold << cyan << "Number of time steps: " << m_nt_divs << "      -n";
                  }
                  else {
                      std::cout << "\n      " << bold << cyan << "Time refinement level: " << m_nt_divs << "        -n";
                  }
                  std::cout << "\n      " << bold << cyan << "Solver:";
                  if (m_iterative_solver_Q) {
                      std::cout << " Iterative" << "               -i" << std::endl;
                  }
                  else {
                      std::cout << " Direct" << "                  -i" << std::endl;
                  }

        std::cout << "\n   " << bold << red << "POSTPROCESS:" << reset
                  << cyan
                  << "\n        " << bold << "Silo files: ";
                  if (m_render_silo_files_Q) {
                      std::cout << " YES" << " -f";
                  }
                  else {
                      std::cout << " NO" << "  -f";
                  }
                  std::cout << "\n        " << bold << "Energy file: ";
                  if (m_report_energy_Q) {
                        std::cout << "YES" << " -e";
                    }
                    else {
                        std::cout << "NO" << "  -e";
                  }
    }          
    
};

class preprocessor {
    
public:
    
    static void PrintHelp() {

        std::cout <<
                "-k <int>:  blabla Face polynomial degree: default 0\n"
                "-l <int>:  Number of uniform space refinements: default 0\n"
                "-s <0-1>:  Stabilization type 0 -> HHO, 1 -> HDG-like: default 0 \n"
                "-r <0-1>:  Scaled stabilization 0 -> without, 1 -> with: default 0 \n"
                "-n <int>:  Number of uniform time refinements: default 0\n"
                "-f <0-1>:  Write silo files: default 0\n"
                "-e <0-1>:  Report (time,energy) pairs: default 0\n"
                "-c <0-1>:  Static Condensation (implicit schemes): default 0 \n"
                "-help:     Show help\n";
        exit(1);
    }

    static simulation_data process_args(int argc, char** argv) {

        const char* const short_opts = "k:s:r:c:p:l:n:i:f:e:?";
        
        const option long_opts[] = {
            // HHO SETTING
            {"degree",    required_argument, nullptr, 'k'},
            {"stab",      required_argument, nullptr, 's'},
            {"scal",      required_argument, nullptr, 'r'},
            {"c",         optional_argument, nullptr, 'c'},
            // SIMULATION PARAMETERS
            {"poly_mesh", required_argument, nullptr, 'p'},
            {"xref",      required_argument, nullptr, 'l'},
            {"tref",      required_argument, nullptr, 'n'},
            {"solv",      required_argument, nullptr, 'i'},
            // POST PROCESSOR
            {"file",      optional_argument, nullptr, 'f'},
            {"energy",    optional_argument, nullptr, 'e'},
            // HELP
            {"help",      no_argument, nullptr,       '?'},            
            {nullptr,     no_argument, nullptr,         0}
        };

        size_t k_degree      = 0;
        size_t n_divs        = 0;
        size_t nt_divs       = 0;
        size_t poly_mesh     = 0;
        bool hdg_Q           = false;
        bool scaled_Q        = false;
        bool sc_Q            = false;
        bool silo_files_Q    = true;
        bool report_energy_Q = true;
        bool it_sol          = false;
        
        while (true) {
            const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
            if (-1 == opt)
                break;
            switch (opt) {
                // HHO SETTING
                case 'k':
                    k_degree = std::stoi(optarg);
                    break;
                case 's':
                    hdg_Q = std::stoi(optarg);
                    break;
                case 'r':
                    scaled_Q = std::stoi(optarg);
                    break;
                case 'c':
                    sc_Q = std::stoi(optarg);
                    break;
                // SIMULATION PARAMETERS
                case 'p':
                    poly_mesh = std::stoi(optarg);
                    break;
                case 'l':
                    n_divs = std::stoi(optarg);
                    break;
                case 'n':
                    nt_divs = std::stoi(optarg);
                    break;
                case 'i':
                    it_sol = std::stoi(optarg);
                    break;
                // POSTPROCESSOR
                case 'f':
                    silo_files_Q = std::stoi(optarg);
                    break;
                case 'e':
                    report_energy_Q = std::stoi(optarg);
                    break;
                // HELP
                case '?': 
                    preprocessor::PrintHelp();
                    break;
            }
        }
        
        // POPULATING SIM DATA
        simulation_data sim_data;
            // HHO SETTING
            sim_data.m_k_degree               = k_degree;
            sim_data.m_hdg_stabilization_Q    = hdg_Q;
            sim_data.m_scaled_stabilization_Q = scaled_Q;
            sim_data.m_sc_Q                   = sc_Q;
            // SIMULATION PARAMETERS
            sim_data.m_polygonal_mesh_Q       = poly_mesh;
            sim_data.m_n_divs                 = n_divs;
            sim_data.m_nt_divs                = nt_divs;
            sim_data.m_iterative_solver_Q     = it_sol;
            // POSTPROCESSOR
            sim_data.m_render_silo_files_Q    = silo_files_Q;
            sim_data.m_report_energy_Q        = report_energy_Q;
            
            return sim_data;
    }
    
    
};

#endif /* preprocessor_hpp */
