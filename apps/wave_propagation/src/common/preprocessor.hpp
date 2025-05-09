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
    
  size_t m_exact_functions;
    
  simulation_data() : m_k_degree(2), m_n_divs(2), m_hdg_stabilization_Q(false),
		      m_scaled_stabilization_Q(false), m_sc_Q(false), m_nt_divs(0),
		      m_render_silo_files_Q(true), m_report_energy_Q(false),
		      m_iterative_solver_Q(true),  m_polygonal_mesh_Q(false), m_exact_functions(1) {
    
  }
  
  simulation_data(const simulation_data & other){
    
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
    m_exact_functions        = other.m_exact_functions;
  }
  
  const simulation_data & operator=(const simulation_data & other){
    
    // check for self-assignment
    if(&other == this){
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
    m_exact_functions        = other.m_exact_functions;  
    return *this;
  }
  
  ~simulation_data(){  
  }
    
  void write_simulation_data(std::ostream & simulation_log = std::cout){
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

  void print_simulation_data(){
    
    std::cout << std::endl << "   ";
    std::cout << bold << red << "SIMULATION PARAMETERS : " << reset;
    
    std::cout << bold << cyan << std::endl;
    std::cout << "      ";
    std::cout << bold << "Iterative solver : " << m_iterative_solver_Q << reset;
    
    std::cout << bold << cyan << std::endl;
    std::cout << "      ";
    std::cout << bold << "HHO setting : " << reset;
    std::cout << yellow;
    std::cout << bold << " - Polynomial degree    -k : " << m_k_degree << "     (Face unknowns)"
	      << std::endl;
    std::cout << "                    ";
    std::cout << " - Stabilization type   -s : " << m_hdg_stabilization_Q
	      << "     (Equal-order=0, Mixed-order=1)" << std::endl;
    std::cout << "                    ";
    std::cout << " - Scaled stabilization -r : " << m_scaled_stabilization_Q
	      << "     (O(1)=0, O(1/h)=1)" << std::endl;
    std::cout << "                    ";
    std::cout << " - Static condensation  -c : " << m_sc_Q;
    
    std::cout << bold << cyan << std::endl;
    std::cout << "      ";
    std::cout << bold << "Mesh :" << reset;
    std::cout << yellow;
    std::cout << bold << "         - Polygonal mesh       -p : " << m_polygonal_mesh_Q << std::endl;
    std::cout << "                    ";
    std::cout << " - Refinement level     -l : " << m_n_divs << std::endl;
    
    std::cout << bold << cyan;
    std::cout << "      ";
    std::cout << bold << "Time refinement level -n : " << m_nt_divs << reset;
    
    std::cout << bold << cyan << std::endl;
    std::cout << "      ";
    std::cout << bold << "Post process :" << reset;
    std::cout << yellow;
    std::cout << bold << " - Silo files           -f : " << m_render_silo_files_Q << std::endl;
    std::cout << "                    ";
    std::cout << bold << " - Energy file          -e : " << m_report_energy_Q << reset;
    
    std::cout << std::endl << std::endl;
        
  }
    
};

class preprocessor {
    
public:
    
    

    static void PrintHelp()
    {
        std::cout <<
                "-k <int>:  Face polynomial degree: default 0\n"
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

    static simulation_data process_args(int argc, char** argv)
    {
        const char* const short_opts = "k:p:l:s:r:n:c:f:e:";
        const option long_opts[] = {
                {"degree", required_argument, nullptr, 'k'},
                {"poly_mesh", required_argument, nullptr, 'p'},
                {"xref", required_argument, nullptr, 'l'},
                {"stab", required_argument, nullptr, 's'},
                {"scal", required_argument, nullptr, 'r'},
                {"tref", required_argument, nullptr, 'n'},
                {"c", optional_argument, nullptr, 'c'},
                {"file", optional_argument, nullptr, 'f'},
                {"energy", optional_argument, nullptr, 'e'},
                {"help", no_argument, nullptr, 'h'},            
                {nullptr, no_argument, nullptr, 0}
        };

        size_t k_degree = 0;
        size_t n_divs = 0;
        size_t nt_divs = 0;
        size_t poly_mesh = 0;
        bool hdg_Q = false;
        bool scaled_Q = false;
        bool sc_Q = false;
        bool silo_files_Q = true;
        bool report_energy_Q = true;
        
        while (true)
        {
            const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

            if (-1 == opt)
                break;

            switch (opt)
            {
            case 'k':
                k_degree = std::stoi(optarg);
                break;

            case 'p':
                poly_mesh = std::stoi(optarg);
                break;

            case 'l':
                n_divs = std::stoi(optarg);
                break;

            case 's':
                hdg_Q = std::stoi(optarg);
                break;
                    
            case 'r':
                scaled_Q = std::stoi(optarg);
                break;
                    
            case 'n':
                nt_divs = std::stoi(optarg);
                break;
                    
            case 'c':
                sc_Q = std::stoi(optarg);
                break;
                    
            case 'f':
                silo_files_Q = std::stoi(optarg);
                break;

            case 'e':
                report_energy_Q = std::stoi(optarg);
                break;
                    

            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                preprocessor::PrintHelp();
                break;
            }
        }
        
        // populating simulation data
        simulation_data sim_data;
        sim_data.m_k_degree = k_degree;
        sim_data.m_polygonal_mesh_Q = poly_mesh;
        sim_data.m_n_divs = n_divs;
        sim_data.m_hdg_stabilization_Q = hdg_Q;
        sim_data.m_scaled_stabilization_Q = scaled_Q;
        sim_data.m_sc_Q = sc_Q;
        sim_data.m_nt_divs = nt_divs;
        sim_data.m_render_silo_files_Q = silo_files_Q;
        sim_data.m_report_energy_Q = report_energy_Q;
        return sim_data;
    }
    
    static void PrintTestHelp()
    {
        std::cout <<
                "-k <int>:  Maximum Face polynomial degree: default 0\n"
                "-l <int>:  Maximum Number of uniform space refinements: default 0\n"
                "-s <0-1>:  Stabilization type 0 -> HHO, 1 -> HDG-like: default 0 \n"
                "-r <0-1>:  Scaled stabilization 0 -> without, 1 -> with: default 0 \n"
                "-q <0-1>:  Quadratic function type 0 -> non-polynomial, 1 -> quadratic: default 0 \n"
                "-f <0-1>:  Write silo files : default 0\n"
                "-c <0-1>:  Static Condensation: default 0 \n"
                "-help:     Show help\n";
        exit(1);
    }
    
    static simulation_data process_convergence_test_args(int argc, char** argv)
    {
        const char* const short_opts = "k:l:s:r:c:q:f:";
        const option long_opts[] = {
                {"degree", required_argument, nullptr, 'k'},
                {"xref", required_argument, nullptr, 'l'},
                {"stab", required_argument, nullptr, 's'},
                {"scal", required_argument, nullptr, 'r'},
                {"file", optional_argument, nullptr, 'f'},
                {"qfunc", optional_argument, nullptr, 'q'},
                {"cond", optional_argument, nullptr, 'c'},
                {"help", no_argument, nullptr, 'h'},
                {nullptr, no_argument, nullptr, 0}
        };

        size_t k_degree = 0;
        size_t n_divs   = 0;
        bool hdg_Q = false;
        bool scaled_Q = false;
        bool sc_Q = false;
        bool quadratic_func_Q = false;
        bool silo_files_Q = false;
        
        while (true)
        {
            const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

            if (-1 == opt)
                break;

            switch (opt)
            {
            case 'k':
                k_degree = std::stoi(optarg);
                break;

            case 'l':
                n_divs = std::stoi(optarg);
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

            case 'q':
                    quadratic_func_Q = std::stoi(optarg);
                break;
                    
            case 'f':
                silo_files_Q = std::stoi(optarg);
                break;

            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                preprocessor::PrintTestHelp();
                break;
            }
        }
        
        // populating simulation data
        simulation_data sim_data;
        sim_data.m_k_degree = k_degree;
        sim_data.m_n_divs = n_divs;
        sim_data.m_hdg_stabilization_Q = hdg_Q;
        sim_data.m_scaled_stabilization_Q = scaled_Q;
        sim_data.m_sc_Q = sc_Q;
        sim_data.m_render_silo_files_Q = silo_files_Q;
        return sim_data;
    }
    
    static simulation_data process_convergence_test_lua_file(char** argv)
    {
        
        sol::state lua;
        lua.open_libraries(sol::lib::math, sol::lib::base);
        
        // expected input tables
        lua["config"] = lua.create_table();

        std::string lua_file = argv[3];
        auto r = lua.do_file(lua_file);
        if ( !r.valid() ){
            PrintConvergeceTestExample();
        }
        
        // populating simulation data
        simulation_data sim_data;
        sim_data.m_k_degree                 = lua["config"]["max_k_deg"].get_or(0);
        sim_data.m_n_divs                   = lua["config"]["max_l_ref"].get_or(0);
        sim_data.m_hdg_stabilization_Q      = lua["config"]["stab_type"].get_or(0);
        sim_data.m_scaled_stabilization_Q   = lua["config"]["stab_scal"].get_or(0);
        sim_data.m_sc_Q                     = lua["config"]["stat_cond"].get_or(0);
        sim_data.m_render_silo_files_Q      = lua["config"]["silo_output"].get_or(0);
        sim_data.m_iterative_solver_Q       = lua["config"]["iter_solv"].get_or(0);
        sim_data.m_exact_functions          = lua["config"]["exac_func"].get_or(0);
        sim_data.m_polygonal_mesh_Q         = lua["config"]["poly_mesh"].get_or(0);
        return sim_data;
    }
    
    static void PrintConvergeceTestExample()
    {
        std::cout << "Please specify a lua configuration file like this: " << std::endl;
        std::cout <<
        "config.max_k_deg = 4 -- <int>:  Maximum face polynomial degree: default 0\n"
        "config.max_l_ref = 4 -- <int>:  Maximum number of uniform spatial refinements: default 0\n"
        "config.stab_type = 0 -- <0-1>:  Stabilization type 0 (HHO), 1 (HDG-like): default 0\n"
        "config.stab_scal = 1 -- <0-1>:  Stabilization scaling 0 O(1), 1 O(1/h_{f}): default 0\n"
        "config.stat_cond = 1 -- <0-1>:  Static condensation: default 0\n"
        "config.iter_solv = 0 -- <0-1>:  Iterative solver : default 0\n"
        "config.exac_func = 0 -- <0-1>:  Manufactured function type 0 (non-polynomial), 1 (quadratic): default 0\n"
        "config.poly_mesh = 0 -- <0-1>:  Use of polynoal meshes : default 0\n"
        "config.silo_output = 0 -- <0-1>:  Write silo files : default 0\n";
        exit(1);
        throw std::invalid_argument("Program will stop.");
    }
    
    static simulation_data process_acoustics_lua_file(char** argv)
    {
        
        sol::state lua;
        lua.open_libraries(sol::lib::math, sol::lib::base);
        
        // expected input tables
        lua["config"] = lua.create_table();

        std::string lua_file = argv[3];
        auto r = lua.do_file(lua_file);
        if ( !r.valid() ){
            PrintSimulationExample();
        }
        
        // populating simulation data
        simulation_data sim_data;
        sim_data.m_k_degree                 = lua["config"]["fac_k_deg"].get_or(0);
        sim_data.m_n_divs                   = lua["config"]["num_l_ref"].get_or(0);
        sim_data.m_nt_divs                  = lua["config"]["num_t_ref"].get_or(0);
        sim_data.m_hdg_stabilization_Q      = lua["config"]["stab_type"].get_or(0);
        sim_data.m_scaled_stabilization_Q   = lua["config"]["stab_scal"].get_or(0);
        sim_data.m_sc_Q                     = lua["config"]["stat_cond"].get_or(0);
        sim_data.m_render_silo_files_Q      = lua["config"]["silo_output"].get_or(0);
        sim_data.m_iterative_solver_Q       = lua["config"]["iter_solv"].get_or(0);
        sim_data.m_polygonal_mesh_Q         = lua["config"]["poly_mesh"].get_or(0);
        sim_data.m_exact_functions          = lua["config"]["exac_func"].get_or(0);
        sim_data.m_report_energy_Q          = lua["config"]["writ_ener"].get_or(0);
        return sim_data;
    }
    
    static void PrintSimulationExample()
    {
        std::cout << "Please specify a lua configuration file like this: " << std::endl;
        std::cout <<
        "config.fac_k_deg = 3 -- <int>:  Face polynomial degree: default 0\n"
        "config.num_l_ref = 3 -- <int>:  Number of uniform spatial refinements: default 0\n"
        "config.num_t_ref = 7 -- <int>:  Number of uniform time refinements: default 0\n"
        "config.stab_type = 1 -- <0-1>:  Stabilization type 0 (HHO), 1 (HDG-like): default 0\n"
        "config.stab_scal = 0 -- <0-1>:  Stabilization scaling 0 O(1), 1 O(1/h_{f}): default 0\n"
        "config.stat_cond = 1 -- <0-1>:  Static condensation: default 0\n"
        "config.iter_solv = 0 -- <0-1>:  Iterative solver : default 0\n"
        "config.poly_mesh = 0 -- <0-1>:  Use of polynoal meshes : default 0\n"
        "config.exac_func = 0 -- <0-2>:  Exact function type {0,1,2} : default 0\n"
        "config.writ_ener = 1 -- <0-1>:  Report energy at each time value : default 0\n"
        "config.silo_output = 0 -- <0-1>:  Write silo files : default 0\n";
        exit(1);
        throw std::invalid_argument("Program will stop.");
    }
    
};

#endif /* preprocessor_hpp */
