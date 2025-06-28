
//  Contributions by Omar Durán and Romain Mottier

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <list>
#include <getopt.h>


#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
using namespace Eigen;

#include "timecounter.h"
#include "methods/hho"
#include "geometry/geometry.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "output/silo.hpp"

// application common sources
#include "../common/display_settings.hpp"
#include "../common/fitted_geometry_builders.hpp"
#include "../common/linear_solver.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/scal_analytic_functions.hpp"
#include "../common/preprocessor.hpp"
#include "../common/postprocessor.hpp"

// implicit RK schemes
#include "../common/dirk_hho_scheme.hpp"
#include "../common/dirk_butcher_tableau.hpp"

// explicit RK schemes
#include "../common/erk_hho_scheme.hpp"
#include "../common/erk_butcher_tableau.hpp"
#include "../common/ssprk_hho_scheme.hpp"
#include "../common/ssprk_shu_osher_tableau.hpp"

/// Enumerate defining available function prototypes
enum EAcousticPrototype {
    OneFieldConvTest = 0,
    TwoFieldsConvTest = 1,
    IOneFieldAcoustic = 2,
    ITwoFieldsAcoustic = 3,
    ETwoFieldsAcoustic = 4,
    Hete1DIOneFieldAcoustic = 5,
    Hete1DITwoFieldsAcoustic = 6,
    Hete1DETwoFieldsAcoustic = 7,
    Hete2DIOneFieldAcoustic = 8,
    Hete2DITwoFieldsAcoustic = 9,
    Hete2DETwoFieldsAcoustic = 10
};

// ----- common data types ------------------------------
using RealType = double;
typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
typedef disk::BoundaryConditions<mesh_type, true> boundary_type;

//Prototype sources
#include "Prototypes/Elliptic_Conv_Test/EllipticOneFieldConvergenceTest.hpp"
#include "Prototypes/Elliptic_Conv_Test/EllipticTwoFieldsConvergenceTest.hpp"
#include "Prototypes/Acoustic_Conv_Test/IHHOSecondOrder.hpp"
#include "Prototypes/Acoustic_Conv_Test/IHHOFirstOrder.hpp"
#include "Prototypes/Acoustic_Conv_Test/EHHOFirstOrder.hpp"
#include "Prototypes/Heterogeneous/HeterogeneousIHHOSecondOrder.hpp"
#include "Prototypes/Heterogeneous/HeterogeneousIHHOFirstOrder.hpp"
#include "Prototypes/Heterogeneous/HeterogeneousEHHOFirstOrder.hpp"
#include "Prototypes/Heterogeneous_Pulse/HeterogeneousPulseIHHOSecondOrder.hpp"
#include "Prototypes/Heterogeneous_Pulse/HeterogeneousPulseIHHOFirstOrder.hpp"
//#include "Prototypes/Heterogeneous_Pulse/HeterogeneousPulseIHHOFirstOrder2.hpp"
#include "Prototypes/Heterogeneous_Pulse/HeterogeneousPulseEHHOFirstOrder.hpp"
#include "Prototypes/Gar6more_2D/HeterogeneousGar6more2DIHHOSecondOrder.hpp"
#include "Prototypes/Gar6more_2D/HeterogeneousGar6more2DIHHOFirstOrder.hpp"

#include "Prototypes/prototype_selector.hpp"

// ----- function prototypes ------------------------------

// Convergece tests for elliptic problem
void EllipticOneFieldConvergenceTest(char **argv);
void EllipticTwoFieldsConvergenceTest(char **argv);

// Convergece tests for the acoustic wave equation
void IHHOSecondOrder(char **argv);
void IHHOFirstOrder(char **argv);
void EHHOFirstOrder(char **argv);

// Simulation for a heterogeneous acoustics problem
void HeterogeneousIHHOSecondOrder(char **argv);
void HeterogeneousIHHOFirstOrder(char **argv);
void HeterogeneousEHHOFirstOrder(char **argv);

// Simulation for a heterogeneous acoustics problem on polygonal meshes
void HeterogeneousPulseIHHOSecondOrder(char **argv);
void HeterogeneousPulseIHHOFirstOrder(char **argv);
void HeterogeneousPulseEHHOFirstOrder(char **argv);

// Comparison with Gar6more2D
void HeterogeneousGar6more2DIHHOSecondOrder(int argc, char **argv);
void HeterogeneousGar6more2DIHHOFirstOrder(int argc, char **argv);

int prototype_selector(char **argv, EAcousticPrototype prototype);

void print_prototype_description();

int main(int argc, char **argv) {
  
  HeterogeneousGar6more2DIHHOFirstOrder(argc,argv);

  //EAcousticPrototype prototype;
  //const char* const short_opts = "p:h:";
  //const option long_opts[] = {
  //  {"prototype", required_argument, nullptr, 'p'},
  //  {"help", no_argument, nullptr, 'h'},
  //  {nullptr, no_argument, nullptr, 0}
  //};
  //
  //while (true) {
  //  const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
  //
  //  if (-1 == opt)
  //    break;
  //  
  //  switch (opt) {
  //  case 'p':
  //    prototype = EAcousticPrototype(std::stoi(optarg));
  //    break;
  //  case 'h':
  //    print_prototype_description();
  //    break;
  //  case '?':
  //  default:
  //    throw std::invalid_argument("Invalid option.");
  //    break;
  //  }
  //}
  //
  //if (argc != 4) {
  //  throw std::invalid_argument("Please specify and option and a lua configuration file.");
  //  return 0;
  //}
  //
  //HeterogeneousGar6more2DIHHOFirstOrder(argc,argv);
    
  //return prototype_selector(argv, prototype);
    
}

