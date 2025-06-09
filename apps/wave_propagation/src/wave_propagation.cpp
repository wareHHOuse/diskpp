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
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
using namespace Eigen;

#include "diskpp/common/timecounter.hpp"
#include "diskpp/methods/hho"
#include "diskpp/geometry/geometry.hpp"
#include "diskpp/boundary_conditions/boundary_conditions.hpp"
#include "diskpp/output/silo.hpp"

// application common sources
// #include "common/display_settings.hpp"
// #include "common/fitted_geometry_builders.hpp"
// #include "common/linear_solver.hpp"
// #include "common/acoustic_material_data.hpp"
// #include "common/elastic_material_data.hpp"
// #include "common/scal_vec_analytic_functions.hpp"
// #include "common/preprocessor.hpp"
// #include "common/postprocessor.hpp"

// RK schemes
// #include "../common/dirk_hho_scheme.hpp"
// #include "../common/dirk_butcher_tableau.hpp"
// #include "../common/erk_butcher_tableau.hpp"
// #include "../common/erk_hho_scheme.hpp"
// #include "../common/erk_coupling_hho_scheme.hpp"

// Prototypes:
// Computation of an empirical CFL criteria                    
// #include "Prototypes/CFL/EAcoustic_CFL.hpp"                          // CFl - Acoustic                      
// #include "Prototypes/CFL/EElasticity_CFL.hpp"                        // CFl - Linear Elasticity  
// #include "Prototypes/CFL/EHHOFirstOrderCFL.hpp"                      // CFl - Elasto-Acoustic Coupling 
// Stability study & Spectral radius computation:
// #include "Prototypes/Stability_Study/EAcoustic_stability.hpp"        // Acoustic
// #include "Prototypes/Stability_Study/EElastic_stability.hpp"         // Linear Elasticity
// #include "Prototypes/Stability_Study/EHHOFirstOrder_stability.hpp"   // Elasto-Acoustic Coupling                   
// Convergence test on sinusoidal analytical solution - Debug Implicit Direct Solver!!!!
// #include "Prototypes/Conv_Test/IAcoustic_conv_test.hpp"              // Implicit Acoustic               
// #include "Prototypes/Conv_Test/IElastic_conv_test.hpp"               // Implicit Elastic 
// #include "Prototypes/Conv_Test/IHHOFirstOrder.hpp"                   // Implicit Coupling                         
// #include "Prototypes/Conv_Test/IHHOFirstOrder_conv_tests.hpp"        // Explicit Coupling    
// #include "Prototypes/Conv_Test/EHHOFirstOrder.hpp"                   // Explicit Coupling    
// #include "Prototypes/Conv_Test/EHHOFirstOrder_conv_tests.hpp"        // Explicit Coupling    
// Pulses for comparison with Gar6more
// #include "Prototypes/Pulses/HeterogeneousIHHOFirstOrder.hpp"         // Implicit Pulse (adimensional)
// #include "Prototypes/Pulses/HeterogeneousEHHOFirstOrder.hpp"         // Explicit Pulse (adimensional)
// #include "Prototypes/Pulses/ConicWavesIHHOFirstOrder.hpp"            // Implicit Pulse (geophysic) 
// #include "Prototypes/Pulses/ConicWavesEHHOFirstOrder.hpp"            // Implicit Pulse (geophysic) 
// Sedimentary Basin
// #include "Prototypes/Basin/BassinIHHOFirstOrder.hpp"                 // Implicit Sedimentary Basin
// #include "Prototypes/Basin/Test_Laurent.hpp"                 // Implicit Sedimentary Basin
//  #include "Prototypes/Basin/BassinEHHOFirstOrder.hpp"              // Explicit Sedimentary Basin

int main(int argc, char **argv){

// CFL tables:
   // EAcoustic_CFL(argc, argv); 
   // EElasticity_CFL(argc, argv);
   // EHHOFirstOrderCFL(argc, argv); 

// Stability study & Spectral radius computation:
   // EAcoustic_stability(argc, argv);
   // EElastic_stability(argc, argv);
   // EHHOFirstOrder_stability(argc, argv); 

// Convergence test:
   // IAcoustic_conv_test(argc, argv);
   // IElastic_conv_test(argc, argv);
   // IHHOFirstOrder(argc, argv);
   // IHHOFirstOrder_conv_tests(argc, argv);
   // EHHOFirstOrder(argc, argv);
   // EHHOFirstOrder_conv_tests(argc, argv);

// Pulse:
   // HeterogeneousIHHOFirstOrder(argc, argv); 
   // HeterogeneousEHHOFirstOrder(argc, argv); 
   // ConicWavesIHHOFirstOrder(argc, argv);
   // ConicWavesEHHOFirstOrder(argc, argv);

// Bassin sédimentaire:
   // BassinIHHOFirstOrder(argc, argv);
   // Test_Laurent(argc, argv);
   // BassinEHHOFirstOrder(argc, argv); Not working 
  
}




