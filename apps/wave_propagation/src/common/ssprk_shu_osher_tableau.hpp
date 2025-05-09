//
//  ssprk_shu_osher_tableau.hpp
//  acoustics
//
//  Created by Omar Durán on 4/21/20.
//

#pragma once
#ifndef ssprk_shu_osher_tableau_hpp
#define ssprk_shu_osher_tableau_hpp

#include <Eigen/Dense>

class ssprk_shu_osher_tableau
{
    public:
    
    static void ossprk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double,Dynamic, Dynamic> &b){
        
        // Optimal Strong-Stability-Preserving Runge-Kutta Time Discretizations for Discontinuous Galerkin Methods
        switch (s) {
            case 1:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    b(0,0) = 1.0;
                }
                break;
            case 2:
                {
                    // DG book (Alexandre)
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,0) = 1.0/2.0;
                    a(1,1) = 1.0/2.0;

                    b(0,0) = 1.0;
                    b(1,1) = 1.0/2.0;
                                        
                }
                break;
            case 3:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,0) = 0.087353119859156;
                    a(1,1) = 0.912646880140844;
                    a(2,0) = 0.344956917166841;
                    a(2,1) = 0.0;
                    a(2,2) = 0.655043082833159;

                    b(0,0) = 0.528005024856522;
                    b(1,0) = 0.0;
                    b(1,1) = 0.481882138633993;
                    b(2,0) = 0.022826837460491;
                    b(2,1) = 0.0;
                    b(2,2) = 0.345866039233415;
        
                    
                }
                break;
            case 4:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,0) = 0.522361915162541;
                    a(1,1) = 0.477638084837459;
                    a(2,0) = 0.368530939472566;
                    a(2,1) = 0.0;
                    a(2,2) = 0.631469060527434;
                    a(3,0) = 0.334082932462285;
                    a(3,1) = 0.006966183666289;
                    a(3,2) = 0.0;
                    a(3,3) = 0.658950883871426;
                    
                    b(0,0) = 0.594057152884440;
                    b(1,0) = 0.0;
                    b(1,1) = 0.283744320787718;
                    b(2,0) = 0.000000038023030;
                    b(2,1) = 0.0;
                    b(2,2) = 0.375128712231540;
                    b(3,0) = 0.116941419604231;
                    b(3,1) = 0.004138311235266;
                    b(3,2) = 0.0;
                    b(3,3) = 0.391454485963345;
                }
                break;
            case 5:
            {
                a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                a(0,0) = 1.0;
                a(1,0) = 0.261216512493821;
                a(1,1) = 0.738783487506179;
                a(2,0) = 0.623613752757655;
                a(2,1) = 0.0;
                a(2,2) = 0.376386247242345;
                a(3,0) = 0.444745181201454;
                a(3,1) = 0.120932584902288;
                a(3,2) = 0.0;
                a(3,3) = 0.434322233896258;
                a(4,0) = 0.213357715199957;
                a(4,1) = 0.209928473023448;
                a(4,2) = 0.063353148180384;
                a(4,3) = 0.0;
                a(4,4) = 0.513360663596212;
                
                b(0,0) = 0.605491839566400;
                b(1,0) = 0.0;
                b(1,1) = 0.447327372891397;
                b(2,0) = 0.000000844149769;
                b(2,1) = 0.0;
                b(2,2) = 0.227898801230261;
                b(3,0) = 0.002856233144485;
                b(3,1) = 0.073223693296006;
                b(3,2) = 0.0;
                b(3,3) = 0.262978568366434;
                b(4,0) = 0.002362549760441;
                b(4,1) = 0.127109977308333;
                b(4,2) = 0.038359814234063;
                b(4,3) = 0.0;
                b(4,4) = 0.310835692561898;
            }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
    static void new_ossprk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double,Dynamic, Dynamic> &b){
        
        // NEW OPTIMAL SSPRK TIME DISCRETIZATION M (s,s)
        switch (s) {
            case 1:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    b(0,0) = 1.0;
                }
                break;
            case 2:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,0) = 1.0/2.0;
                    a(1,1) = 1.0/2.0;

                    b(0,0) = 1.0;
                    b(1,0) = 0.0;
                    b(1,1) = 1.0/2.0;
                                        
                }
                break;
            case 3:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    
                    a(0,0) = 1.0;
                    a(1,0) = 3.0/4.0;
                    a(1,1) = 1.0/4.0;
                    a(2,0) = 1.0/3.0;
                    a(2,1) = 0.0;
                    a(2,2) = 2.0/3.0;

                    b(0,0) = 1.0;
                    b(1,0) = 0.0;
                    b(1,1) = 1.0/4.0;
                    b(2,0) = 0.0;
                    b(2,1) = 0.0;
                    b(2,2) = 2.0/3.0;
                    
                }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
    static void new_ossprk_s_p_one_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double,Dynamic, Dynamic> &b){
        
        // NEW OPTIMAL SSPRK TIME DISCRETIZATION M (s+1,s)
        switch (s) {
            case 1:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    b(0,0) = 1.0;
                }
                break;
            case 2:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,1) = 1.0;

                    b(0,0) = 1.0/2.0;
                    b(1,1) = 1.0/2.0;
                                        
                }
                break;
            case 3:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    
                    a(0,0) = 1.0;
                    a(1,0) = 0.0;
                    a(1,1) = 1.0;
                    a(2,0) = 1.0/3.0;
                    a(2,1) = 0.0;
                    a(2,2) = 2.0/3.0;

                    b(0,0) = 1.0/2.0;
                    b(1,1) = 1.0/2.0;
                    b(2,2) = 1.0/3.0;
        
                    
                }
                break;
            case 4:
                {
                    a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                    a(0,0) = 1.0;
                    a(1,0) = 0.0;
                    a(1,1) = 1.0;
                    a(2,0) = 2.0/3.0;
                    a(2,1) = 0.0;
                    a(2,2) = 1.0/3.0;
                    a(3,0) = 0.0;
                    a(3,1) = 0.0;
                    a(3,2) = 0.0;
                    a(3,3) = 1.0;
                    
                    b(0,0) = 1.0/2.0;
                    b(1,0) = 0.0;
                    b(1,1) = 1.0/2.0;
                    b(2,0) = 0.0;
                    b(2,1) = 0.0;
                    b(2,2) = 1.0/6.0;
                    b(3,0) = 0.0;
                    b(3,1) = 0.0;
                    b(3,2) = 0.0;
                    b(3,3) = 1.0/2.0;
                }
                break;
            case 5:
            {
                a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                b = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
                a(0,0) = 1.0;
                a(1,0) = 0.44437049406734;
                a(1,1) = 0.55562950593266;
                a(2,0) = 0.62010185138540;
                a(2,1) = 0.0;
                a(2,2) = 0.37989814861460;
                a(3,0) = 0.17807995410773;
                a(3,1) = 0.0;
                a(3,2) = 0.0;
                a(3,3) = 0.82192004589227;
                a(4,0) = 0.00683325884039;
                a(4,1) = 0.0;
                a(4,2) = 0.51723167208978;
                a(4,3) = 0.12759831133288;
                a(4,4) = 0.34833675773694;
                
                b(0,0) = 0.39175222700392;
                b(1,0) = 0.0;
                b(1,1) = 0.36841059262959;
                b(2,0) = 0.0;
                b(2,1) = 0.0;
                b(2,2) = 0.25189177424738;
                b(3,0) = 0.0;
                b(3,1) = 0.0;
                b(3,2) = 0.0;
                b(3,3) = 0.54497475021237;
                b(4,0) = 0.0;
                b(4,1) = 0.0;
                b(4,2) = 0.0;
                b(4,3) = 0.08460416338212;
                b(4,4) = 0.22600748319395;
            }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
  
};

#endif /* ssprk_shu_osher_tableau_hpp */
