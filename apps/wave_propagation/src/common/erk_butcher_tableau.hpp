//
//  erk_butcher_tableau.hpp
//  acoustics
//
//  Created by Omar Durán on 5/8/20.
//

#pragma once
#ifndef erk_butcher_tableau_hpp
#define erk_butcher_tableau_hpp

#include <Eigen/Dense>

class erk_butcher_tableau
{
    public:
    
    static void erk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        // Mathematical Aspects of Discontinuous Galerkin Methods
        // Authors: Di Pietro, Daniele Antonio, Ern, Alexandre
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        switch (s) {
            case 1:
                {
                    a(0,0) = 0.0;
                    b(0,0) = 1.0;
                    c(0,0) = 0.0;
                }
                break;
            case 2:
                {
//                    a(0,0) = 0.0;
//                    a(1,0) = 1.0;
//                    a(1,1) = 0.0;
//
//                    b(0,0) = 0.5;
//                    b(1,0) = 0.5;
//
//                    c(0,0) = 0.0;
//                    c(1,0) = 1.0;
                    
                    a(0,0) = 0.0;
                    a(1,0) = 0.5;
                    a(1,1) = 0.0;

                    b(0,0) = 0.0;
                    b(1,0) = 1.0;

                    c(0,0) = 0.0;
                    c(1,0) = 0.5;
                    
                }
                break;
            case 3:
                {

//                    a(0,0) = 0.0;
//                    a(1,0) = 1.0/3.0;
//                    a(1,1) = 0.0;
//                    a(2,0) = 0.0;
//                    a(2,1) = 2.0/3.0;
//                    a(2,2) = 0.0;
//
//                    b(0,0) = 1.0/4.0;
//                    b(1,0) = 0.0;
//                    b(2,0) = 3.0/4.0;
//
//                    c(0,0) = 0.0;
//                    c(1,0) = 1.0/3.0;
//                    c(2,0) = 2.0/3.0;
                    
                    a(0,0) = 0.0;
                    a(1,0) = 1.0/2.0;
                    a(1,1) = 0.0;
                    a(2,0) = -1.0;
                    a(2,1) = 2.0;
                    a(2,2) = 0.0;
                    
                    b(0,0) = 1.0/6.0;
                    b(1,0) = 2.0/3.0;
                    b(2,0) = 1.0/6.0;
                    
                    c(0,0) = 0.0;
                    c(1,0) = 1.0/2.0;
                    c(2,0) = 1.0;
                    
                }
                break;
                case 4:
                {

                    a(0,0) = 0.0;
                    a(1,0) = 1.0/2.0;
                    a(2,0) = 0.0;
                    a(3,0) = 0.0;
                    a(1,1) = 0.0;
                    a(2,1) = 1.0/2.0;
                    a(3,1) = 0.0;
                    a(2,2) = 0.0;
                    a(3,2) = 1.0;
                    a(3,3) = 0.0;
                    
                    b(0,0) = 1.0/6.0;
                    b(1,0) = 1.0/3.0;
                    b(2,0) = 1.0/3.0;
                    b(3,0) = 1.0/6.0;
                    
                    c(0,0) = 0.0;
                    c(1,0) = 1.0/2.0;
                    c(2,0) = 1.0/2.0;
                    c(3,0) = 1.0;
                    
                }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
    static void ssprk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
            
            // A New Class of Optimal High-Order Strong-Stability-Preserving Time Discretization Methods
            // Authors:  Raymond J. Spiteri and Steven J. Ruuth
            // Apendix B
            a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
            b = Matrix<double, Dynamic, 1>::Zero(s, 1);
            c = Matrix<double, Dynamic, 1>::Zero(s, 1);
            
            switch (s) {
                case 1: // Order 1
                    {
                        a(0,0) = 0.0;
                        b(0,0) = 1.0;
                        c(0,0) = 0.0;
                    }
                    break;
                case 2: // Order 1
                    {
                        a(0,0) = 0.0;
                        a(1,0) = 0.5;
                        a(1,1) = 0.0;
    
                        b(0,0) = 0.5;
                        b(1,0) = 0.5;
    
                        c(0,0) = 0.0;
                        c(1,0) = 0.5;
                        
                        
                    }
                    break;
                case 3: // Order 2
                    {

                        a(0,0) = 0.0;
                        a(1,0) = 0.5;
                        a(1,1) = 0.0;
                        a(2,0) = 0.5;
                        a(2,1) = 0.5;
                        a(2,2) = 0.0;
    
                        b(0,0) = 1.0/3.0;
                        b(1,0) = 1.0/3.0;
                        b(2,0) = 1.0/3.0;
    
                        c(0,0) = 0.0;
                        c(1,0) = 1.0/2.0;
                        c(2,0) = 1.0;
                        
                        
                    }
                    break;
                    case 4: // Order 3
                    {
                        
                        a(0,0) = 0.0;
                        a(1,0) = 1.0/2.0;
                        a(2,0) = 1.0/2.0;
                        a(3,0) = 1.0/6.0;
                        a(1,1) = 0.0;
                        a(2,1) = 1.0/2.0;
                        a(3,1) = 1.0/6.0;
                        a(2,2) = 0.0;
                        a(3,2) = 1.0/6.0;
                        a(3,3) = 0.0;
                        
                        b(0,0) = 1.0/6.0;
                        b(1,0) = 1.0/6.0;
                        b(2,0) = 1.0/6.0;
                        b(3,0) = 1.0/2.0;
                        
                        c(0,0) = 0.0;
                        c(1,0) = 1.0/2.0;
                        c(2,0) = 1.0;
                        c(3,0) = 1.0/2.0;
                        
                    }
                    break;
                    case 5: // Order 4
                    {
                        
                        a(0,0) = 0.0;
                        a(1,0) = 0.39175222700392;
                        a(2,0) = 0.21766909633821;
                        a(3,0) = 0.08269208670950;
                        a(4,0) = 0.06796628370320;
                        
                        a(1,1) = 0.0;
                        a(2,1) = 0.36841059262959;
                        a(3,1) = 0.13995850206999;
                        a(4,1) = 0.11503469844438;
                        
                        a(2,2) = 0.0;
                        a(3,2) = 0.2518917742473;
                        a(4,2) = 0.20703489864929;
                        
                        a(3,3) = 0.0;
                        a(4,3) = 0.54497475021237;
                        
                        b(0,0) = 0.14681187618661;
                        b(1,0) = 0.24848290924556;
                        b(2,0) = 0.10425883036650;
                        b(3,0) = 0.27443890091960;
                        b(4,0) = 0.22600748319395;
                        
                        c(0,0) = 0.0;
                        c(1,0) = 0.39175222700392;
                        c(2,0) = 0.58607968896779;
                        c(3,0) = 0.47454236302687;
                        c(4,0) = 0.93501063100924;
                        
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

#endif /* erk_butcher_tableau_hpp */
