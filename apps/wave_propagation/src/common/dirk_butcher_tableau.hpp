//
//  dirk_butcher_tableau.hpp
//  acoustics
//
//  Created by Omar Durán on 4/21/20.
//

#pragma once
#ifndef dirk_butcher_tableau_hpp
#define dirk_butcher_tableau_hpp

#include <Eigen/Dense>

class dirk_butcher_tableau
{
    public:
    
    static void dirk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        switch (s) { // Implicit midpoint DIRK(1, 2) or SDIRK(1, 2)
            case 1:
                {
                    a(0,0) = 0.5;
                    b(0,0) = 1.0;
                    c(0,0) = 0.5;
                }
                break;
            case 2:
                {
                    // https://www.jstor.org/stable/pdf/43692558.pdf?refreqid=excelsior%3A6b91257ea1427dcb8b00a78353c96ef4
                    // Construction Of High-order symplectic Runge-Kutta Methods
                    // SDIRK(2, 2)
                    a(0,0) = 0.25;
                    a(1,0) = 0.5;
                    a(1,1) = 0.25;

                    b(0,0) = 0.5;
                    b(1,0) = 0.5;

                    c(0,0) = 1.0/4.0;
                    c(1,0) = 3.0/4.0;
                    
                }
                break;
            case 3:
                {
//                    // DIRK(3, 4)
                    a(0,0) = 0.956262348020;
                    a(1,0) = -0.676995728936;
                    a(1,1) = 1.092920059741;
                    a(2,0) = 4.171447220367;
                    a(2,1) = -5.550750999686;
                    a(2,2) = 1.189651889660;
                    
                    b(0,0) = 0.228230955547;
                    b(1,0) = 0.706961029433;
                    b(2,0) = 0.064808015020;
                    
                    c(0,0) = 0.956262348020;
                    c(1,0) = 0.415924330804;
                    c(2,0) = -0.189651889660;
                    
                }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
    static void odirk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        // Optimized diagonally implicit Runge-Kutta schemes for time-dependent wave propagation problems
        switch (s) {
            case 1:
                {
                    a(0,0) = 0.5;
                    b(0,0) = 1.0;
                    c(0,0) = 0.5;
                }
                break;
            case 2:
                {
                    a(0,0) = 0.780078125000;
                    a(1,0) = -0.595072059507;
                    a(1,1) = 0.797536029754;

                    b(0,0) = 0.515112081837;
                    b(1,0) = 0.484887918163;

                    c(0,0) = 0.780078125;
                    c(1,0) = 0.202463970246;
                    

                    
                }
                break;
            case 3:
                {
                    
                    a(0,0) = 0.956262348020;
                    a(1,0) = -0.676995728936;
                    a(1,1) = 1.092920059741;
                    a(2,0) = 4.171447220367;
                    a(2,1) = -5.550750999686;
                    a(2,2) = 1.189651889660;
                    
                    b(0,0) = 0.228230955547;
                    b(1,0) = 0.706961029433;
                    b(2,0) = 0.064808015020;
                    
                    c(0,0) = 0.956262348020;
                    c(1,0) = 0.415924330804;
                    c(2,0) = -0.189651889660;
                    
                }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
    static void sdirk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        // A−stable SDIRK method
        switch (s) { // Implicit midpoint (1, 2)
            case 1:
                {
                    a(0,0) = 0.5;
                    b(0,0) = 1.0;
                    c(0,0) = 0.5;
                }
                break;
            case 2:
                {
                    // The only 2-stages A-stable SDIRK scheme of order 3 has been given by Crouzeix
                    // Crouzeix/LDD-A2(2, 3)
                    double gamma = 1.0/std::sqrt(3.0);
                    a(0,0) = 0.5 + 0.5*gamma;
                    a(1,0) = -1.0*gamma;
                    a(1,1) = 0.5 + 0.5*gamma;

                    b(0,0) = 0.5;
                    b(1,0) = 0.5;

                    c(0,0) = 0.5 + 0.5*gamma;
                    c(1,0) = 0.5 - 0.5*gamma;
                    
                }
                break;
            case 3:
                {
                    
                    // Crouzeix & Raviart (3, 4)
                    double gamma = (1.0/std::sqrt(3.0))*std::cos(M_PI/18.0)+0.5;
                    double delta = 1.0/(6.0*(2.0*gamma-1.0)*(2.0*gamma-1.0));
                    a(0,0) = gamma;
                    a(1,0) = (1.0/2.0)-gamma;
                    a(1,1) = gamma;
                    a(2,0) = 2.0*gamma;
                    a(2,1) = 1.0-4.0*gamma;
                    a(2,2) = gamma;
                    
                    b(0,0) = delta;
                    b(1,0) = 1.0-2.0*delta;
                    b(2,0) = delta;
                    
                    c(0,0) = gamma;
                    c(1,0) = 1.0/2.0;
                    c(2,0) = 1.0-gamma;
                    
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

#endif /* dirk_butcher_tableau_hpp */
