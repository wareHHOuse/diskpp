//
//  lsrk_butcher_tableau.hpp
//  acoustics
//
//  Created by Omar Durán on 4/25/20.
//

#ifndef lsrk_butcher_tableau_hpp
#define lsrk_butcher_tableau_hpp

#include <Eigen/Dense>

class lsrk_butcher_tableau
{
    public:
    
    static void lsrk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        
        // Nodal Discontinuous Galerkin Methods
        // Fourth-order 2n-storage Runge- Kutta schemes
        switch (s) {
            case 1:
                {
                    a(0,0) = 0.0;
                    b(0,0) = 1.0;
                    c(0,0) = 1.0;
                }
                break;
            case 2:
                {
                    a(0,0) = 1.0/3.0;
                    a(1,0) = 3.0/4.0;
                    a(1,1) = 1.0/4.0;
                    
                    b(0,0) = 3.0/4.0;
                    b(1,0) = 1.0/4.0;
                    
                    c(0,0) = 1.0/3.0;
                    c(1,0) = 1.0;
                    
                }
                break;
            case 3:
                {
                    
                    a(0,0) = 1.0;
                    a(1,0) = -1.0/12.0;
                    a(1,1) = 5.0/12.0;
                    a(2,0) = 0.0;
                    a(2,1) = 3.0/4.0;
                    a(2,2) = 1.0/4.0;
                    
                    b(0,0) = 0.0;
                    b(1,0) = 3.0/4.0;
                    b(2,0) = 3.0/4.0;
                    
                    c(0,0) = 1.0;
                    c(1,0) = 1.0/3.0;
                    c(2,0) = 1.0;
                    
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

#endif /* lsrk_butcher_tableau_hpp */
