//
//  vec_analytic_functions.hpp
//  acoustics
//
//  Created by Omar Durán on 4/14/20.
//

#pragma once
#ifndef vec_analytic_functions_hpp
#define vec_analytic_functions_hpp

#include <stdio.h>

class vec_analytic_functions
{
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { EFunctionNonPolynomial = 0, EFunctionQuadraticInTime = 1, EFunctionQuadraticInSpace = 2};
    
    
    vec_analytic_functions(){
        m_function_type = EFunctionNonPolynomial;
    }
    
    ~vec_analytic_functions(){
        
    }
    
    void set_function_type(EFunctionType function_type){
        m_function_type = function_type;
    }

    std::function<static_vector<double, 2>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_u(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,ux,uy;
                            x = pt.x();
                            y = pt.y();
                            ux = -std::cos(M_PI*y)*std::sin(std::sqrt(2)*M_PI*t)*std::sin(M_PI*x);
                            uy =  std::cos(M_PI*x)*std::sin(std::sqrt(2)*M_PI*t)*std::sin(M_PI*y);
                            static_vector<double, 2> u{ux,uy};
                            return u;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,ux,uy;
                            x = pt.x();
                            y = pt.y();
                            ux = (1.0 - x)*x*(1.0 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            uy = (1.0 - x)*x*(1.0 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            static_vector<double, 2> u{ux,uy};
                            return u;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,ux,uy;
                            x = pt.x();
                            y = pt.y();
                            ux = -t*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                            uy =  t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                            static_vector<double, 2> u{ux,uy};
                            return u;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                        static_vector<double, 2> u;
                        return u;
                    };
            }
                break;
        }
        
    }
    
    std::function<static_vector<double, 2>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_v(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,vx,vy;
                            x = pt.x();
                            y = pt.y();
                            vx = -(std::sqrt(2)*M_PI*std::cos(std::sqrt(2)*M_PI*t)*std::cos(M_PI*y)*std::sin(M_PI*x));
                            vy = std::sqrt(2)*M_PI*std::cos(std::sqrt(2)*M_PI*t)*std::cos(M_PI*x)*std::sin(M_PI*y);
                            static_vector<double, 2> v{vx,vy};
                            return v;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {

                            double x,y,vx,vy;
                            x = pt.x();
                            y = pt.y();
                            vx = std::sqrt(2.0)*M_PI*(1.0 - x)*x*(1.0 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                            vy = std::sqrt(2.0)*M_PI*(1.0 - x)*x*(1.0 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                            static_vector<double, 2> v{vx,vy};
                            return v;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,vx,vy;
                            x = pt.x();
                            y = pt.y();
                            vx = -2*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                            vy =  2*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                            static_vector<double, 2> v{vx,vy};
                            return v;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                        double x,y;
                        x = pt.x();
                        y = pt.y();
                        static_vector<double, 2> v;
                        return v;
                    };
            }
                break;
        }
        
    }
    
    std::function<static_vector<double, 2>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_a(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,ax,ay;
                             x = pt.x();
                             y = pt.y();
                             ax = 2*M_PI*M_PI*std::cos(M_PI*y)*std::sin(std::sqrt(2)*M_PI*t)*std::sin(M_PI*x);
                             ay = -2*M_PI*M_PI*std::cos(M_PI*x)*std::sin(std::sqrt(2)*M_PI*t)*std::sin(M_PI*y);
                             static_vector<double, 2> a{ax,ay};
                            return a;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,ax,ay;
                            x = pt.x();
                            y = pt.y();
                            ax = -2.0*M_PI*M_PI*(1.0 - x)*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            ay = -2.0*M_PI*M_PI*(1.0 - x)*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            static_vector<double, 2> a{ax,ay};
                            return a;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                        double x,y,ax,ay;
                        x = pt.x();
                        y = pt.y();
                        ax = -2*std::cos(M_PI*y)*std::sin(M_PI*x);
                        ay = 2*std::cos(M_PI*x)*std::sin(M_PI*y);
                        static_vector<double, 2> a{ax,ay};
                        return a;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                        static_vector<double, 2> f;
                        return f;
                    };
            }
                break;
        }
        
    }
    
    std::function<static_vector<double, 2>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_f(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            static_vector<double, 2> f{0,0};
                            return f;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                            double x,y,fx,fy;
                            x = pt.x();
                            y = pt.y();
                            fx = -2.0*(1.0 + (-3.0 + x)*x - 5.0*y + (4.0 - M_PI*M_PI*(-1.0 + x))*x*y + (3.0 + M_PI*M_PI*(-1.0 + x)*x)*y*y)*std::sin(std::sqrt(2.0)*M_PI*t);
                            fy = -2*(1 + (-3 + y)*y + x*(-5 + (4 - M_PI*M_PI*(-1 + y))*y) + x*x*(3 + M_PI*M_PI*(-1 + y)*y))*std::sin(std::sqrt(2.0)*M_PI*t);
                            static_vector<double, 2> f{fx,fy};
                            return f;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                        double x,y,fx,fy;
                        x = pt.x();
                        y = pt.y();
                        fx = -2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*y)*std::sin(M_PI*x);
                        fy = 2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*x)*std::sin(M_PI*y);
                        static_vector<double, 2> f{fx,fy};
                        return f;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_vector<double, 2> {
                    static_vector<double, 2> f;
                    return f;
                    };
            }
                break;
        }
        
    }
    
    
    std::function<static_matrix<double,2,2>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_sigma(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_matrix<double,2,2> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            static_matrix<double,2,2> sigma = static_matrix<double,2,2>::Zero(2,2);
                            sigma(0,0) = -2*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(std::sqrt(2)*M_PI*t);
                            sigma(1,1) = 2*M_PI*std::cos(M_PI*x)*std::cos(M_PI*y)*std::sin(std::sqrt(2)*M_PI*t);
                            return sigma;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_matrix<double,2,2> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            static_matrix<double,2,2> sigma = static_matrix<double,2,2>::Zero(2,2);
                            sigma(0,0) = 2*(1 - x)*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t) - 2*x*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t) +
                            (2*(1 - x)*x*(1 - y)*std::sin(std::sqrt(2)*M_PI*t) - 2*(1 - x)*x*y*std::sin(std::sqrt(2)*M_PI*t))/2. +
                            (2*(1 - x)*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t) - 2*x*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t))/2.;
                            sigma(0,1) = (1 - x)*x*(1 - y)*std::sin(std::sqrt(2)*M_PI*t) - (1 - x)*x*y*std::sin(std::sqrt(2)*M_PI*t) + (1 - x)*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t) -
                            x*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t);
                            sigma(1,0) = sigma(0,1);
                            sigma(1,1) = 2*(1 - x)*x*(1 - y)*std::sin(std::sqrt(2)*M_PI*t) - 2*(1 - x)*x*y*std::sin(std::sqrt(2)*M_PI*t) +
                            (2*(1 - x)*x*(1 - y)*std::sin(std::sqrt(2)*M_PI*t) - 2*(1 - x)*x*y*std::sin(std::sqrt(2)*M_PI*t))/2. +
                            (2*(1 - x)*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t) - 2*x*(1 - y)*y*std::sin(std::sqrt(2)*M_PI*t))/2.;
                             return sigma;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_matrix<double,2,2> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            static_matrix<double,2,2> sigma = static_matrix<double,2,2>::Zero(2,2);
                            sigma(0,0) = -2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                            sigma(1,1) = 2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                            return sigma;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> static_matrix<double,2,2> {
                        static_matrix<double,2,2> sigma(2,2);
                         return sigma;
                    };
            }
                break;
        }
        
    }
    
    private:
    
    EFunctionType m_function_type;
  
};

#endif /* vec_analytic_functions_hpp */
