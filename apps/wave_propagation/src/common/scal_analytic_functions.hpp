//
//  scal_analytic_functions.hpp
//  acoustics
//
//  Created by Omar Durán on 4/9/20.
//

#pragma once
#ifndef scal_analytic_functions_hpp
#define scal_analytic_functions_hpp

#define contrast 3.0
#define n_terms 200

class scal_analytic_functions {

public:
  
  /// Enumerate defining the function type
  enum EFunctionType { EFunctionNonPolynomial        = 0, EFunctionQuadraticInTime = 1,
                       EFunctionQuadraticInSpace     = 2, EFunctionQuadraticInSpaceTime = 3,
		       EFunctionInhomogeneousInSpace = 4, Fluid_pulse = 5};
    
    
  scal_analytic_functions(){
    m_function_type = EFunctionNonPolynomial;
  }
  
  ~scal_analytic_functions(){  
  }
  
  void set_function_type(EFunctionType function_type){
    m_function_type = function_type;
  }
  
  std::function<double(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_u(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return (1.0/(std::sqrt(2.0)*M_PI))*std::sin(std::sqrt(2.0)*M_PI*t) * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return (1 - pt.x())*pt.x()*(1 - pt.y())*pt.y()*std::sin(std::sqrt(2.0)*M_PI*t);
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return t*t*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpaceTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return t*t*(1 - pt.x())*pt.x()*(1 - pt.y())*pt.y();
                        };
                }
                break;
            case EFunctionInhomogeneousInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        
                        size_t n = n_terms;
                        double x,y,c1,c2;
                        c1 = contrast;
                        c2 = 1.0;
                        x = pt.x();
                        y = pt.y();

                        double p = 0.0;
                        double c = 1.0;
                        for (int k = 0; k <= n; k++) {

                            if (x < 0.5) {
                                double plvp, plvm;
                                plvp = (1.0/100.0)*exp(-(20.0*((k+x-c1*t)-0.2))*(20.0*((k+x-c1*t)-0.2)));
                                plvm = (1.0/100.0)*exp(-(20.0*((k-x-c1*t)-0.2))*(20.0*((k-x-c1*t)-0.2)));
                                p+=c*(plvp-plvm);
                            }else{
                                double pr;
                                pr = (1.0/100.0)*exp(-(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2))*(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2)));
                                p+=((2.0*c1)/(c2+c1))*c*pr;
                            }
                            c*=(c2-c1)/(c2+c1);
                        }
                        return p;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 0;
                    };
            }
                break;
        }
        
    }
    
    std::function<double(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_v(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return std::cos(std::sqrt(2.0)*M_PI*t) * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return std::sqrt(2.0)*M_PI*(1 - pt.x())*pt.x()*(1 - pt.y())*pt.y()*std::cos(std::sqrt(2.0)*M_PI*t);
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 2.0*t*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpaceTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 2.0*t*(1 - pt.x())*pt.x()*(1 - pt.y())*pt.y();
                        };
                }
                break;
            case EFunctionInhomogeneousInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        
                            size_t n = n_terms;
                            double x,y,c1,c2;
                            c1 = contrast;
                            c2 = 1.0;
                            x = pt.x();
                            y = pt.y();

                            double v = 0.0;
                            double c = 1.0;
                            for (int k = 0; k <= n; k++) {

                                if (x < 0.5) {
                                    double vlvp, vlvm;
                                    vlvp = (8.0*c1)* exp(-(20.0*((k+x-c1*t)-0.2))*(20.0*((k+x-c1*t)-0.2)))*((k+x-c1*t)-0.2);
                                    vlvm = (8.0*c1)* exp(-(20.0*((k-x-c1*t)-0.2))*(20.0*((k-x-c1*t)-0.2)))*((k-x-c1*t)-0.2);
                                    v+=c*(vlvp-vlvm);
                                }else{
                                    double vl;
                                    vl = (8.0*c1)*exp(-(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2))*(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2)))*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2);
                                    v+=((2.0*c1)/(c2+c1))*c*vl;
                                }
                                c*=(c2-c1)/(c2+c1);
                            }
                            return v;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 0;
                    };
            }
                break;
        }
        
    }
    
    std::function<double(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_a(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return -std::sqrt(2.0) * M_PI * std::sin(std::sqrt(2.0)*M_PI*t) * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return -2*M_PI*M_PI*(1 - pt.x())*pt.x()*(1 - pt.y())*pt.y()*std::sin(std::sqrt(2)*M_PI*t);
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 2.0*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpaceTime:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 2.0*(1 - pt.x())*pt.x()*(1 - pt.y())*pt.y();
                        };
                }
                break;
            case EFunctionInhomogeneousInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        
                            size_t n = n_terms;
                            double x,y,c1,c2;
                            c1 = contrast;
                            c2 = 1.0;
                            x = pt.x();
                            y = pt.y();

                            double a = 0.0;
                            double c = 1.0;
                            for (int k = 0; k <= n; k++) {

                                if (x < 0.5) {
                                    double alvp, alvm, xip, xim;
                                    xip=((k+x-c1*t)-0.2);
                                    alvp = (8.0*c1*c1)* exp(-(20.0*((k+x-c1*t)-0.2))*(20.0*((k+x-c1*t)-0.2)))*(800.0*xip*xip-1.0);
                                    
                                    xim=((k-x-c1*t)-0.2);
                                    alvm = (8.0*c1*c1)* exp(-(20.0*((k-x-c1*t)-0.2))*(20.0*((k-x-c1*t)-0.2)))*(800.0*xim*xim-1.0);
                                    a+=c*(alvp-alvm);
                                }else{
                                    double al, xi;
                                    xi= (((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2);
                                    al = (8.0*c1*c1)*exp(-(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2))*(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2)))*(800.0*xi*xi-1.0);
                                    a+=((2.0*c1)/(c2+c1))*c*al;
                                }
                                c*=(c2-c1)/(c2+c1);
                            }
                        return a;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 0;
                    };
            }
                break;
        }
        
    }
    
    std::function<double(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_f(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return 0;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            double x,y,f;
                            x = pt.x();
                            y = pt.y();
                            f = 2*(x - x*x + y - M_PI*M_PI*(-1 + x)*x*(-1 + y)*y - y*y)*std::sin(std::sqrt(2.0)*M_PI*t);
                            return f;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 2.0*(1.0 + M_PI*M_PI*t*t)*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                        };
                }
                break;
            case EFunctionQuadraticInSpaceTime:
                 {
                     return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            double x,y,f;
                            x = pt.x();
                            y = pt.y();
                            f = 2.0*((-1.0 + x)*x*(-1.0 + y)*y + t*t*(x - x*x + y - y*y));
                            return f;
                         };
                 }
                 break;
            case EFunctionInhomogeneousInSpace:
                {
                    return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                            return 0;
                        };
                }
                break;
            case Fluid_pulse:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {    
                        double x, y, xc, yc, tau, freq, Ricker, alpha, sigma, sigma2, sigma3, eps, f;

                        // Source term ponctuel
                        x   = pt.x(); y  = pt.y();
                        xc  = 0.5;    yc = 0.5;
                        
                        if (x==xc && y==yc) {    
                            
                            tau   = 0.4;
                            freq  = 5.0; 
                            alpha = - M_PI*M_PI * freq*freq;
                            
                            if ( t < 0.5*tau ) {
                                sigma  = alpha*(t-tau)*(t-tau);
                                Ricker = 2.0 * alpha * (1+2*sigma) * exp(sigma);
                            }
                            else {
                                Ricker = 0; 
                            }

                            sigma  = M_PI*freq*(t-tau);
                            sigma2 = sigma*sigma;
                            sigma3 = M_PI*freq;
                            
                            f = -4.0*sigma*sigma3*exp(-sigma2)-2.0*sigma*sigma3*Ricker;

                            return f;
                        }

                        else {
                            f = 0; 
                            return f;
                        }                  

                    return f;

                    };
                }
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                        return 0;
                    };
            }
                break;
        }
        
    }
    
    
    std::function<std::vector<double>(const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_q(double & t){
        
        switch (m_function_type) {
            case EFunctionNonPolynomial:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            std::vector<double> flux(2);
                            flux[0] = (std::sin(std::sqrt(2)*M_PI*t)*std::cos(M_PI*x)*std::sin(M_PI*y))/std::sqrt(2.0);
                            flux[1] = (std::sin(std::sqrt(2)*M_PI*t)*std::sin(M_PI*x)*std::cos(M_PI*y))/std::sqrt(2.0);
                            return flux;
                        };
                }
                break;
            case EFunctionQuadraticInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            std::vector<double> flux(2);
                            flux[0] = (1 - x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t) - x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            flux[1] = (1 - x)*x*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t) - (1 - x)*x*y*std::sin(std::sqrt(2.0)*M_PI*t);
                            return flux;
                        };
                }
                break;
            case EFunctionQuadraticInTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            std::vector<double> flux(2);
                            flux[0] = M_PI*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                            flux[1] = M_PI*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                            return flux;
                        };
                }
                break;
            case EFunctionQuadraticInSpaceTime:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                            double x,y;
                            x = pt.x();
                            y = pt.y();
                            std::vector<double> flux(2);
                            flux[0] = t*t*(1 - x)*(1 - y)*y - t*t*x*(1 - y)*y;
                            flux[1] = t*t*(1 - x)*x*(1 - y) - t*t*(1 - x)*x*y;
                            return flux;
                        };
                }
                break;
            case EFunctionInhomogeneousInSpace:
                {
                    return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                        
                            size_t n = n_terms;
                            double x,y,c1,c2;
                            c1 = contrast;
                            c2 = 1.0;
                            x = pt.x();
                            y = pt.y();

                            double qx = 0.0;
                            double c = 1.0;
                            for (int k = 0; k <= n; k++) {

                                if (x < 0.5) {
                                    double qxlvp, qxlvm;
                                    qxlvp = -(8.0)* exp(-(20.0*((k-x-c1*t)-0.2))*(20.0*((k-x-c1*t)-0.2)))*((k-x-c1*t)-0.2);
                                    qxlvm = -(8.0)* exp(-(20.0*((k+x-c1*t)-0.2))*(20.0*((k+x-c1*t)-0.2)))*((k+x-c1*t)-0.2);
                                    qx+=c1*c1*c*(qxlvp+qxlvm);
                                }else{
                                    double qxl;
                                    qxl = -(8.0*(c1/c2))*exp(-(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2))*(20.0*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2)))*(((c1/c2)*(x-0.5)+0.5+k-c1*t)-0.2);
                                    qx+=c2*c2*((2.0*c1)/(c2+c1))*c*qxl;
                                }
                                c*=(c2-c1)/(c2+c1);
                            }
                        
                            std::vector<double> flux(2);
                            flux[0] = qx;
                            flux[1] = 0.0;
                            return flux;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> std::vector<double> {
                        std::vector<double> f;
                        return f;
                    };
            }
                break;
        }
        
    }
    
    private:
    
    EFunctionType m_function_type;
  
};

#endif /* scal_analytic_functions_hpp */
