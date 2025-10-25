
//  Created by Omar Durán
//  Contributor: Romain Mottier

#pragma once
#ifndef scal_vec_analytic_functions_hpp
#define scal_vec_analytic_functions_hpp

class scal_vec_analytic_functions {
    
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { 
      EFunctionNonPolynomial            = 0,  
      EFunctionQuadraticInTime          = 1,
      EFunctionQuadraticInSpace         = 2, 
      reproduction_acoustic             = 3,
      reproduction_elastic              = 4,
      EFunctionNonPolynomial_paper      = 5, 
      EFunctionCubicInTimeAcoustic      = 6,
      EFunctionQuarticInTimeAcoustic    = 7,
      EFunctionQuadraticInSpaceAcoustic = 8,
      EFunctionPlaneWaveAcoustic = 9};
      
      scal_vec_analytic_functions() {
          m_function_type = EFunctionNonPolynomial;
      }
      
      ~scal_vec_analytic_functions() {  
      }
          
      void set_function_type(EFunctionType function_type) {
          m_function_type = function_type;
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_u(double & t) {
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = x*x*std::sin(M_PI*y)*std::cos((M_PI*x)/2.)*std::cos(std::sqrt(2.0)*M_PI*t);
                      uy = x*x*std::sin(M_PI*y)*std::cos((M_PI*x)/2.)*std::cos(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = x*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      uy = x*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = (1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      uy = (1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case reproduction_elastic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = -t*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                      uy =  t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ux,uy, w, theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      ux = x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      uy = x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> u;
                      return u;
                  };
              }
              break;
          }   
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_v(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = -(std::sqrt(2.0)*M_PI*x*x*std::cos((M_PI*x)/2.)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y));
                      vy = -(std::sqrt(2.0)*M_PI*x*x*std::cos((M_PI*x)/2.)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y));
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = 2.0*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      vy = 2.0*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = std::sqrt(2.0)*M_PI*(1 + x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                      vy = std::sqrt(2.0)*M_PI*(1 + x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case reproduction_elastic:{
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = -2*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                      vy =  2*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper:{
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,vx,vy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      vx = -(theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos((M_PI*w*x)/2.));
                      vy = -(theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos((M_PI*w*x)/2.));
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_vector<double, 2> v;
                      return v;
                  };
              }
              break;    
          }
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_a(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*M_PI*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y);
                      ay = -2*M_PI*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = 2.0*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      ay = 2.0*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*M_PI*M_PI*(1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      ay = -2*M_PI*M_PI*(1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case reproduction_elastic:{
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*std::cos(M_PI*y)*std::sin(M_PI*x);
                      ay = 2*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      ax = -theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      ay = -theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> f;
                      return f;
                  };
              }
              break;   
          }   
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_f(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -(std::cos(std::sqrt(2.0)*M_PI*t)*(-4*M_PI*x*std::sin((M_PI*x)/2.)*(M_PI*x*std::cos(M_PI*y) + 6*std::sin(M_PI*y)) + std::cos((M_PI*x)/2.)*(16*M_PI*x*std::cos(M_PI*y) + (24 + M_PI*M_PI*x*x)*std::sin(M_PI*y))))/4.0;
                      fy = (std::cos(std::sqrt(2.0)*M_PI*t)*(4*M_PI*x*std::sin((M_PI*x)/2.)*(M_PI*x*std::cos(M_PI*y) + 2*std::sin(M_PI*y)) + std::cos((M_PI*x)/2.)*(-16*M_PI*x*std::cos(M_PI*y) + (-8 + 5*M_PI*M_PI*x*x)*std::sin(M_PI*y))))/4.;
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -2 * ( M_PI*t*t*std::cos(M_PI*x)*( M_PI*x*std::cos(M_PI*y) + 3*std::sin(M_PI*y) ) +
                      std::sin(M_PI*x)*(M_PI*t*t*std::cos(M_PI*y) - (1+2*M_PI*M_PI*t*t)*x*std::sin(M_PI*y)));
                      fy = (1 + M_PI*M_PI*t*t)*x*std::cos(M_PI*(x - y)) - (1+3*M_PI*M_PI*t*t)*x*std::cos(M_PI*(x+y))
                      - 2*M_PI*t*t*std::sin(M_PI*(x+y));
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = 2 * ( x * (-2+(-2+x)*x) -3*y -x*(5+x*(-6+M_PI*M_PI*(1+x)))*y
                           + (3+x*(9+M_PI*M_PI*x*(1+x)))*y*y)*std::sin(std::sqrt(2.0)*M_PI*t);
                      fy = 2 * ( x*x*(6 + M_PI*M_PI*(-1 + y))*y + (-1 + y)*y + x*x*x*(3 + M_PI*M_PI*(-1 + y)*y)
                           + x*(-2 + y + 3*y*y))*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case reproduction_elastic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*y)*std::sin(M_PI*x);
                      fy = 2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      fx = - theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0) 
                           + M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*x/2)*std::cos(theta*t)*std::cos(M_PI*w*y) 
                           + 7*M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*cos(M_PI*w*x/2)/4.0  
                           + 6*M_PI*w*x*std::sin(M_PI*w*x/2)*std::sin(M_PI*w*y)*std::cos(theta*t) 
                           - 4*M_PI*w*x*std::cos(theta*t)*std::cos(M_PI*w*x/2)*std::cos(M_PI*w*y) 
                           - 6*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2);
                           
                      fy = - theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0) 
                           + M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*x/2)*std::cos(theta*t)*std::cos(M_PI*w*y) 
                           + 13*M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*cos(M_PI*w*x/2)/4.0 
                           + 2*M_PI*w*x*std::sin(M_PI*w*x/2)*std::sin(M_PI*w*y)*std::cos(theta*t) 
                           - 4*M_PI*w*x*std::cos(theta*t)*std::cos(M_PI*w*x/2)*std::cos(M_PI*w*y) 
                           - 2*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              default:
              {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> f;
                      return f;
                  };
              }
              break;
          }    
      }
      
      std::function<disk::static_matrix<double,2,2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_sigma(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) + 4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) -
                      M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y) +
                      (4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) - M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.0;
                           
                      sigma(0,1) = M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) + 2*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) -
                      (M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.;
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = 3*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) +
                      (4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) - M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.;
                      return sigma;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = t*t*(M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + 3*(M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      
                      sigma(0,1) = t*t*(M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + (M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = t*t*(3*M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + (M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      return sigma;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) =   2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1+x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + (2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0
                      + (2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0;
                      
                      sigma(0,1) =   x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 2*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) =   2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + (2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0
                      + (2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0;
                      return sigma;
                  };
              }
              break;
              
              case reproduction_elastic:
              {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = -2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                      sigma(1,1) = 2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                      return sigma;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = - 3.0*M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0               
                      + M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y)   
                      + 6.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      
                      sigma(0,1) = - M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0 
                      + M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y)                     
                      + 2.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = - M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0 
                      + 3*M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y) 
                      + 2.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      return sigma;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      disk::static_matrix<double,2,2> sigma(2,2);
                      return sigma;
                  };
              }
              break;
          }
          
      }
      
      // ACOUSTIC
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_u(double & t) {
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return x*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return (1 - x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return y*(1-x*x)*(1-y)*std::sin(M_PI*t);
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      return x*x*std::sin(theta*t)*std::sin(w*M_PI*x)*std::sin(w*M_PI*y);
                  };
              }
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return -(1.0/std::sqrt(2)) * std::sin(theta);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_v(double & t){
                
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return std::sqrt(2.0)*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
                    
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return std::sqrt(2.0)*M_PI*(1 - x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return M_PI*y*(1-x*x)*(1-y)*std::cos(M_PI*t);
                  };
              }
              break;

              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 2.0*t*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      return theta*x*x*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)*std::cos(theta*t);
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 4.0*t*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 3.0*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return std::cos(theta);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
          
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_a(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -2*M_PI*M_PI*x*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
                    
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -2*M_PI*M_PI*(1 - x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -M_PI*M_PI*y*(1-x*x)*(1-y)*std::sin(M_PI*t);
                  };
              }
              break;

              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2.0*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
                    
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      return -theta*theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y);
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 12*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 6*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return std::sqrt(2) * std::sin(theta);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }    
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_f(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = -2*std::sin(std::sqrt(2.0)*M_PI*t)*(2*M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*(-(M_PI*t*t*std::cos(M_PI*x)) + (1 + M_PI*M_PI*t*t)*x*std::sin(M_PI*x))*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*((-1 + y)*y - 3*x*(-1 + y)*y + x*x*x*(-1 - M_PI*M_PI*(-1 + y)*y)
                      + x*x*(1 + M_PI*M_PI*(-1 + y)*y))*std::sin(std::sqrt(2.0)*M_PI*t);
                      return f;
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = - M_PI*M_PI*y*(1-x*x)*(1-y)*std::sin(M_PI*t)
                          + 2*y*(1-y)*std::sin(M_PI*t)
                          + 2*(1-x*x)*std::sin(M_PI*t);
                      return f;
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      return 2.0*(1.0 + M_PI*M_PI*t*t)*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      f = - theta*theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y) 
                      + 2*M_PI*M_PI*theta*w*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta 
                      - 4*M_PI*theta*w*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos(M_PI*w*x)/theta 
                      - 2*theta*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta;
                      return f;
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*M_PI*M_PI*t*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y) + 12*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*M_PI*M_PI*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y) + 6*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 0.0;
                      return f;
                  };
              }
              break;

              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_q(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*x*x*std::cos(M_PI*x)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y) + 2*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*x*x*std::cos(M_PI*y)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = M_PI*t*t*x*std::cos(M_PI*x)*std::sin(M_PI*y) + t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*x*std::cos(M_PI*y)*std::sin(M_PI*x);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = 2*(1 - x)*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                           - x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      qy = (1 - x)*x*x*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                           - (1 - x)*x*x*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = -2*x*y*(1-y)*std::sin(M_PI*t);
                      qy = - y*(1-x*x)*std::sin(M_PI*t) 
                           + (1-x*x)*(1-y)*std::sin(M_PI*t);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = 1.0; 
                      theta = 10*M_PI; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      qx = M_PI*theta*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos(M_PI*w*x)/theta + 2*theta*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta;
                      qy = M_PI*theta*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::cos(M_PI*w*y)/theta;
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x, y, qx, qy, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      qx = (1.0/std::sqrt(2.0)) * std::sin(theta);
                      qy = (1.0/std::sqrt(2.0)) * std::sin(theta);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;

              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2>{
                      disk::static_vector<double, 2> f{0,0};
                      return f;
                  };
              }
              break;
          }    
      }
      
      private:
      
      EFunctionType m_function_type;
      
  };
  
  #endif /* scal_vec_analytic_functions_hpp */
  