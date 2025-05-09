//
//  elastic_material_data.hpp
//  acoustics
//
//  Created by Omar Durán on 4/14/20.
//

#ifndef elastic_material_data_hpp
#define elastic_material_data_hpp

#include <stdio.h>

template<typename T = double>
class elastic_material_data {
        
    /// Fluid density
    T m_rho;
    
    /// Compressional P-wave velocity
    T m_vp;
    
    /// Shear S-wave velocity
    T m_vs;
    
public:
    
    /// Default constructor
    elastic_material_data(T rho, T vp, T vs){
        m_rho = rho;
        m_vp = vp;
        m_vs = vs;
    }
    
    /// Copy constructor
    elastic_material_data(const elastic_material_data & other){
        m_rho       = other.m_rho;
        m_vp         = other.m_vp;
        m_vs         = other.m_vs;
    }
    
    /// Assignement constructor
    const elastic_material_data & operator=(const elastic_material_data & other){
        
        // check for self-assignment
        if(&other == this){
            return *this;
        }
        
        m_rho       = other.m_rho;
        m_vp         = other.m_vp;
        m_vs         = other.m_vs;
        return *this;
        
    }
    
    /// Desconstructor
    virtual ~elastic_material_data(){
        
    }
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const{
        out << "\n density = " << m_rho;
        out << "\n p-wave velocity = " << m_vp;
        out << "\n s-wave velocity = " << m_vs;
    }
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const elastic_material_data & material ){
        material.Print(out);
        return out;
    }
    
    void Set_rho(T rho)
    {
        m_rho = rho;
    }
    
    T rho()
    {
        return m_rho;
    }
    
    
    void Set_vp(T vp)
    {
        m_vp = vp;
    }
    
    T vp()
    {
        return m_vp;
    }
    
    void Set_vs(T vs)
    {
        m_vs = vs;
    }
    
    T vs()
    {
        return m_vs;
    }
    
};

#endif /* elastic_material_data_hpp */
