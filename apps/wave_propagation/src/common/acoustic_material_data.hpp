//
//  acoustic_material_data.hpp
//  acoustics
//
//  Created by Omar Durán on 4/14/20.
//

#pragma once
#ifndef acoustic_material_data_hpp
#define acoustic_material_data_hpp

#include <stdio.h>

template<typename T = double>
class acoustic_material_data {
        
    /// Fluid density
    T m_rho;
    
    /// Compressional P-wave velocity
    T m_vp;
    
public:
    
    /// Default constructor
    acoustic_material_data(T rho, T vp){
        m_rho = rho;
        m_vp = vp;
    }
    
    /// Copy constructor
    acoustic_material_data(const acoustic_material_data & other){
        m_rho       = other.m_rho;
        m_vp         = other.m_vp;
    }
    
    /// Assignement constructor
    const acoustic_material_data & operator=(const acoustic_material_data & other){
        
        // check for self-assignment
        if(&other == this){
            return *this;
        }
        
        m_rho       = other.m_rho;
        m_vp         = other.m_vp;
        return *this;
        
    }
    
    /// Desconstructor
    virtual ~acoustic_material_data(){
        
    }
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const{
        out << "\n density = " << m_rho;
        out << "\n p-wave velocity = " << m_vp;
    }
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const acoustic_material_data & material ){
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
    
};

#endif /* acoustic_material_data_hpp */
