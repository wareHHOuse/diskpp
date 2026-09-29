/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <iostream>

namespace disk {

namespace mechanics {

class InitInfo {
  public:
    double m_time_rigi, m_time_solve, m_time_dyna;

    InitInfo() : m_time_rigi( 0. ), m_time_solve( 0. ), m_time_dyna( 0. ) {}
};

class AssemblyInfo {
  public:
    size_t m_linear_system_size;
    double m_time_assembly;
    double m_time_gradrec, m_time_statcond, m_time_stab, m_time_postpro;
    double m_time_elem, m_time_law, m_time_contact, m_time_dyna, m_time_assembler;
    double m_time_rigi, m_time_load, m_time_fint;

    AssemblyInfo()
        : m_linear_system_size( 0 ),
          m_time_assembly( 0.0 ),
          m_time_gradrec( 0.0 ),
          m_time_statcond( 0.0 ),
          m_time_stab( 0.0 ),
          m_time_elem( 0.0 ),
          m_time_law( 0.0 ),
          m_time_contact( 0.0 ),
          m_time_postpro( 0.0 ),
          m_time_dyna( 0.0 ),
          m_time_assembler( 0.0 ),
          m_time_rigi( 0. ),
          m_time_load( 0. ),
          m_time_fint( 0. ) {}

    AssemblyInfo &operator+=( const AssemblyInfo &other ) {
        m_linear_system_size = other.m_linear_system_size;
        m_time_assembly += other.m_time_assembly;
        m_time_gradrec += other.m_time_gradrec;
        m_time_statcond += other.m_time_statcond;
        m_time_stab += other.m_time_stab;
        m_time_elem += other.m_time_elem;
        m_time_law += other.m_time_law;
        m_time_contact += other.m_time_contact;
        m_time_dyna += other.m_time_dyna;
        m_time_postpro += other.m_time_postpro;
        m_time_assembler += other.m_time_assembler;
        m_time_rigi += other.m_time_rigi;
        m_time_fint += other.m_time_fint;
        m_time_load += other.m_time_load;
        return *this;
    }
};

class SolveInfo {
  public:
    size_t m_linear_system_size, m_nonzeros;
    double m_time_solve;

    SolveInfo() : m_linear_system_size( 0 ), m_nonzeros( 0 ), m_time_solve( 0.0 ) {}
    SolveInfo( const size_t linear_system_size, const size_t nonzeros, const double time_solve )
        : m_linear_system_size( linear_system_size ),
          m_nonzeros( nonzeros ),
          m_time_solve( time_solve ) {}
};

class NewtonSolverInfo {
  public:
    AssemblyInfo m_assembly_info;
    SolveInfo m_solve_info;
    double m_time_newton;
    size_t m_iter;

    NewtonSolverInfo() : m_assembly_info(), m_solve_info(), m_time_newton( 0.0 ), m_iter( 0 ) {}

    void updateAssemblyInfo( const AssemblyInfo &assembly_info ) {
        m_assembly_info += assembly_info;
    }

    void updateSolveInfo( const SolveInfo &solve_info ) {
        m_solve_info.m_linear_system_size = solve_info.m_linear_system_size;
        m_solve_info.m_nonzeros = solve_info.m_nonzeros;
        m_solve_info.m_time_solve += solve_info.m_time_solve;
    }

    void updateInitInfo( const InitInfo &init_info ) {
        m_assembly_info.m_time_dyna += init_info.m_time_dyna;
        m_assembly_info.m_time_rigi += init_info.m_time_rigi;
        m_solve_info.m_time_solve += init_info.m_time_solve;
    }

    void printInfo() const {
        std::cout << "** Time in this step " << m_time_newton << " sec" << std::endl;
        std::cout << "**** Assembly time: " << m_assembly_info.m_time_assembly << " sec"
                  << std::endl;
        std::cout << "****** Gradient reconstruction: " << m_assembly_info.m_time_gradrec << " sec"
                  << std::endl;
        std::cout << "****** Stabilisation: " << m_assembly_info.m_time_stab << " sec" << std::endl;
        std::cout << "****** Mechanical computation: " << m_assembly_info.m_time_elem << " sec"
                  << std::endl;
        std::cout << "       *** Behavior computation: " << m_assembly_info.m_time_law << " sec"
                  << std::endl;
        std::cout << "       *** Contact computation: " << m_assembly_info.m_time_contact << " sec"
                  << std::endl;
        std::cout << "       *** Dynamic computation: " << m_assembly_info.m_time_dyna << " sec"
                  << std::endl;
        std::cout << "       *** Rigidity computation: " << m_assembly_info.m_time_rigi << " sec"
                  << std::endl;
        std::cout << "       *** Force int computation: " << m_assembly_info.m_time_fint << " sec"
                  << std::endl;
        std::cout << "       *** Load computation: " << m_assembly_info.m_time_load << " sec"
                  << std::endl;
        std::cout << "****** Static condensation: " << m_assembly_info.m_time_statcond << " sec"
                  << std::endl;
        std::cout << "****** Postprocess time: " << m_assembly_info.m_time_postpro << " sec"
                  << std::endl;
        std::cout << "****** Assembling time: " << m_assembly_info.m_time_assembler << " sec"
                  << std::endl;
        std::cout << "**** Solver time: " << m_solve_info.m_time_solve << " sec" << std::endl;
    }

    std::map< std::string, double > getValues() const {
        std::map< std::string, double > vals;

        vals["nb_iter"] = m_iter;
        vals["total_time"] = m_time_newton;
        vals["assembly_time"] = m_assembly_info.m_time_assembly;
        vals["gradrec_time"] = m_assembly_info.m_time_gradrec;
        vals["stab_time"] = m_assembly_info.m_time_stab;
        vals["law_time"] = m_assembly_info.m_time_law;
        vals["contact_time"] = m_assembly_info.m_time_contact;
        vals["dyna_time"] = m_assembly_info.m_time_dyna;
        vals["rigi_time"] = m_assembly_info.m_time_rigi;
        vals["fint_time"] = m_assembly_info.m_time_fint;
        vals["load_time"] = m_assembly_info.m_time_load;
        vals["stat_cond"] = m_assembly_info.m_time_statcond;
        vals["post_time"] = m_assembly_info.m_time_postpro;
        vals["assembler_time"] = m_assembly_info.m_time_assembler;
        vals["solver_time"] = m_solve_info.m_time_solve;

        return vals;
    }
};

class SolverInfo {
  public:
    NewtonSolverInfo m_newton_info;
    double m_time_solver;
    size_t m_iter, m_time_step;

    SolverInfo() : m_newton_info(), m_time_solver( 0.0 ), m_iter( 0 ), m_time_step( 0 ) {}

    void updateInfo( const NewtonSolverInfo &newton_solver_info ) {
        m_iter += newton_solver_info.m_iter;
        m_newton_info.updateAssemblyInfo( newton_solver_info.m_assembly_info );
        m_newton_info.updateSolveInfo( newton_solver_info.m_solve_info );
        m_newton_info.m_time_newton += newton_solver_info.m_time_newton;
    }

    void printInfo() const {
        std::cout << " " << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << "Summaring: " << std::endl;
        std::cout << "Total Newton's iterations: " << m_iter << " in " << m_time_step
                  << " load increments" << std::endl;
        std::cout << "Total time to solve the problem: " << m_time_solver << " sec" << std::endl;
        std::cout << "**** Assembly time: " << m_newton_info.m_assembly_info.m_time_assembly
                  << " sec" << std::endl;
        std::cout << "****** Gradient reconstruction: "
                  << m_newton_info.m_assembly_info.m_time_gradrec << " sec" << std::endl;
        std::cout << "****** Stabilisation: " << m_newton_info.m_assembly_info.m_time_stab << " sec"
                  << std::endl;
        std::cout << "****** Elementary computation: " << m_newton_info.m_assembly_info.m_time_elem
                  << " sec" << std::endl;
        std::cout << "       *** Behavior computation: " << m_newton_info.m_assembly_info.m_time_law
                  << " sec" << std::endl;
        std::cout << "       *** Contact computation: "
                  << m_newton_info.m_assembly_info.m_time_contact << " sec" << std::endl;
        std::cout << "       *** Dynamic computation: " << m_newton_info.m_assembly_info.m_time_dyna
                  << " sec" << std::endl;
        std::cout << "       *** Rigidity computation: "
                  << m_newton_info.m_assembly_info.m_time_rigi << " sec" << std::endl;
        std::cout << "       *** Force int computation: "
                  << m_newton_info.m_assembly_info.m_time_fint << " sec" << std::endl;
        std::cout << "       *** Load computation: " << m_newton_info.m_assembly_info.m_time_load
                  << " sec" << std::endl;
        std::cout << "****** Static condensation: " << m_newton_info.m_assembly_info.m_time_statcond
                  << " sec" << std::endl;
        std::cout << "****** Postprocess time: " << m_newton_info.m_assembly_info.m_time_postpro
                  << " sec" << std::endl;
        std::cout << "****** Assembling time: " << m_newton_info.m_assembly_info.m_time_assembler
                  << " sec" << std::endl;
        std::cout << "**** Solver time: " << m_newton_info.m_solve_info.m_time_solve << " sec"
                  << std::endl;
        std::cout << "------------------------------------------------------- " << std::endl;
        std::cout << " " << std::endl;
    }
};

} // namespace mechanics
} // namespace disk