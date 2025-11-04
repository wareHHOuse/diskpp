/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018, 2024               nicolas.pignet@enpc.fr
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

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

enum StabilizationType : int
{
    HDG = 0,
    HHO = 1,
    HHO_SYM = 4,
    NO = 2,
    DG = 3
};

enum FrictionType : int
{
    NO_FRICTION = 0,
    TRESCA      = 1,
    COULOMB     = 2,
};

enum DynamicType : int
{
    STATIC = 0,
    HHT = 1,
    NEWMARK = 2,
};

template<typename T>
class NewtonSolverParameter
{
  public:
    int m_face_degree; // face degree
    int m_cell_degree; // cell_degree
    int m_grad_degree; // grad degree

    std::vector<std::pair<T, int>> m_time_step;         // number of time time_step
    bool                           m_has_user_end_time; // final time is given
    T                              m_user_end_time;     // final time of the simulation
    int                            m_sublevel;          // number of sublevel if there are problems
    int                            m_iter_max;          // maximun nexton iteration
    T                              m_epsilon;           // stop criteria

    bool m_verbose; // some printing

    bool m_precomputation; // to compute the gradient before (it's memory consuption)

    int  m_stab_type;  // type of stabilization
    T    m_beta;       // stabilization parameter
    bool m_stab;       // stabilization yes or no
    bool m_adapt_stab; // adaptative stabilization

    DynamicType m_dyna_type;              // type of dyna
    bool m_dynamic;                       // dynamic or static simulation
    std::map<std::string, T> m_dyna_para; // list of parameters

    int          m_n_time_save; // number of saving
    std::list<T> m_time_save;   // list of time where we save result;

    T   m_theta;     // theta-parameter for contact
    T   m_gamma_0;   // parameter for Nitsche
    T   m_threshold; // threshol for Tesca friction
    int m_frot_type; // Friction type ?

    NewtonSolverParameter() : m_face_degree(1), m_cell_degree(1), m_grad_degree(1), m_sublevel(5), m_iter_max(20), m_epsilon(T(1E-6)),
                              m_verbose(false), m_precomputation(false), m_stab(true), m_beta(1), m_stab_type(HHO), m_n_time_save(0),
                              m_user_end_time(1.0), m_has_user_end_time(false), m_adapt_stab(false), m_dynamic(false), m_theta(1), m_gamma_0(1),
                              m_threshold(0), m_frot_type(NO_FRICTION), m_dyna_type(DynamicType::STATIC)
    {
        m_time_step.push_back(std::make_pair(m_user_end_time, 1));
    }

    void
    infos()
    {
        std::cout << "Newton Solver's parameters:" << std::endl;
        std::cout << " - Face degree: " << m_face_degree << std::endl;
        std::cout << " - Cell degree: " << m_cell_degree << std::endl;
        std::cout << " - Grad degree: " << m_grad_degree << std::endl;
        std::cout << " - Stabilization ?: " << m_stab << std::endl;
        std::cout << " - AdaptativeStabilization ?: " << m_adapt_stab << std::endl;
        std::cout << " - Type: " << m_stab_type << std::endl;
        std::cout << " - Beta: " << m_beta << std::endl;
        std::cout << " - Verbose: " << m_verbose << std::endl;
        std::cout << " - Sublevel: " << m_sublevel << std::endl;
        std::cout << " - IterMax: " << m_iter_max << std::endl;
        std::cout << " - Epsilon: " << m_epsilon << std::endl;
        std::cout << " - Precomputation: " << m_precomputation << std::endl;
        std::cout << " - Dynamic: " << m_dynamic << std::endl;
        std::cout << " - Friction ?: " << m_frot_type << std::endl;
        std::cout << " - Threshold: " << m_threshold << std::endl;
        std::cout << " - Gamma_0: " << m_gamma_0 << std::endl;
        std::cout << " - Theta: " << m_theta << std::endl;
    }

    bool
    readParameters(const std::string& filename)
    {
        std::ifstream ifs(filename);
        std::string   keyword;
        int           line(0);

        if (!ifs.is_open())
        {
            std::cout << "Error opening " << filename << std::endl;
            return false;
        }

        ifs >> keyword;
        line++;
        if (keyword != "BeginParameters")
        {
            std::cout << "Expected keyword \"BeginParameters\" line: " << line << std::endl;
            return false;
        }

        ifs >> keyword;
        line++;
        while (keyword != "EndParameters")
        {
            if (keyword == "FaceDegree")
            {
                ifs >> m_face_degree;
                line++;
            }
            else if (keyword == "CellDegree")
            {
                ifs >> m_cell_degree;
                line++;
            }
            else if (keyword == "GradDegree")
            {
                ifs >> m_grad_degree;
                line++;
            }
            else if (keyword == "Sublevel")
            {
                ifs >> m_sublevel;
                line++;
            }
            else if (keyword == "TimeStep")
            {
                int n_time_step(0);
                ifs >> n_time_step;
                line++;

                m_time_step.clear();
                m_time_step.reserve(n_time_step);
                for (int i = 0; i < n_time_step; i++)
                {
                    T   time(0.0);
                    int time_step(0);
                    ifs >> time >> time_step;
                    m_time_step.push_back(std::make_pair(time, time_step));
                    line++;
                }
            }
            else if (keyword == "FinalTime")
            {
                ifs >> m_user_end_time;
                line++;

                m_has_user_end_time = true;
            }
            else if (keyword == "TimeSave")
            {
                ifs >> m_n_time_save;
                line++;

                m_time_save.clear();
                for (int i = 0; i < m_n_time_save; i++)
                {
                    T time(0.0);
                    ifs >> time;
                    m_time_save.push_back(time);
                    line++;
                }
            }
            else if (keyword == "Stabilization")
            {
                std::string logical;
                ifs >> logical;
                line++;
                if (logical == "true" || logical == "True")
                    m_stab = true;
                else
                {
                    m_stab      = false;
                    m_stab_type = NO;
                }
            }
            else if (keyword == "AdaptativeStabilization")
            {
                std::string logical;
                ifs >> logical;
                line++;
                m_adapt_stab = false;
                if (logical == "true" || logical == "True")
                    m_adapt_stab = true;
            }
            else if (keyword == "StabType")
            {
                std::string type;
                ifs >> type;
                line++;
                if (type == "HDG")
                    m_stab_type = HDG;
                else if (type == "HHO")
                    m_stab_type = HHO;
                else if (type == "HHO_SYM")
                    m_stab_type = HHO_SYM;
                else if (type == "DG")
                    m_stab_type = DG;
                else if (type == "NO")
                    m_stab_type = NO;
            }
            else if (keyword == "Beta")
            {
                ifs >> m_beta;
                line++;
            }
            else if (keyword == "Verbose")
            {
                std::string logical;
                ifs >> logical;
                line++;
                if (logical == "true" || logical == "True")
                    m_verbose = true;
                else
                    m_verbose = false;
            }
            else if (keyword == "IterMax")
            {
                ifs >> m_iter_max;
                line++;
            }
            else if (keyword == "Epsilon")
            {
                ifs >> m_epsilon;
                line++;
            }
            else if (keyword == "Precomputation")
            {
                std::string logical;
                ifs >> logical;
                line++;
                if (logical == "true" || logical == "True")
                    m_precomputation = true;
                else
                    m_precomputation = false;
            }
            else if (keyword == "Theta")
            {
                ifs >> m_theta;
                line++;
            }
            else if (keyword == "Gamma0")
            {
                ifs >> m_gamma_0;
                line++;
            }
            else if (keyword == "Friction")
            {
                std::string type;
                ifs >> type;
                line++;
                if (type == "NO")
                    m_frot_type = NO_FRICTION;
                else if (type == "TRESCA")
                    m_frot_type = TRESCA;
                else if (type == "COULOMB")
                    m_frot_type = COULOMB;
            }
            else if (keyword == "Threshold")
            {
                ifs >> m_threshold;
            }
            else if (keyword == "Dynamic")
            {
                std::string logical;
                ifs >> logical;
                line++;
                if (logical == "true" || logical == "True")
                    m_dynamic = true;
                else
                    m_dynamic = false;
            }
            else
            {
                std::cout << "Error parsing Parameters file:" << keyword << " line: " << line << std::endl;
                return false;
            }

            ifs >> keyword;
            line++;
        }

        ifs.close();
        return true;
    }

    void
    setFaceDegree(const int face_degree)
    {
        m_face_degree = face_degree;
    }

    int
    getFaceDegree() const
    {
        return m_face_degree;
    }

    void
    setCellDegree(const int cell_degree)
    {
        m_cell_degree = cell_degree;
    }

    int
    getCellDegree() const
    {
        return m_cell_degree;
    }

    void
    setGradDegree(const int grad_degree)
    {
        m_grad_degree = grad_degree;
    }

    int
    getGradDegree() const
    {
        return m_face_degree;
    }

    void
    setStabilizationParameter(const T stab_para)
    {
        m_beta = stab_para;
    }

    T
    getStabilizationParameter() const
    {
        return m_beta;
    }

    void
    setVerbose(const bool verbose)
    {
        m_verbose = verbose;
    }

    bool
    getVerbose() const
    {
        return m_verbose;
    }

    void
    setPrecomputation(const bool precomp)
    {
        m_precomputation = precomp;
    }

    bool
    getPrecomputation() const
    {
        return m_precomputation;
    }

    bool
    isUnsteady() const
    {
        return m_dynamic;
    }

    DynamicType getUnsteadyScheme() const
    {
        return m_dyna_type;
    }

    void
    setUnsteadyScheme(const DynamicType &scheme)
    {
        m_dyna_type = scheme;
    }

    void
    isUnsteady(const bool dyna)
    {
        m_dynamic = dyna;
    }

    auto
    getUnsteadyParameters() const
    {
        return m_dyna_para;
    }

    void
    setUnsteadyParameters(const std::map<std::string, T> dyna_para)
    {
        m_dyna_para = dyna_para;
    }

    void
    setTimeStep(const T end_time, const int n_time_step)
    {
        m_time_step.clear();
        m_time_step.push_back(std::make_pair(end_time, n_time_step));
    }
};
