/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#pragma once

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

enum StabilizationType : size_t
{
   HDG  = 0,
   PIKF = 1,
   HHO  = 2,
   NO   = 3
};

template<typename T>
class ParamRun
{
 public:
   size_t m_face_degree; // face degree
   size_t m_cell_degree; // cell_degree
   size_t m_grad_degree; // grad degree

   size_t m_leray_param; // leray_paramater (p)

   std::vector<std::pair<T, size_t>> m_time_step; // number of time time_step
   size_t                            m_sublevel;  // number of sublevel if there are problems
   size_t                            m_iter_max;  // maximun nexton iteration
   T                                 m_epsilon;   // stop criteria

   bool m_verbose; // some printing

   bool m_precomputation; // to compute the gradient before (it's memory consuption)

   size_t m_stab_type; // type of stabilization
   T      m_beta;      // stabilization parameter
   bool   m_stab;      // stabilization yes or no
   bool   m_pred;      // prediction (with p =2) yes or no

   size_t       m_n_time_save; // number of saving
   std::list<T> m_time_save;   // list of time where we save result;

   ParamRun() :
     m_face_degree(1), m_cell_degree(1), m_grad_degree(1), m_leray_param(2), m_sublevel(1),
     m_iter_max(20), m_epsilon(T(1E-6)), m_verbose(false), m_precomputation(false), m_stab(true),
     m_beta(1), m_stab_type(HHO), m_n_time_save(0), m_pred(false)
   {
      m_time_step.push_back(std::make_pair(1.0, 10));
   }

   void
   infos()
   {
      std::cout << "Parameters:" << std::endl;
      std::cout << " - Face degree: " << m_face_degree << std::endl;
      std::cout << " - Cell degree: " << m_cell_degree << std::endl;
      std::cout << " - Grad degree: " << m_grad_degree << std::endl;
      std::cout << " - LerayParam: " << m_leray_param << std::endl;
      std::cout << " - Prediction (LerayParam =2): " << m_pred << std::endl;
      std::cout << " - Stabilization ?: " << m_stab << std::endl;
      std::cout << " - Type: " << m_stab_type << std::endl;
      std::cout << " - Beta: " << m_beta << std::endl;
      std::cout << " - Verbose: " << m_verbose << std::endl;
      std::cout << " - Sublevel: " << m_sublevel << std::endl;
      std::cout << " - IterMax: " << m_iter_max << std::endl;
      std::cout << " - Epsilon: " << m_epsilon << std::endl;
      std::cout << " - Precomputation: " << m_precomputation << std::endl;
   }

   bool
   readParameters(const std::string& filename)
   {
      std::ifstream ifs(filename);
      std::string   keyword;
      size_t        line(0);

      if (!ifs.is_open()) {
         std::cout << "Error opening " << filename << std::endl;
         return false;
      }

      ifs >> keyword;
      line++;
      if (keyword != "BeginParameters") {
         std::cout << "Expected keyword \"BeginParameters\" line: " << line << std::endl;
         return false;
      }

      ifs >> keyword;
      line++;
      while (keyword != "EndParameters") {
         if (keyword == "FaceDegree") {
            ifs >> m_face_degree;
            line++;
         } else if (keyword == "CellDegree") {
            ifs >> m_cell_degree;
            line++;
         } else if (keyword == "GradDegree") {
            ifs >> m_grad_degree;
            line++;
         } else if (keyword == "LerayParam") {
            ifs >> m_leray_param;
            line++;
         } else if (keyword == "Prediction") {
            std::string logical;
            ifs >> logical;
            line++;
            if (logical == "true")
               m_pred = true;
            else {
               m_pred = false;
            }
         } else if (keyword == "Sublevel") {
            ifs >> m_sublevel;
            line++;
         } else if (keyword == "TimeStep") {
            size_t n_time_step(0);
            ifs >> n_time_step;
            line++;

            m_time_step.clear();
            m_time_step.reserve(n_time_step);
            for (size_t i = 0; i < n_time_step; i++) {
               T      time(0.0);
               size_t time_step(0);
               ifs >> time >> time_step;
               m_time_step.push_back(std::make_pair(time, time_step));
               line++;
            }
         } else if (keyword == "TimeSave") {
            ifs >> m_n_time_save;
            line++;

            m_time_save.clear();
            for (size_t i = 0; i < m_n_time_save; i++) {
               T time(0.0);
               ifs >> time;
               m_time_save.push_back(time);
               line++;
            }
         } else if (keyword == "Stabilization") {
            std::string logical;
            ifs >> logical;
            line++;
            if (logical == "true")
               m_stab = true;
            else {
               m_stab      = false;
               m_stab_type = NO;
            }
         } else if (keyword == "StabType") {
            std::string type;
            ifs >> type;
            line++;

            if (type == "HDG")
               m_stab_type = HDG;
            else if (type == "PIKF")
               m_stab_type = PIKF;
            else if (type == "HHO")
               m_stab_type = HHO;
            else if (type == "NO")
               m_stab_type = NO;
         } else if (keyword == "Beta") {
            ifs >> m_beta;
            line++;
         } else if (keyword == "Verbose") {
            std::string logical;
            ifs >> logical;
            line++;
            if (logical == "true")
               m_verbose = true;
            else
               m_verbose = false;
         } else if (keyword == "IterMax") {
            ifs >> m_iter_max;
            line++;
         } else if (keyword == "Epsilon") {
            ifs >> m_epsilon;
            line++;
         } else if (keyword == "Precomputation") {
            std::string logical;
            ifs >> logical;
            line++;
            if (logical == "true")
               m_precomputation = true;
            else
               m_precomputation = false;
         } else {
            std::cout << "Error parsing Parameters file:" << keyword << " line: " << line
                      << std::endl;
            return false;
         }

         ifs >> keyword;
         line++;
      }

      ifs.close();
      return true;
   }
};
