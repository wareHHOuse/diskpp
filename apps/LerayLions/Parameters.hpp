/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite
*/

#pragma once

#include <fstream>
#include <iostream>
#include <string>

enum StabilizationType : size_t
{
   L2 = 0,
   PIKF = 1,
   HHO = 2,
   NOTHING = 3
};

template< typename T>
class ParamRun
{
public:
   size_t  m_face_degree;     //face degree
   size_t  m_cell_degree;    //cell_degree = face_degree + l
   size_t  m_grad_degree;     // grad degree

   std::vector<std::pair<T, size_t>>  m_time_step;   //number of time time_step
   size_t  m_sublevel;      //number od sublevel if there are problems

   bool    m_stab;          //stabilization yes or no
   bool    m_adapt_coeff;   //adapts automatically the stabilisation coefficient
   bool    m_adapt_stab;    //use the adpatative stabilization
   bool    m_verbose;       //some printing
   bool    m_compute_energy; // to compute intere energy
   bool    m_precomputation; // to compute the gradient before (it's memory consuption)

   size_t  m_stab_init;      //type of stabilization
   size_t  m_stab_obj;      //type of stabilization
   size_t  m_iter_max;        //maximun nexton iteration
   T       m_epsilon;         //stop criteria

   T       m_beta_max;
   T       m_beta_init;      // minium of stabilization constant
   T       m_beta_obj;      // maximum of stabilization constant

   size_t  m_n_time_save;  // number of saving
   std::list<T>   m_time_save; //list of time where we save result;


   ParamRun() : m_face_degree(1), m_cell_degree(1), m_grad_degree(1),
                m_sublevel(1), m_stab(true),
                m_verbose(false), m_adapt_coeff(false), m_adapt_stab(false),
                m_compute_energy(false), m_precomputation(false),
                m_iter_max(10), m_epsilon(T(1E-6)),
                m_beta_init(T(1.0)), m_beta_obj(T(1.0)), m_beta_max(T(1)),
                m_stab_init(L2), m_stab_obj(HHO),
                m_n_time_save(0)
                {
                   m_time_step.push_back(std::make_pair(1.0, 10));
                }

   void infos()
   {
      std::cout << "Parameters:" << std::endl;
      std::cout << " - Face degree: "  << m_face_degree << std::endl;
      std::cout << " - Cell degree: "  << m_cell_degree << std::endl;
      std::cout << " - Grad degree: "  << m_grad_degree << std::endl;
      std::cout << " - Sublevel: "  << m_sublevel << std::endl;
      std::cout << " - Stabilization ?: "  << m_stab << std::endl;
      std::cout << " - Adapt coefficient: "  << m_adapt_coeff << std::endl;
      std::cout << "    - Beta init: "  << m_beta_init << std::endl;
      std::cout << "    - Beta obj: "  << m_beta_obj << std::endl;
      std::cout << " - Adapt stabilization: "  << m_adapt_stab << std::endl;
      std::cout << "    - Stab init : "  << m_stab_init << std::endl;
      std::cout << "    - Stab obj: "  << m_stab_obj << std::endl;
      std::cout << " - Verbose: "  << m_verbose << std::endl;
      std::cout << " - Compute Energy: "  << m_compute_energy << std::endl;
      std::cout << " - IterMax: "  << m_iter_max << std::endl;
      std::cout << " - Epsilon: "  << m_epsilon << std::endl;
      std::cout << " - Precomputation: "  << m_precomputation << std::endl;
   }

   bool readParameters(const std::string& filename)
   {
      std::ifstream   ifs(filename);
      std::string     keyword;
      size_t line(0);

      if (!ifs.is_open())
      {
         std::cout << "Error opening " << filename << std::endl;
         return false;
      }

      ifs >> keyword;
      line++;
      if ( keyword != "BeginParameters" )
      {
         std::cout << "Expected keyword \"BeginParameters\" line: " << line << std::endl;
         return false;
      }


      ifs >> keyword;
      line++;
      while( keyword != "EndParameters")
      {
         if ( keyword == "FaceDegree" )
         {
            ifs >> m_face_degree;
            line++;
         }
         else if ( keyword == "CellDegree" )
         {
            ifs >> m_cell_degree;
            line++;
         }
         else if ( keyword == "GradDegree" )
         {
            ifs >> m_grad_degree;
            line++;
         }
         else if ( keyword == "Sublevel" )
         {
            ifs >> m_sublevel;
            line++;
         }
         else if ( keyword == "TimeStep" )
         {
            size_t n_time_step(0);
            ifs >> n_time_step;
            line++;

            m_time_step.clear();
            m_time_step.reserve(n_time_step);
            for (size_t i = 0; i < n_time_step; i++) {
               T  time(0.0);
               size_t   time_step(0);
               ifs >> time >> time_step;
               m_time_step.push_back(std::make_pair(time, time_step));
               line++;
            }
         }
         else if ( keyword == "TimeSave" )
         {
            ifs >> m_n_time_save;
            line++;

            m_time_save.clear();
            for (size_t i = 0; i < m_n_time_save; i++) {
               T  time(0.0);
               ifs >> time;
               m_time_save.push_back(time);
               line++;
            }
         }
         else if ( keyword == "Stabilization" )
         {
            std::string logical;
            ifs >> logical;
            line++;
            if(logical == "true")
               m_stab = true;
            else{
               m_stab = false;
               m_stab_init = NOTHING;
               m_stab_obj = NOTHING;
            }
         }
         else if ( keyword == "StabInit" )
         {
            std::string type;
            ifs >> type;
            line++;

            if(type == "L2")
               m_stab_init = L2;
            else if(type == "PIKF")
               m_stab_init = PIKF;
            else if(type == "HHO")
               m_stab_init = HHO;
            else if(type == "NOTHING")
               m_stab_init = NOTHING;
         }
         else if ( keyword == "StabObj" )
         {
            std::string type;
            ifs >> type;
            line++;

            if(type == "L2")
               m_stab_obj = L2;
            else if(type == "PIKF")
               m_stab_obj = PIKF;
            else if(type == "HHO")
               m_stab_obj = HHO;
            else if(type == "NOTHING")
               m_stab_obj = NOTHING;
         }
         else if ( keyword == "AdaptativeCoefficient" )
         {
            std::string logical;
            ifs >> logical;
            line++;

            if(logical == "true")
               m_adapt_coeff = true;
            else
               m_adapt_coeff = false;
         }
         else if ( keyword == "AdaptativeStabilization" )
         {
            std::string logical;
            ifs >> logical;
            line++;

            if(logical == "true")
               m_adapt_stab = true;
            else
               m_adapt_stab = false;
         }
         else if ( keyword == "BetaInit" )
         {
            ifs >> m_beta_init;
            line++;
         }
         else if ( keyword == "BetaMax" )
         {
            ifs >> m_beta_max;
            line++;
         }
         else if ( keyword == "BetaObj" )
         {
            ifs >> m_beta_obj;
            line++;
         }
         else if ( keyword == "Verbose" )
         {
            std::string logical;
            ifs >> logical;
            line++;
            if(logical == "true")
               m_verbose = true;
            else
               m_verbose = false;
         }
         else if ( keyword == "ComputeEnergy" )
         {
            std::string logical;
            ifs >> logical;
            line++;

            if(logical == "true")
               m_compute_energy = true;
            else
               m_compute_energy = false;
         }
         else if ( keyword == "IterMax" )
         {
            ifs >> m_iter_max;
            line++;
         }
         else if ( keyword == "Epsilon" )
         {
            ifs >> m_epsilon;
            line++;
         }
         else if ( keyword == "Precomputation" )
         {
            std::string logical;
            ifs >> logical;
            line++;
            if(logical == "true")
               m_precomputation = true;
            else
               m_precomputation = false;
         }
         else
         {
            std::cout << "Error parsing Parameters file:" << keyword << " line: "  << line << std::endl;
            return false;
         }

         ifs >> keyword;
         line++;
      }

      ifs.close();
      return true;
   }
};
