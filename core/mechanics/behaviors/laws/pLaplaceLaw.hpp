/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018                     nicolas.pignet@enpc.fr
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

#include "common/eigen.hpp"

namespace disk
{

/* Material: LaplaceLaw
 * Stress :  S(G) = norm(G)^{p-2} * G
 * Module :  A(G) = d lambda / dG \times G  + lambda(G) * I_4
 */

template<typename scalar_type>
class pLaplaceLaw
{
    size_t m_p;

    template<int DIM>
    static_matrix<scalar_type, DIM, DIM>
    compute_A(const static_vector<scalar_type, DIM>& G) const
    {
        if (m_p == 2)
            return static_matrix<scalar_type, DIM, DIM>::Identity();

        const scalar_type norm_G = G.norm();

        if (norm_G <= 10E-16 && m_p < 4)
            return static_matrix<scalar_type, DIM, DIM>::Identity();

        const int pm2 = m_p - 2;
        const int pm4 = m_p - 4;

        const scalar_type norm_Gm2 = std::pow(norm_G, pm2);

        scalar_type norm_Gm4 = 0;
        if (pm4 == 0)
        {
            norm_Gm4 = 1.0;
        }
        else if (pm4 == -1)
        {
            norm_Gm4 = 1.0 / norm_G;
        }
        else
        {
            norm_Gm4 = std::pow(norm_G, pm4);
        }

        return norm_Gm2 * static_matrix<scalar_type, DIM, DIM>::Identity() + (pm2)*norm_Gm4 * G * G.transpose();
    }

  public:
    pLaplaceLaw() : m_p(2) {}

    pLaplaceLaw(const size_t p) : m_p(p)
    {
        if (m_p < 2)
        {
            throw std::invalid_argument("pLaplaceLaw: p have to be >= 2");
        }
    }

    template<int DIM>
    static_vector<scalar_type, DIM>
    compute_stress(const static_vector<scalar_type, DIM>& G) const
    {
        if (m_p == 2)
            return G;
        else
            return std::pow(G.norm(), m_p - 2) * G;
    }

    scalar_type
    compute_stress(const scalar_type& G) const
    {
        if (m_p == 2)
            return G;
        else
            return std::pow(std::abs(G), m_p - 2) * G;
    }

    scalar_type
    compute_tangent_moduli(const scalar_type& G) const
    {
        if (m_p == 2)
            return 1;
        else
        {
            return (m_p - 1) * std::pow(std::abs(G), m_p - 2);
        }
    }

    template<int DIM>
    static_matrix<scalar_type, DIM, DIM>
    compute_tangent_moduli(const static_vector<scalar_type, DIM>& G) const
    {
        if (m_p == 2)
            return static_matrix<scalar_type, DIM, DIM>::Identity();
        else
            return compute_A(G);
    }

    template<int DIM>
    std::pair<static_vector<scalar_type, DIM>, static_matrix<scalar_type, DIM, DIM>>
    compute_whole(const static_vector<scalar_type, DIM>& G) const
    {
        const auto P = compute_stress(G);
        const auto A = compute_tangent_moduli(G);

        return std::make_pair(P, A);
    }

    std::pair<scalar_type, scalar_type>
    compute_whole(const scalar_type& G) const
    {
        const auto P = compute_stress(G);
        const auto A = compute_tangent_moduli(G);

        return std::make_pair(P, A);
    }
};

} // end LerayLions