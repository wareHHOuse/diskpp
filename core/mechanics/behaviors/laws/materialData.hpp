/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
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

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

template<typename T>
class curve_point
{
  private:
    T m_p;
    T m_Rp;

  public:
    curve_point() : m_p(T(0)), m_Rp(T(0)) {}

    curve_point(const T p, const T Rp) : m_p(p), m_Rp(Rp) {}

    T
    getP() const
    {
        return m_p;
    }

    T
    getRp() const
    {
        return m_Rp;
    }
};

template<typename scalar_type>
class MaterialData
{
  private:
    scalar_type                           m_lambda;
    scalar_type                           m_mu;
    scalar_type                           m_H;
    scalar_type                           m_K;
    scalar_type                           m_sigma_y0;
    size_t                                m_type;
    std::vector<curve_point<scalar_type>> m_Rp_curve;

  public:
    MaterialData() :
      m_lambda(1.0), m_mu(1.0), m_H(0), m_K(0), m_sigma_y0(std::numeric_limits<scalar_type>::max()), m_type(0)
    {
    }

    MaterialData(const scalar_type& lambda,
                 const scalar_type& mu,
                 const scalar_type& H,
                 const scalar_type& K,
                 const scalar_type& sigma_y0) :
      m_lambda(lambda),
      m_mu(mu), m_H(H), m_K(K), m_sigma_y0(sigma_y0), m_type(0)
    {
    }

    MaterialData(const scalar_type& lambda, const scalar_type& mu) :
      m_lambda(lambda), m_mu(mu), m_H(0), m_K(0), m_sigma_y0(std::numeric_limits<scalar_type>::max()), m_type(0)
    {
    }

    void
    setMu(const scalar_type mu)
    {
        m_mu = mu;
    }

    void
    setMu(const scalar_type E, const scalar_type nu)
    {
        m_mu = E / (2.0 * (1.0 + nu));
    }

    void
    setLambda(const scalar_type lambda)
    {
        m_lambda = lambda;
    }

    void
    setLambda(const scalar_type E, const scalar_type nu)
    {
        m_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    }

    void
    setH(const scalar_type H)
    {
        m_H = H;
    }

    void
    setH(const scalar_type E, const scalar_type ET, const scalar_type K)
    {
        m_H = E * ET / (E - ET) - 1.5 * K;
    }

    void
    setK(const scalar_type K)
    {
        m_K = K;
    }

    void
    setSigma_y0(const scalar_type sigma_y0)
    {
        m_sigma_y0 = sigma_y0;
    }

    void
    setType(const size_t type)
    {
        m_type = type;
    }

    void
    setRpCurve(const std::vector<curve_point<scalar_type>> RpCurve)
    {
        m_Rp_curve = RpCurve;
    }

    void
    addCurvePoint(const curve_point<scalar_type> point)
    {
        m_Rp_curve.push_back(point);
    }

    void
    addCurvePoint(const scalar_type p, const scalar_type Rp)
    {
        m_Rp_curve.push_back(curve_point<scalar_type>(p, Rp));
    }

    scalar_type
    getE() const
    {
        return m_mu * (3 * m_lambda + 2 * m_mu) / (m_lambda + m_mu);
    }

    scalar_type
    getNu() const
    {
        return m_lambda / (2 * (m_lambda + m_mu));
    }

    scalar_type
    getET() const
    {
        const scalar_type E = getE();
        return E * (m_H + 1.5 * m_K) / (m_H + 1.5 * m_K + E);
    }

    scalar_type
    getLambda() const
    {
        return m_lambda;
    }

    scalar_type
    getMu() const
    {
        return m_mu;
    }

    scalar_type
    getH() const
    {
        return m_H;
    }

    scalar_type
    getK() const
    {
        return m_K;
    }

    scalar_type
    getSigma_y0() const
    {
        return m_sigma_y0;
    }

    size_t
    getType() const
    {
        return m_type;
    }

    std::vector<curve_point<scalar_type>>
    getRpCurve() const
    {
        return m_Rp_curve;
    }

    void
    checkRpCurve()
    {
        if (m_Rp_curve.size() > 0)
        {
            std::sort(m_Rp_curve.begin(), m_Rp_curve.end(), [](const auto& lhs, const auto& rhs) {
                return lhs.getP() < rhs.getP();
            });

            for (size_t i = 0; i < m_Rp_curve.size() - 1; i++)
            {
                if (std::abs(m_Rp_curve[i].getP() - m_Rp_curve[i + 1].getP()) <
                    std::numeric_limits<scalar_type>::epsilon())
                {
                    throw std::invalid_argument("RpCurve: You have two values with the same p");
                }
            }
        }
    }

    void
    print() const
    {
        std::cout << "Material parameters: " << std::endl;
        std::cout << "* E: " << getE() << std::endl;
        std::cout << "* Nu: " << getNu() << std::endl;
        std::cout << "* ET: " << getET() << std::endl;
        std::cout << "* H: " << getH() << std::endl;
        std::cout << "* K: " << getK() << std::endl;
        std::cout << "* Sy0: " << getSigma_y0() << std::endl;
        std::cout << "* Lambda: " << getLambda() << std::endl;
        std::cout << "* Mu: " << getMu() << std::endl;
        std::cout << "* Traction Curve:" << std::endl;
        std::cout << "(p, R(p))" << std::endl;
        for (auto& pt : m_Rp_curve)
            std::cout << "( " << pt.getP() << ", " << pt.getRp() << " )" << std::endl;
    }
};
}