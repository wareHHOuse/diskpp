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
 * cite it.
 */

 #pragma once


/*

Fichier pour gérer les lois de comportements

-faire une boite qui prends en entrée eps --> behaviorbox --> stress
-penser aux variables internes pour sauvegarder partie elastique
- utilise t on des tenseur pour le module tangent
*/


template<typename scalar_type>
class Material_parameters
{
    scalar_type m_Young;
    scalar_type m_Poisson;

    scalar_type m_mu;
    scalar_type m_lambda;

    void
    converttoLame()
    {
      m_lambda = m_Young * m_Poisson / (( 1.0 + m_Poisson)*(1.0 -2.0 * m_Poisson));
      m_mu = m_Young  / (2.0*( 1.0 + m_Poisson));
    }

    void
    converttoHooke()
    {
      m_Young = m_mu * (3.0 * m_lambda + 2.0 * m_mu) / ( m_lambda + m_mu);
      m_Poisson = m_lambda  / (2.0*( m_lambda + m_mu));

    }

public:
     Material_parameters()
     : m_Young(0.0), m_Poisson(0.0), m_mu(0.0), m_lambda(0.0) {}

     Material_parameters(const scalar_type mu, const scalar_type lambda)
     : m_mu(mu), m_lambda(lambda)
     {
         converttoHooke();
     }

    void
    setYoung(const scalar_type young)
    {
        m_Young = young;
        converttoLame();
    }

    void
    setPoisson(const scalar_type poisson)
    {
        m_Poisson = poisson;
        converttoLame();
    }

    void
    setMu(const scalar_type mu)
    {
        m_mu = mu;
        converttoHooke();
    }

    void
    setLambda(const scalar_type lambda)
    {
        m_lambda = lambda;
        converttoHooke();
    }

    void
    setHooke(const scalar_type young, const scalar_type poisson)
    {
        m_Young = young;
        m_Poisson = poisson;
        converttoLame();
    }

    void
    setLame(const scalar_type mu, const scalar_type lambda)
    {
        m_mu = mu;
        m_lambda = lambda;
        converttoHooke();
    }

    scalar_type
    giveYoung() const {return m_Young;}

    scalar_type
    givePoisson() const {return m_Poisson;}

    scalar_type
    giveMu() const {return m_mu;}

    scalar_type
    giveLambda() const {return m_lambda;}

};
