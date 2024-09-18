/*
 *       /\        Matteo Cicuttin (C) 2016-2021
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2021                     nicolas.pignet@enpc.fr
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
#include <list>
#include <vector>

namespace disk
{

namespace mechanics
{

/**
 * @brief Represent a time step on the time interval [start_time, end_time]
 *
 * @tparam T scalar type
 */
template<typename T>
class StabCoef
{
  private:
    T m_coeff;

  public:
    StabCoef() : m_coeff(T(1)) {}

    /**
     * @brief Construct a new StabCoef object
     *
     * @param coeff value of the stabilization coefficient
     */
    StabCoef(const T coeff) : m_coeff(coeff) {}

    /**
     * @brief Return the value of the stabiliztion coefficient
     *
     * @return T  value of the stabiliztion coefficient
     */
    T
    getValue(void) const
    {
        return m_coeff;
    }

    void
    setValue(const T& value)
    {
        m_coeff = value;
    }
};

/**
 * @brief This class StabCoeffManager allows to manage the different time step for the simulation
 *
 * @tparam T scalar type
 */
template<typename T>
class StabCoeffManager
{
  private:
    std::vector<StabCoef<T>> m_stab_coeff;
    std::vector<StabCoef<T>> m_stab_coeff_new;

  public:
    template<typename Mesh>
    StabCoeffManager(const Mesh& mesh, const T value)
    {
        m_stab_coeff.clear();
        m_stab_coeff.reserve(mesh.cells_size());
        m_stab_coeff_new.clear();
        m_stab_coeff_new.reserve(mesh.cells_size());

        for (auto& cl : mesh)
        {
            m_stab_coeff.push_back(value);
            m_stab_coeff_new.push_back(value);
        }
    }

    template<typename Mesh>
    T
    getValue(const Mesh& mesh, const typename Mesh::cell& cl) const
    {
        return m_stab_coeff[mesh.lookup(cl)].getValue();
    }

    template<typename Mesh>
    void
    setValue(const Mesh& mesh, const typename Mesh::cell& cl, const T& value)
    {
        m_stab_coeff[mesh.lookup(cl)].setValue(value);
    }

    template<typename Mesh>
    void
    setValueNext(const Mesh& mesh, const typename Mesh::cell& cl, const T& value)
    {
        m_stab_coeff_new[mesh.lookup(cl)].setValue(value);
    }

    void
    update()
    {
        m_stab_coeff = m_stab_coeff_new;
    }
};
}
}