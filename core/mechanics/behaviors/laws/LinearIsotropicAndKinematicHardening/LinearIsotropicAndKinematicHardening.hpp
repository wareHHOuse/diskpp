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

#include <vector>

#include "common/eigen.hpp"
#include "mechanics/behaviors/laws/LinearIsotropicAndKinematicHardening/LinearIsotropicAndKinematicHardening_cell.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "mechanics/behaviors/maths_utils.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law forLinearIsotropicAndKinematicHardening model in small deformations

template<typename MeshType>
class LinearIsotropicAndKinematicHardening
{
  private:
    typedef MeshType                        mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell        cell_type;

    size_t                                                            m_nb_qp;
    std::vector<LinearIsotropicAndKinematicHardening_cell<mesh_type>> m_list_cell_qp;
    LinearIsotropicAndKinematicHardening_Data<scalar_type>            m_data;

  public:
    LinearIsotropicAndKinematicHardening() : m_nb_qp(0){};

    LinearIsotropicAndKinematicHardening(const mesh_type& msh, const int degree)
    {
        m_nb_qp = 0;
        m_list_cell_qp.clear();
        m_list_cell_qp.reserve(msh.cells_size());

        for (auto& cl : msh)
        {
            LinearIsotropicAndKinematicHardening_cell<mesh_type> cell_qp(msh, cl, degree);

            m_list_cell_qp.push_back(cell_qp);
            m_nb_qp += cell_qp.getNumberOfQP();
        }
    }

    void
    addMaterialData(const scalar_type& lambda,
                    const scalar_type& mu,
                    const scalar_type& H,
                    const scalar_type& K,
                    const scalar_type& sigma_y0)
    {
        m_data = LinearIsotropicAndKinematicHardening_Data<scalar_type>(lambda, mu, H, K, sigma_y0);
    }

    LinearIsotropicAndKinematicHardening_Data<scalar_type>
    getMaterialData() const
    {
        return m_data;
    }

    int
    getNumberOfQP() const
    {
        return m_nb_qp;
    }

    void
    update()
    {
        for (auto& qp_cell : m_list_cell_qp)
        {
            qp_cell.update();
        }
    }

    LinearIsotropicAndKinematicHardening_cell<mesh_type>&
    getCellQPs(const int cell_id)
    {
        return m_list_cell_qp.at(cell_id);
    }

    LinearIsotropicAndKinematicHardening_cell<mesh_type>
    getCellIVs(const int cell_id) const
    {
        return m_list_cell_qp.at(cell_id);
    }
};
}
