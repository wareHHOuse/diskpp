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
#include "mechanics/behaviors/laws/LinearLaw/LinearLaw_cell.hpp"
#include "mechanics/behaviors/maths_tensor.hpp"
#include "mechanics/behaviors/maths_utils.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

namespace disk
{

// Law for LinearLaw model (to test the finite deformations)

template<typename MeshType>
class LinearLaw
{
  public:
    typedef MeshType                            mesh_type;
    typedef typename mesh_type::coordinate_type scalar_type;
    typedef typename mesh_type::cell            cell_type;
    typedef LinearLaw_Data<scalar_type>         data_type;
    typedef LinearLaw_cell<mesh_type>           law_cell_type;

  private:
    size_t                     m_nb_qp;
    std::vector<law_cell_type> m_list_cell_qp;
    data_type                  m_data;

  public:
    LinearLaw() : m_nb_qp(0){};

    LinearLaw(const mesh_type& msh, const size_t degree)
    {
        m_nb_qp = 0;
        m_list_cell_qp.clear();
        m_list_cell_qp.reserve(msh.cells_size());

        for (auto& cl : msh)
        {
            law_cell_type cell_qp(msh, cl, degree);

            m_list_cell_qp.push_back(cell_qp);
            m_nb_qp += cell_qp.getNumberOfQP();
        }
    }

    void
    addMaterialData(const scalar_type& lambda,
                    const scalar_type& mu)
    {
        m_data = data_type(lambda, mu);
    }

    data_type
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

    law_cell_type&
    getCellQPs(const int cell_id)
    {
        return m_list_cell_qp.at(cell_id);
    }

    law_cell_type
    getCellIVs(const int cell_id) const
    {
        return m_list_cell_qp.at(cell_id);
    }
};
}
