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

#include "common/eigen.hpp"
#include "core/mechanics/behaviors/laws/behaviorlaws.hpp"
#include "core/mechanics/behaviors/logarithmic_strain/logarithmic_tools.hpp"
#include "core/mechanics/behaviors/maths_tensor.hpp"
#include "core/mechanics/behaviors/maths_utils.hpp"
#include "core/mechanics/behaviors/tensor_conversion.hpp"
#include "core/mechanics/deformation_tensors.hpp"
#include "core/mechanics/stress_tensors.hpp"
#include "mesh/point.hpp"

namespace disk
{

namespace mechanics
{

// Routine for Logarithmic Stain

/* For details see the paper:
 *   Anisotropic additive plasticity in the logaritmic strain: modular
 *   kinematic formulation and implementation based on incremental minimization
 *   principles for standard materials
 *   C. Miehe, N. Apel, M. Lambrecht
 *   Comput. Methods Appl. Mech. Engrg. (2002)
 */

template<typename LawTypeQp>
class LogarithmicStrain_qp
{
  public:
    typedef LawTypeQp                             law_hpp_qp_type;
    typedef typename law_hpp_qp_type::data_type   data_type;
    typedef typename law_hpp_qp_type::scalar_type scalar_type;

    const static size_t DIM = LawTypeQp::dimension;

    typedef static_matrix<scalar_type, DIM, DIM> static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>     static_matrix_type3D;
    typedef static_tensor<scalar_type, DIM>      tensor_type;

  private:
    // law hpp at qp
    law_hpp_qp_type m_law_hpp_qp;

    static_tensor<scalar_type, 3> Pn; // Projector

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        return computeContractedProduct(Pn, this->compute_stress3D_T(data));
    }

    static_matrix_type3D
    compute_stress3D_T(const data_type& data) const
    {
        return m_law_hpp_qp.compute_stress3D(data);
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        std::cout << "F" << std::endl;
        std::cout << F_curr << std::endl;

        const scalar_type J = F_curr.determinant();
        if (J <= 0.01)
        {
            const std::string mess = "J= " + std::to_string(J) + " <= 0";
            throw std::invalid_argument(mess);
        }

        const auto C    = convertFtoCauchyGreenRight(F_curr);
        const auto ev_C = compute_eigenvalues(C);
        const auto Elog = compute_Elog(ev_C.first, ev_C.second);

        std::cout << "C" << std::endl;
        std::cout << C << std::endl;
        std::cout << ev_C.first << std::endl;
        std::cout << ev_C.second << std::endl;
        std::cout << "Elog" << std::endl;
        std::cout << Elog << std::endl;

        const auto behavior3D_hpp = m_law_hpp_qp.compute_whole3D(Elog, data, tangentmodulus);
        const auto projector      = compute_projector(F_curr, behavior3D_hpp.first, ev_C.first, ev_C.second);

        Pn             = projector.first;
        const auto PK1 = this->compute_stress3D(data);

        std::cout << "T" << std::endl;
        std::cout << behavior3D_hpp.first << std::endl;

        std::cout << "PK1" << std::endl;
        std::cout << PK1 << std::endl;

        std::cout << "PK2" << std::endl;
        std::cout << convertPK1toPK2(PK1, F_curr) << std::endl;

        if(!tangentmodulus)
        {
            return std::make_pair(PK1, behavior3D_hpp.second);
        }

        const auto A = transpose(Pn) * computeContractedProduct(behavior3D_hpp.second, Pn) + projector.second;
        // std::cout << "A" << std::endl;
        // std::cout << A << std::endl;

        return std::make_pair(PK1, A);
    }

  public:
    LogarithmicStrain_qp(const point<scalar_type, DIM>& point, const scalar_type& weight)
    {
        m_law_hpp_qp = law_hpp_qp_type(point, weight);
        Pn           = static_tensor<scalar_type, 3>::Zero();
    }

    point<scalar_type, DIM>
    point() const
    {
        return m_law_hpp_qp.point();
    }

    scalar_type
    weight() const
    {
        return m_law_hpp_qp.weight();
    }

    bool
    is_plastic() const
    {
        return m_law_hpp_qp.is_plastic();
    }

    static_matrix_type3D
    getElasticStrain() const
    {
        return m_law_hpp_qp.getElasticStrain();
    }

    static_matrix_type3D
    getPlasticStrain() const
    {
        return m_law_hpp_qp.getPlasticStrain();
    }

    static_matrix_type
    getTotalStrain() const
    {
        return m_law_hpp_qp.getTotalStrain();
    }

    static_matrix_type
    getTotalStrainPrev() const
    {
        return m_law_hpp_qp.getTotalStrainPrev();
    }

    scalar_type
    getAccumulatedPlasticStrain() const
    {
        return m_law_hpp_qp.getAccumulatedPlasticStrain();
    }

    void
    update()
    {
        m_law_hpp_qp.update();
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(this->compute_stress3D(data));
    }

    static_matrix_type
    compute_stress_T(const data_type& data) const
    {
        return m_law_hpp_qp.compute_stress(data);
    }

    std::pair<static_matrix_type, tensor_type>
    compute_whole(const static_matrix_type& F_curr, const data_type& data, bool tangentmodulus = true)
    {
        const static_matrix_type   ID          = 2 * static_matrix_type::Identity();
        const static_matrix_type3D F_curr_3D   = convertMatrix3DwithOne(F_curr);
        const auto                 behaviors3D = compute_whole3D(F_curr_3D, data, tangentmodulus);

        const static_matrix_type              PK1 = convertMatrix<scalar_type, DIM>(behaviors3D.first);
        const static_tensor<scalar_type, DIM> A   = convertTensor<scalar_type, DIM>(behaviors3D.second);

        return std::make_pair(PK1,A);
        //return m_law_hpp_qp.compute_whole(F_curr, data, tangentmodulus);
    }
};
}
}