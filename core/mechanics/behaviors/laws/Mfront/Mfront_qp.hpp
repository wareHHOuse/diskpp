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

#ifdef HAVE_MGIS

#include <iostream>

#include "common/eigen.hpp"
#include "core/mechanics/behaviors/laws/law_qp_bones.hpp"
#include "core/mechanics/behaviors/laws/materialData.hpp"
#include "core/mechanics/behaviors/maths_tensor.hpp"
#include "core/mechanics/behaviors/maths_utils.hpp"
#include "core/mechanics/behaviors/tensor_conversion.hpp"
#include "mesh/point.hpp"

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace disk
{

// Law developped with Front interface
// see http://tfel.sourceforge.net/index.html

template<typename T, int DIM>
class Mfront_qp : public law_qp_bones<T, DIM>
{
  public:
    const static size_t                                 dimension = DIM;
    typedef T                                           scalar_type;
    typedef static_matrix<scalar_type, DIM, DIM>        static_matrix_type;
    typedef static_matrix<scalar_type, 3, 3>            static_matrix_type3D;
    typedef MaterialData<scalar_type>                   data_type;
    typedef std::shared_ptr<mgis::behaviour::Behaviour> BehaviourPtr;

  private:
    bool        l_small_def;
    std::string StressName;

    // shared pointer to mgis behaviour
    BehaviourPtr m_behav;

    // BehaviourDatat for Mfront computation
    mgis::behaviour::BehaviourData     m_behavData;
    mgis::behaviour::BehaviourDataView m_behavDataView;

  public:
    Mfront_qp() :
      law_qp_bones<T, DIM>(), m_behav(nullptr), m_behavData(*m_behav),
      m_behavDataView(mgis::behaviour::make_view(m_behavData))
    {
    }

    Mfront_qp(const point<scalar_type, DIM>& point, const scalar_type& weight, const BehaviourPtr& behav) :
      law_qp_bones<T, DIM>(point, weight), m_behav(behav), m_behavData(*m_behav),
      m_behavDataView(mgis::behaviour::make_view(m_behavData))
    {
        if ((*m_behav).kinematic == mgis::behaviour::Behaviour::SMALLSTRAINKINEMATIC)
        {
            l_small_def = true;
            StressName  = "Stress";
        }
        else if ((*m_behav).kinematic == mgis::behaviour::Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)
        {
            l_small_def = false;
            StressName  = "FirstPiolaKirchhoffStress";
        }
        else
            throw std::runtime_error("Error");
    }

    void
    update()
    {
        law_qp_bones<T, DIM>::update();
        mgis::behaviour::update(m_behavData);
    }

    scalar_type
    getEquivalentPlasticStrain() const
    {
        try
        {
            const auto ePSPtr = mgis::behaviour::getInternalStateVariable(m_behavData.s1, "EquivalentPlasticStrain");
            return *ePSPtr;
        }
        catch(...)
        {}

        return 0.0;
    }

    bool
    is_plastic() const
    {
        if (getEquivalentPlasticStrain() != 0.0)
            return true;

        return false;
    }

    void
    addMaterialParameters(const data_type& data)
    {
        const auto& mdata = data.getMfrontParameters();
        for (auto& [param, value] : mdata)
        {
            mgis::behaviour::setMaterialProperty(m_behavData.s1, param, value);
        }
    }

    void
    addInitialMaterialParameters(const data_type& data)
    {
        const auto& mdata = data.getMfrontParameters();
        for (auto& [param, value] : mdata)
        {
            mgis::behaviour::setMaterialProperty(m_behavData.s0, param, value);
        }
    }

    static_matrix_type3D
    compute_stress3D(const data_type& data) const
    {
        const auto stressPtr = mgis::behaviour::getThermodynamicForce(m_behavData.s1, StressName);
        if (stressPtr != &(m_behavData.s1.thermodynamic_forces[0]))
            throw std::runtime_error("We assume that the stress is the first variable");

        static_matrix_type3D stress;
        convertMatrixFromMgis(m_behavData.s1.thermodynamic_forces, stress);

        return stress;
    }

    static_matrix_type
    compute_stress(const data_type& data) const
    {
        return convertMatrix<scalar_type, DIM>(compute_stress3D(data));
    }

    std::pair<static_matrix_type3D, static_tensor<scalar_type, 3>>
    compute_whole3D(const static_matrix_type3D& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        this->m_estrain_curr = strain_curr;

        if (tangentmodulus)
        {
            m_behavData.K[0] = 4;
        }
        else
        {
            m_behavData.K[0] = 1;
        }

        // Output: PK1, d PK1 / d F
        if (!l_small_def)
        {
            m_behavData.K[1] = 2;
            m_behavData.K[2] = 2;
        }

        // std::cout << "K: " << m_behavData.K[0] << ", " << m_behavData.K[1] << ", " << m_behavData.K[2] << std::endl;

        // mfront
        convertMatrixToMgis(strain_curr, m_behavData.s1.gradients);
        auto v = mgis::behaviour::make_view(m_behavData);
        mgis::behaviour::integrate(v, *m_behav);

        // compute stress tensor depending on the choice
        const static_matrix_type3D stress = compute_stress3D(data);
        // compute tangent module (consistent with the stress tensor)

        static_tensor<scalar_type, 3> Aep;
        convertTensorFromMgis(m_behavData.K, Aep);

        // std::cout << "MGIS" << std::endl;
        // for (auto& elem : m_behavData.K)
        //     std::cout << elem << ", ";
        // std::cout << std::endl;
        // std::cout << "Aep" << std::endl;
        // std::cout << Aep << std::endl;

        // printing
        // print_markdown();

        return std::make_pair(stress, Aep);
    }

    // strain_curr is the symetric gradient for small deformation and F for finite def.
    std::pair<static_matrix_type, static_tensor<scalar_type, DIM>>
    compute_whole(const static_matrix_type& strain_curr, const data_type& data, bool tangentmodulus = true)
    {
        static_matrix_type3D strain3D_curr;
        if (l_small_def)
        {
            strain3D_curr = convertMatrix3D(strain_curr);
        }
        else
        {
            strain3D_curr = convertMatrix3DwithOne(strain_curr);
        }
        const auto behaviors3D = this->compute_whole3D(strain3D_curr, data, tangentmodulus);

        const static_matrix_type              stress = convertMatrix<scalar_type, DIM>(behaviors3D.first);
        const static_tensor<scalar_type, DIM> Cep    = convertTensor<scalar_type, DIM>(behaviors3D.second);

        return std::make_pair(stress, Cep);
    }

    void
    print_markdown() const
    {
        mgis::behaviour::print_markdown(std::cout, *m_behav, m_behavData, 1);
    }
};
#endif
}
