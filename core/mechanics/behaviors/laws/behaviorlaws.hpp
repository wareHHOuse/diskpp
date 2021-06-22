/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018 - 2021                nicolas.pignet@enpc.fr
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

#include "Cavitation/Cavitation_qp.hpp"
#include "HenckyMises/HenckyMises_qp.hpp"
#include "IsotropicHardeningVMis/IsotropicHardeningVMis_qp.hpp"
#include "LinearElasticityLaw/LinearElasticityLaw.hpp"
#include "LinearIsotropicAndKinematicHardening/LinearIsotropicAndKinematicHardening_qp.hpp"
#include "LinearLaw/LinearLaw_qp.hpp"
#include "Neohookean/Neohookean_qp.hpp"
#include "mechanics/behaviors/logarithmic_strain/LogarithmicStrain.hpp"
#include "Mfront/Mfront_law.hpp"
#include "core/methods/hho"
#include "law_bones.hpp"
#include "behaviorlaws_names.hpp"

#ifdef HAVE_MGIS
#include "MGIS/Behaviour/Behaviour.hxx"
#endif

namespace disk
{

template<typename MeshType>
using Cavitation =
  LawTypeBones<MeshType, Cavitation_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using Neohookean =
  LawTypeBones<MeshType, Neohookean_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using HenckyMises =
  LawTypeBones<MeshType, HenckyMises_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using LinearElasticityLaw =
  LawTypeBones<MeshType, LinearElasticity_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using LinearLaw = LawTypeBones<MeshType, LinearLaw_qp<typename MeshType::coordinate_type, MeshType::dimension>, false>;

template<typename MeshType>
using LinearIsotropicAndKinematicHardening =
  LawTypeBones<MeshType,
               LinearIsotropicAndKinematicHardening_qp<typename MeshType::coordinate_type, MeshType::dimension>,
               true>;

template<typename MeshType>
using IsotropicHardeningVMis =
  LawTypeBones<MeshType, IsotropicHardeningVMis_qp<typename MeshType::coordinate_type, MeshType::dimension>, true>;

#ifdef HAVE_MGIS
template<typename MeshType>
using Mfront = Mfront_law<MeshType>;
#endif

template<typename MeshType>
class Behavior
{
  private:
    typedef typename MeshType::coordinate_type scalar_type;
    typedef typename MeshType::cell            cell_type;
    typedef MaterialData<scalar_type>          material_type;
    typedef dynamic_vector<scalar_type>        vector_type;

    typedef static_matrix<scalar_type, MeshType::dimension, MeshType::dimension>  static_matrix_type;
    typedef static_tensor<scalar_type, MeshType::dimension>                       static_tensor_type;

    size_t m_deformation;
    size_t m_law;
    size_t m_id;

    material_type m_data;

#ifdef HAVE_MGIS
    typedef std::shared_ptr<mgis::behaviour::Behaviour> BehaviourPtr;
    BehaviourPtr                                        m_behav;
#endif

    Cavitation<MeshType>                           m_cavitation;
    Neohookean<MeshType>                           m_neohokean;
    HenckyMises<MeshType>                          m_henckymises;
    LinearElasticityLaw<MeshType>                  m_elastic;
    LinearIsotropicAndKinematicHardening<MeshType> m_linearHard;
    IsotropicHardeningVMis<MeshType>               m_nonlinearHard;
//
    mechanics::LogarithmicStrain<LinearElasticityLaw<MeshType>>                  m_log_elastic;
    mechanics::LogarithmicStrain<LinearIsotropicAndKinematicHardening<MeshType>> m_log_linearHard;
    mechanics::LogarithmicStrain<IsotropicHardeningVMis<MeshType>>               m_log_nonlinearHard;

#ifdef HAVE_MGIS
    Mfront<MeshType>                               m_mfront;
#endif

    void
    select_law(void)
    {
        m_id = 0;
        switch (m_deformation)
        {
            case DeformationMeasure::SMALL_DEF:
                switch (m_law)
                {
                    case LawType::ELASTIC: m_id = 100; break;
                    case LawType::LINEAR_HARDENING: m_id = 101; break;
                    case LawType::NONLINEAR_HARDENING: m_id = 102; break;
                    case LawType::HENCKY_MISES: m_id = 103; break;
#ifdef HAVE_MGIS
                    case LawType::MFRONT: m_id = 500; break;
#endif
                    default: throw std::invalid_argument("Incompatible law with SMALL_DEF");
                }
                break;
            case DeformationMeasure::F_DEF:
                switch (m_law)
                {
                    case LawType::NEOHOKEAN: m_id = 200; break;
                    case LawType::CAVITATION: m_id = 201; break;
#ifdef HAVE_MGIS
                    case LawType::MFRONT: m_id = 500; break;
#endif
                    default: throw std::invalid_argument("Incompatible law with F_DEF");
                }
                break;
            case DeformationMeasure::LOGARITHMIC_DEF:
                switch (m_law)
                {
                    case LawType::ELASTIC: m_id = 300; break;
                    case LawType::LINEAR_HARDENING: m_id = 301; break;
                    case LawType::NONLINEAR_HARDENING: m_id = 302; break;
#ifdef HAVE_MGIS
                    case LawType::MFRONT: m_id = 500; break;
#endif
                    default: throw std::invalid_argument("Incompatible law with LOGARITHMIC_DEF");
                }
                break;

            default: throw std::invalid_argument("Unknown deformation");
        }
    }

  public:
    Behavior() : m_deformation(DeformationMeasure::SMALL_DEF), m_law(LawType::ELASTIC) { select_law(); }

    Behavior(const size_t deformation, const size_t law) : m_deformation(deformation), m_law(law) { select_law(); }

    Behavior(const MeshType& msh, const size_t degree, const size_t deformation, const size_t law) :
      m_deformation(deformation), m_law(law)
    {
        select_law();
        switch (m_id)
        {
            case 100: m_elastic = LinearElasticityLaw<MeshType>(msh, degree); break;
            case 101: m_linearHard = LinearIsotropicAndKinematicHardening<MeshType>(msh, degree); break;
            case 102: m_nonlinearHard = IsotropicHardeningVMis<MeshType>(msh, degree); break;
            case 103: m_henckymises = HenckyMises<MeshType>(msh, degree); break;
            case 200: m_neohokean = Neohookean<MeshType>(msh, degree); break;
            case 201: m_cavitation = Cavitation<MeshType>(msh, degree); break;
//
            case 300: m_log_elastic = mechanics::LogarithmicStrain<LinearElasticityLaw<MeshType>>(msh, degree); break;
            case 301:
                m_log_linearHard =
                  mechanics::LogarithmicStrain<LinearIsotropicAndKinematicHardening<MeshType>>(msh, degree);
                break;
            case 302:
                m_log_nonlinearHard =
                  mechanics::LogarithmicStrain<IsotropicHardeningVMis<MeshType>>(msh, degree);
                break;

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

#ifdef HAVE_MGIS

    Behavior(const MeshType&                   msh,
             const size_t                      degree,
             const std::string&                filename,
             const std::string&                law,
             const mgis::behaviour::Hypothesis h) :
      m_law(LawType::MFRONT)
    {
        using namespace mgis::behaviour;
        std::cout << "Loading MFRONT law: " << law << std::endl;
        if (isStandardFiniteStrainBehaviour(filename, law))
        {
            m_deformation         = DeformationMeasure::F_DEF;
            auto opts             = FiniteStrainBehaviourOptions{};
            opts.stress_measure   = FiniteStrainBehaviourOptions::PK1;
            opts.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
            m_behav               = std::make_shared<Behaviour>(load(opts, filename, law, h));
        }
        else
        {
            m_deformation = DeformationMeasure::SMALL_DEF;
            m_behav       = std::make_shared<Behaviour>(load(filename, law, h));
        }

        std::cout << "Behaviour type: " << (*m_behav).btype << std::endl;
        std::cout << "Kinematic: " << (*m_behav).kinematic << std::endl;

        std::cout << "Material properties: (name, type)" << std::endl;
        for (const auto& mp : (*m_behav).mps)
            std::cout << mp.name << ", " << MfrontVariableTypeName(mp.type) << std::endl;

        std::cout << "Internal State Variables: (name, type)" << std::endl;
        for (const auto& is : (*m_behav).isvs)
            std::cout << is.name << ", " << MfrontVariableTypeName(is.type) << std::endl;

        std::cout << "Thermodynamic forces: (name, type)" << std::endl;
        for (const auto& fc : (*m_behav).thermodynamic_forces)
            std::cout << fc.name << ", " << MfrontVariableTypeName(fc.type) << std::endl;

        select_law();
        switch (m_id)
        {
            case 500: m_mfront = Mfront<MeshType>(msh, degree, m_behav); break;
            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }
#endif

    size_t getDeformation(void) const
    {
        return m_deformation;
    }

    std::string
    getDeformationName(void) const
    {
        return DeformationMeasureName( m_deformation);
    }

    std::string
    getLawName(void) const
    {
        return LawTypeName(m_law);
    }

    void
    addMaterialData(const material_type& materialData)
    {
        m_data = materialData;
        switch (m_id)
        {
            case 100: m_elastic.addMaterialData(materialData); break;
            case 101: m_linearHard.addMaterialData(materialData); break;
            case 102: m_nonlinearHard.addMaterialData(materialData); break;
            case 103: m_henckymises.addMaterialData(materialData); break;
            case 200: m_neohokean.addMaterialData(materialData); break;
            case 201: m_cavitation.addMaterialData(materialData); break;
            case 300: m_log_elastic.addMaterialData(materialData); break;
            case 301: m_log_linearHard.addMaterialData(materialData); break;
            case 302: m_log_nonlinearHard.addMaterialData(materialData); break;
#ifdef HAVE_MGIS

            case 500: m_mfront.addMaterialData(materialData); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    material_type
    getMaterialData(void)
    {
        return m_data;
    }

    const material_type&
    getMaterialData(void) const
    {
        return m_data;
    }

    /**
     * @brief Get number of quadrature points for the used law
     *
     * @return size_t number of quadrature law
     */
    size_t
    numberOfQP(void) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getNumberOfQP(); break;
            case 101: return m_linearHard.getNumberOfQP(); break;
            case 102: return m_nonlinearHard.getNumberOfQP(); break;
            case 103: return m_henckymises.getNumberOfQP(); break;
            case 200: return m_neohokean.getNumberOfQP(); break;
            case 201: return m_cavitation.getNumberOfQP(); break;
            case 300: return m_log_elastic.getNumberOfQP(); break;
            case 301: return m_log_linearHard.getNumberOfQP(); break;
            case 302: return m_log_nonlinearHard.getNumberOfQP(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getNumberOfQP(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    size_t
    numberOfQP(const size_t& cell_id) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getNumberOfQP(); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getNumberOfQP(); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getNumberOfQP(); break;
            case 103: return m_henckymises.getCellQPs(cell_id).getNumberOfQP(); break;
            case 200: return m_neohokean.getCellQPs(cell_id).getNumberOfQP(); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getNumberOfQP(); break;
            case 300: return m_log_elastic.getCellQPs(cell_id).getNumberOfQP(); break;
            case 301: return m_log_linearHard.getCellQPs(cell_id).getNumberOfQP(); break;
            case 302: return m_log_nonlinearHard.getCellQPs(cell_id).getNumberOfQP(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getCellQPs(cell_id).getNumberOfQP(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    size_t
    numberOfQP(const MeshType& msh, const cell_type& cl) const
    {
        const auto cell_id = msh.lookup(cl);

        return numberOfQP(cell_id);
    }

    auto
    quadrature_point(const size_t& cell_id, const size_t& qp_id) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 103: return m_henckymises.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 200: return m_neohokean.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 300: return m_log_elastic.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 301: return m_log_linearHard.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
            case 302: return m_log_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getCellQPs(cell_id).getQP(qp_id).quadrature_point(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    std::pair<static_matrix_type, static_tensor_type>
    compute_whole(const size_t& cell_id, const size_t& qp_id, const static_matrix_type& RkT_iqn, bool tangent = true)
    {
        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent); break;
            case 103:
                return m_henckymises.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
            case 200: return m_neohokean.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
            case 300:
                return m_log_elastic.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
            case 301:
                return m_log_linearHard.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
            case 302:
                return m_log_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
#ifdef HAVE_MGIS
            case 500:
                return m_mfront.getCellQPs(cell_id).getQP(qp_id).compute_whole(RkT_iqn, m_data, tangent);
                break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    static_matrix<scalar_type, 3, 3>
    compute_stress3D(const size_t& cell_id, const size_t& qp_id) const
    {
        switch (m_id)
        {
            case 100:
                return m_elastic.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data);
                break;
            case 101:
                return m_linearHard.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data);
                break;
            case 102:
                return m_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data);
                break;
            case 103:
                return m_henckymises.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data);
                break;
            case 200: return m_neohokean.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data); break;
            case 201:
                return m_cavitation.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data);
                break;
            case 300: return m_log_elastic.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data); break;
            case 301: return m_log_linearHard.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data); break;
            case 302: return m_log_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getCellQPs(cell_id).getQP(qp_id).compute_stress3D(m_data); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    bool
    is_plastic(const size_t& cell_id, const size_t& qp_id) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 103: return m_henckymises.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 200: return m_neohokean.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 300: return m_log_elastic.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 301: return m_log_linearHard.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
            case 302: return m_log_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getCellQPs(cell_id).getQP(qp_id).is_plastic(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    scalar_type
    equivalentPlasticStrain(const size_t& cell_id, const size_t& qp_id) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 103: return m_henckymises.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 200: return m_neohokean.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 300: return m_log_elastic.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 301: return m_log_linearHard.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
            case 302: return m_log_nonlinearHard.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.getCellQPs(cell_id).getQP(qp_id).getEquivalentPlasticStrain(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    void update(void)
    {
        switch (m_id)
        {
            case 100: return m_elastic.update(); break;
            case 101: return m_linearHard.update(); break;
            case 102: return m_nonlinearHard.update(); break;
            case 103: return m_henckymises.update(); break;
            case 200: return m_neohokean.update(); break;
            case 201: return m_cavitation.update(); break;
            case 300: return m_log_elastic.update(); break;
            case 301: return m_log_linearHard.update(); break;
            case 302: return m_log_nonlinearHard.update(); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.update(); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    vector_type
    projectStressOnCell(const MeshType& msh, const cell_type& cl, const hho_degree_info& hdi) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.projectStressOnCell(msh, cl, hdi, m_data); break;
            case 101: return m_linearHard.projectStressOnCell(msh, cl, hdi, m_data); break;
            case 102: return m_nonlinearHard.projectStressOnCell(msh, cl, hdi, m_data); break;
            case 103: return m_henckymises.projectStressOnCell(msh, cl, hdi, m_data); break;
            case 200: return m_neohokean.projectStressOnCell(msh, cl, hdi, m_data); break;
            case 201: return m_cavitation.projectStressOnCell(msh, cl, hdi, m_data); break;
#ifdef HAVE_MGIS
            case 500: return m_mfront.projectStressOnCell(msh, cl, hdi, m_data); break;
#endif

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    // vector_type
    // projectPOnCell(const MeshType& msh, const cell_type& cl, const hho_degree_info& hdi) const
    // {
    //     switch (m_id)
    //     {
    //         case 100: return m_elastic.projectPOnCell(msh, cl, hdi, mate); break;
    //         case 101: return m_linearHard.projectPOnCell(msh, cl, hdi, mate); break;
    //         case 102: return m_nonlinearHard.projectPOnCell(msh, cl, hdi, mate); break;
    //         case 103: return m_henckymises.projectPOnCell(msh, cl, hdi, mate); break;
    //         case 200: return m_neohokean.projectPOnCell(msh, cl, hdi, mate); break;
    //         case 201: return m_cavitation.projectPOnCell(msh, cl, hdi, mate); break;

    //         default: throw std::invalid_argument("Behavior error: Unknown id law");
    //     }
    // }
};
}
