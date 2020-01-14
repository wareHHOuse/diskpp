/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2018, 2019                nicolas.pignet@enpc.fr
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
#include "core/methods/hho"
#include "law_bones.hpp"

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

enum DeformationMeasure : size_t
{
    SMALL_DEF       = 0,
    LOGARITHMIC_DEF = 1,
    F_DEF           = 2,
};

enum LawType : size_t
{
    ELASTIC             = 0,
    LINEAR_HARDENING    = 1,
    NONLINEAR_HARDENING = 2,
    HENCKY_MISES        = 3,
    NEOHOKEAN           = 4,
    CAVITATION          = 5
};

template<typename MeshType>
class Behavior
{
  private:
    typedef typename MeshType::coordinate_type scalar_type;
    typedef typename MeshType::cell            cell_type;
    typedef MaterialData<scalar_type>          material_type;
    typedef dynamic_vector<scalar_type>        vector_type;

    size_t m_deformation;
    size_t m_law;
    size_t m_id;

    Cavitation<MeshType>                           m_cavitation;
    Neohookean<MeshType>                           m_neohokean;
    HenckyMises<MeshType>                          m_henckymises;
    LinearElasticityLaw<MeshType>                  m_elastic;
    LinearIsotropicAndKinematicHardening<MeshType> m_linearHard;
    IsotropicHardeningVMis<MeshType>               m_nonlinearHard;

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
                    default: throw std::invalid_argument("Incompatible law with SMALL_DEF");
                }
                break;
            case DeformationMeasure::F_DEF:
                switch (m_law)
                {
                    case LawType::NEOHOKEAN: m_id = 200; break;
                    case LawType::CAVITATION: m_id = 201; break;
                    default: throw std::invalid_argument("Incompatible law with F_DEF");
                }
                break;
            case DeformationMeasure::LOGARITHMIC_DEF:
                switch (m_law)
                {
                    // case LawType::ELASTIC: m_id = 300; break;
                    // case LawType::LINEAR_HARDENING: m_id = 301; break;
                    // case LawType::NONLINEAR_HARDENING: m_id = 302; break;
                    // case LawType::HENCKY_MISES: m_id = 303; break;
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

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    size_t
    getDeformation(void) const
    {
        return m_deformation;
    }

    void
    addMaterialData(const material_type& materialData)
    {
        switch (m_id)
        {
            case 100: m_elastic.addMaterialData(materialData); break;
            case 101: m_linearHard.addMaterialData(materialData); break;
            case 102: m_nonlinearHard.addMaterialData(materialData); break;
            case 103: m_henckymises.addMaterialData(materialData); break;
            case 200: m_neohokean.addMaterialData(materialData); break;
            case 201: m_cavitation.addMaterialData(materialData); break;

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    material_type
    getMaterialData(void) const
    {
        switch (m_id)
        {
            case 100: return m_elastic.getMaterialData(); break;
            case 101: return m_linearHard.getMaterialData(); break;
            case 102: return m_nonlinearHard.getMaterialData(); break;
            case 103: return m_henckymises.getMaterialData(); break;
            case 200: return m_neohokean.getMaterialData(); break;
            case 201: return m_cavitation.getMaterialData(); break;

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    // const material_type&
    // getMaterialData(void) const
    // {
    //     switch (m_id)
    //     {
    //         case 100: return m_linear.getMaterialData(); break;
    //         case 101: return m_linearHard.getMaterialData(); break;
    //         case 102: return m_nonlinearHard.getMaterialData(); break;
    //         case 103: return m_henckymises.getMaterialData(); break;
    //         case 200: return m_neohokean.getMaterialData(); break;
    //         case 201: return m_cavitation.getMaterialData(); break;

    //         default: throw std::invalid_argument("Behavior error: Unknown id law");
    //     }
    // }

    auto&
    law()
    {
        return m_elastic;
    }

    std::vector<law_qp_bones<scalar_type, MeshType::dimension>>&
    getQPs(const MeshType& msh, const cell_type& cl)
    {
        const auto cell_id = msh.lookup(cl);

        return m_elastic.getCellQPs(cell_id).getQPs();
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

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    size_t
    numberOfQP(const MeshType& msh, const cell_type& cl) const
    {
        const auto cell_id = msh.lookup(cl);

        switch (m_id)
        {
            case 100: return m_elastic.getCellQPs(cell_id).getNumberOfQP(); break;
            case 101: return m_linearHard.getCellQPs(cell_id).getNumberOfQP(); break;
            case 102: return m_nonlinearHard.getCellQPs(cell_id).getNumberOfQP(); break;
            case 103: return m_henckymises.getCellQPs(cell_id).getNumberOfQP(); break;
            case 200: return m_neohokean.getCellQPs(cell_id).getNumberOfQP(); break;
            case 201: return m_cavitation.getCellQPs(cell_id).getNumberOfQP(); break;

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    void
    update(void)
    {
        switch (m_id)
        {
            case 100: return m_elastic.update(); break;
            case 101: return m_linearHard.update(); break;
            case 102: return m_nonlinearHard.update(); break;
            case 103: return m_henckymises.update(); break;
            case 200: return m_neohokean.update(); break;
            case 201: return m_cavitation.update(); break;

            default: throw std::invalid_argument("Behavior error: Unknown id law");
        }
    }

    vector_type
    projectStressOnCell(const MeshType& msh, const cell_type& cl, const hho_degree_info& hdi) const
    {
        const material_type& mate = this->getMaterialData();
        switch (m_id)
        {
            case 100: return m_elastic.projectStressOnCell(msh, cl, hdi, mate); break;
            case 101: return m_linearHard.projectStressOnCell(msh, cl, hdi, mate); break;
            case 102: return m_nonlinearHard.projectStressOnCell(msh, cl, hdi, mate); break;
            case 103: return m_henckymises.projectStressOnCell(msh, cl, hdi, mate); break;
            case 200: return m_neohokean.projectStressOnCell(msh, cl, hdi, mate); break;
            case 201: return m_cavitation.projectStressOnCell(msh, cl, hdi, mate); break;

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
