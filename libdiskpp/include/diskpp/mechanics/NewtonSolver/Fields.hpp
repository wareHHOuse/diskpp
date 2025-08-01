/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet (C) 2024
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

#include <vector>

#include "diskpp/adaptivity/adaptivity.hpp"

namespace disk
{

    namespace mechanics
    {

        enum FieldName
        {
            DEPL,
            DEPL_CELLS,
            DEPL_FACES,
            VITE_CELLS,
            ACCE_CELLS,
        };

        std::string
        getFieldName(const FieldName &name)
        {
            switch (name)
            {
            case FieldName::DEPL:
                return "DEPL";
                break;
            case FieldName::DEPL_CELLS:
                return "DEPL_CELLS";
                break;
            case FieldName::DEPL_FACES:
                return "DEPL_FACES";
                break;
            case FieldName::VITE_CELLS:
                return "VITE_CELLS";
                break;
            case FieldName::ACCE_CELLS:
                return "ACCE_CELLS";
                break;
            default:
                break;
            }

            return "Unknown type";
        }

        /**
         * @brief Field at one time.
         *
         * @tparam scalar_type
         */
        template <typename scalar_type>
        class TimeField
        {
            typedef dynamic_vector<scalar_type> vector_type;

            std::map<FieldName, std::vector<vector_type>> m_fields;

            scalar_type m_time;

        public:
            TimeField() : m_time(0) {}

            /**
             * @brief Get the current time
             *
             */
            scalar_type
            getTime(void) const
            {
                return m_time;
            }

            /**
             * @brief Set the current time
             *
             */
            void
            setTime(const scalar_type &time)
            {
                m_time = time;
            }

            void
            setField(const FieldName &name, const std::vector<vector_type> &field)
            {
                m_fields[name] = field;
            }

            auto
            getField(const FieldName &name) const
            {
                // std::cout << getFieldName(name) << std::endl;
                return m_fields.at(name);
            }

            void
            clearField(const FieldName &name)
            {
                if (auto search = m_fields.find(name); search != m_fields.end())
                {
                    m_fields.erase(search);
                }
            }

            template <typename Mesh>
            void createZeroField(const FieldName &name, const Mesh &mesh, const MeshDegreeInfo<Mesh> &degree_infos)
            {
                std::vector<vector_type> field;

                if (name == DEPL || name == DEPL_CELLS || name == VITE_CELLS || name == ACCE_CELLS)
                {
                    field.reserve(mesh.cells_size());

                    for (auto &cl : mesh)
                    {
                        const auto di = degree_infos.cellDegreeInfo(mesh, cl);
                        const auto num_cell_dofs = vector_basis_size(di.cell_degree(), Mesh::dimension, Mesh::dimension);
                        size_t num_faces_dofs = 0;

                        if (name == DEPL)
                        {
                            const auto fcs = faces(mesh, cl);
                            const auto fcs_di = di.facesDegreeInfo();

                            for (auto &fc_di : fcs_di)
                            {
                                if (fc_di.hasUnknowns())
                                {
                                    num_faces_dofs += vector_basis_size(fc_di.degree(), Mesh::dimension - 1, Mesh::dimension);
                                }
                            }
                        }

                        field.push_back(vector_type::Zero(num_cell_dofs + num_faces_dofs));
                    }
                }
                else if (name == DEPL_FACES)
                {
                    field.reserve(mesh.faces_size());

                    for (auto itor = mesh.faces_begin(); itor != mesh.faces_end(); itor++)
                    {
                        const auto fc = *itor;
                        const auto di = degree_infos.degreeInfo(mesh, fc);

                        size_t num_face_dofs = 0;
                        if (di.hasUnknowns())
                        {
                            num_face_dofs = vector_basis_size(di.degree(), Mesh::dimension - 1, Mesh::dimension);
                        }

                        field.push_back(vector_type::Zero(num_face_dofs));
                    }
                }
                else
                {
                    throw std::runtime_error("Unknown field");
                }

                this->setField(name, field);
            }

            template <typename Mesh>
            void createField(const FieldName &name, const Mesh &mesh, const MeshDegreeInfo<Mesh> &degree_infos, const vector_rhs_function<Mesh> func)
            {
                std::vector<vector_type> field;

                if (name == DEPL || name == DEPL_CELLS || name == VITE_CELLS || name == ACCE_CELLS)
                {
                    field.reserve(mesh.cells_size());

                    for (auto &cl : mesh)
                    {
                        if (name == DEPL)
                        {
                            field.push_back(project_function(mesh, cl, degree_infos, func, 2));
                        }
                        else
                        {
                            const auto di = degree_infos.cellDegreeInfo(mesh, cl);

                            field.push_back(project_function(mesh, cl, di.cell_degree(), func, 2));
                        }
                    }
                }
                else if (name == DEPL_FACES)
                {
                    field.reserve(mesh.faces_size());

                    for (auto itor = mesh.faces_begin(); itor != mesh.faces_end(); itor++)
                    {
                        const auto fc = *itor;
                        const auto fdi = degree_infos.degreeInfo(mesh, fc);

                        field.push_back(project_function(mesh, fc, fdi.degree(), func, 2));
                    }
                }
                else
                {
                    throw std::runtime_error("Unknown field");
                }

                this->setField(name, field);
            }
        };

        template <typename scalar_type>
        class MultiTimeField
        {
        private:
            typedef dynamic_vector<scalar_type> vector_type;
            typedef TimeField<scalar_type> field_type;

            std::vector<field_type> m_fields;

        public:
            MultiTimeField() {};

            MultiTimeField(const int n_fields)
            {
                m_fields.resize(n_fields);
            };

            void
            setCurrentTime(const scalar_type& time)
            {
                m_fields.at(0).setTime(time);
            }

            field_type getCurrentTimeField() const
            {
                return this->getTimeField(0);
            }

            field_type getPreviousTimeField() const
            {
                return this->getTimeField(-1);
            }

            field_type getTimeField(const int &relative_index) const
            {
                return m_fields.at(-relative_index);
            }

            void setTimeField(const int &relative_index, const field_type &field)
            {
                m_fields.at(-relative_index) = field;
            }

            void setCurrentTimeField(const field_type &field)
            {
                return this->setTimeField(0, field);
            }

            auto
            getField(const int &relative_index, const FieldName &name) const
            {
                return m_fields.at(-relative_index).getField(name);
            }

            auto
            getCurrentField(const FieldName &name) const
            {
                return m_fields.at(0).getField(name);
            }

            void
            setField(const int &relative_index, const FieldName &name, const std::vector<vector_type> &field)
            {
                return m_fields.at(-relative_index).setField(name, field);
            }

            void setCurrentField(const FieldName &name, const std::vector<vector_type> &field)
            {
                return this->setField(0, name, field);
            }

            template <typename Mesh>
            void createField(const int &relative_index, const FieldName &name, const Mesh &mesh, const MeshDegreeInfo<Mesh> &degree_infos, const vector_rhs_function<Mesh> func)
            {
                m_fields.at(-relative_index).createField(name, mesh, degree_infos, func);
            }

            auto getNumberOfTimeField() const
            {
                return m_fields.size();
            }

            void update()
            {

                int n_fields = this->getNumberOfTimeField();

                for (int i = n_fields - 1; i > 0; i--)
                {
                    m_fields.at(i) = m_fields.at(i - 1);
                }
            }
        };
    }

} // end disk
