/*
 *       /\        Matteo Cicuttin (C) 2016, 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
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

#include <iostream>
#include <cassert>

#include <QtGui>

#include "mesh/mesh.hpp"

template<typename Mesh>
class MeshWidget : public QWidget
{
    typedef Mesh                            mesh_type;
    typedef typename mesh_type::point_type  point_type;

    mesh_type                       m_msh;
    disk::bounding_box<mesh_type>   m_bbox;

    int size_x, size_y;

    QPoint transform(const point_type& pt)
    {
        auto pos = pt + m_bbox.min();
        auto dist = m_bbox.max() - m_bbox.min();

        int xs = int((pos.x()*(size_x/2))/dist.x())+10;
        int ys = size_x - int((pos.y()*(size_y/2))/dist.y())-10;

        return QPoint(xs, ys);
    }

public:
    MeshWidget(QWidget *parent = nullptr) : QWidget(parent)
    {
        size_x = 600;
        size_y = 600;

        resize(size_x, size_y);
    }

    void
    setMesh(const mesh_type& msh)
    {
        m_msh = msh;
        m_bbox = disk::bounding_box<mesh_type>(m_msh);
        std::cout << m_bbox << std::endl;
    }

    virtual void
    paintEvent(QPaintEvent *event)
    {
        QPainter    painter(this);

        for (auto itor = m_msh.faces_begin(); itor != m_msh.faces_end(); itor++)
        {
            auto fc = *itor;
            auto pts = points(m_msh, fc);
            assert(pts.size() == 2);
            auto qpt1 = transform(pts[0]);
            auto qpt2 = transform(pts[1]);
            painter.drawLine(qpt1, qpt2);
        }
    }
};
