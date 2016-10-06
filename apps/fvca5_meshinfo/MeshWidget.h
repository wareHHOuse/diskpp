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

        int xs = int((pos.x()*(size_x/4))/dist.x());
        int ys = int((pos.y()*(size_y/4))/dist.y());

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
