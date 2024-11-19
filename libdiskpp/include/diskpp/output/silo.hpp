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

/*
 * This source file is part of EMT, the ElectroMagneticTool.
 *
 * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#pragma once


#include <Eigen/src/Core/util/Constants.h>
#include <silo.h>
#include <complex>
#include <algorithm>
#include <vector>

namespace disk {

namespace priv {
template<template<typename, size_t, typename> class Mesh,
    typename T, typename Storage>
int
silo_face_orientation(const Mesh<T, 3, Storage>& msh,
                      const typename Mesh<T, 3, Storage>::cell& cl,
                      const typename Mesh<T, 3, Storage>::face& fc)
{
    auto pts = points(msh, fc);
    assert(pts.size() >= 3);

    auto v0 = (pts[1] - pts[0]).to_vector();
    auto v1 = (pts[2] - pts[1]).to_vector();
    auto n = v0.cross(v1);

    auto cell_bar = barycenter(msh, cl);
    auto face_bar = barycenter(msh, fc);
    auto outward_vector = (face_bar - cell_bar).to_vector();

    if ( n.dot(outward_vector) < T(0) )
        return -1;

    return 1;
}
}


enum variable_centering_t
{
    nodal_variable_t,
    zonal_variable_t
};

template<typename T, variable_centering_t centering>
class silo_variable
{
    std::string     m_name;
    std::vector<T>  m_values;

public:

    typedef typename std::vector<T>::iterator           iterator;
    typedef typename std::vector<T>::const_iterator     const_iterator;

    silo_variable()
    {}

    silo_variable(const std::string& name, const std::vector<T>& values)
    {
        m_name = name;
        m_values.resize( values.size() );
        std::copy(values.begin(), values.end(), m_values.begin());
    }

    silo_variable(const std::string& name, std::vector<T>&& values)
    {
        m_name = name;
        m_values = std::move(values);
    }

    size_t size(void) const {
        return m_values.size();
    }

    std::string name(void) const {
        return m_name;
    }

    T* data(void) {
        return m_values.data();
    }

    iterator begin() {
        return m_values.begin();
    }

    iterator end() {
        return m_values.end();
    }

    const_iterator begin() const {
        return m_values.begin();
    }

    const_iterator end() const {
        return m_values.end();
    }

};

template<typename T>
using silo_nodal_variable = silo_variable<T, nodal_variable_t>;

template<typename T>
using silo_zonal_variable = silo_variable<T, zonal_variable_t>;

class silo_database
{
    DBfile          *m_siloDb;

public:
    silo_database()
        : m_siloDb(nullptr)
    {}

    bool create(const std::string& db_name)
    {
        m_siloDb = DBCreate(db_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        if (m_siloDb)
            return true;

        std::cout << "Error creating database" << std::endl;
        return false;
    }

    bool open(const std::string& db_name)
    {
        m_siloDb = DBOpen(db_name.c_str(), DB_PDB, DB_APPEND);
        if (m_siloDb)
            return true;

        std::cout << "Error opening database" << std::endl;
        return false;
    }

    bool close()
    {
        if (m_siloDb)
            DBClose(m_siloDb);
        m_siloDb = NULL;
        return true;
    }

    ~silo_database()
    {
        if (m_siloDb)
            DBClose(m_siloDb);
    }

    template<typename T>
    bool add_mesh(const simplicial_mesh<T,3>& msh, const std::string& name)
    {
        std::vector<double> x_coords, y_coords, z_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());
        z_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            /* explicit cast to support rational<>. silo works in double anyways */
            x_coords.push_back( double(pt.x()) );
            y_coords.push_back( double(pt.y()) );
            z_coords.push_back( double(pt.z()) );
        }

        double *coords[] = {x_coords.data(), y_coords.data(), z_coords.data()};

        std::vector<int> nodelist;
        nodelist.reserve( 4*msh.cells_size() );

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            assert(ptids.size() == 4);

            for (auto& ptid : ptids)
                nodelist.push_back( ptid + 1 );
        }

        int lnodelist = nodelist.size();

        int shapetype[] = { DB_ZONETYPE_TET };
        int shapesize[] = {4};
        int shapecounts[] = { static_cast<int>(msh.cells_size()) };
        int nshapetypes = 1;
        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 3;

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        DBPutZonelist2(m_siloDb, zonelist_name.c_str(), nzones, ndims,
            nodelist.data(), lnodelist, 1, 0, 0, shapetype, shapesize,
            shapecounts, nshapetypes, NULL);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            zonelist_name.c_str(), NULL, DB_DOUBLE, NULL);

        return true;
    }

    template<typename T>
    bool add_mesh(const std::vector<point<T,2>>& pts, const std::string& name)
    {
        std::cout << "Point mesh: " << pts.size() << std::endl;
        std::vector<double> x_coords, y_coords;
        x_coords.reserve(pts.size());
        y_coords.reserve(pts.size());

        for (auto& pt : pts) 
        {
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
        }

        double *coords[] = { x_coords.data(), y_coords.data() };

        int ndims = 2;
        int nels = pts.size();

        auto err = DBPutPointmesh(m_siloDb, name.c_str(), ndims, coords, nels, DB_DOUBLE, nullptr);
        return err == 0;
    }

    template<typename T>
    bool add_mesh(const simplicial_mesh<T,2>& msh, const std::string& name)
    {
        std::vector<double> x_coords, y_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            /* explicit cast to support rational<>. silo works in double anyways */
            x_coords.push_back( double(pt.x()) );
            y_coords.push_back( double(pt.y()) );
        }

        double *coords[] = {x_coords.data(), y_coords.data()};

        std::vector<int> nodelist;
        nodelist.reserve( 3*msh.cells_size() );

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            assert(ptids.size() == 3);

            for (auto& ptid : ptids)
                nodelist.push_back( ptid + 1 );
        }

        int lnodelist = nodelist.size();

        int shapetype[] = { DB_ZONETYPE_TRIANGLE };
        int shapesize[] = {3};
        int shapecounts[] = { static_cast<int>(msh.cells_size()) };
        int nshapetypes = 1;
        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 2;

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        DBPutZonelist2(m_siloDb, zonelist_name.c_str(), nzones, ndims,
            nodelist.data(), lnodelist, 1, 0, 0, shapetype, shapesize,
            shapecounts, nshapetypes, NULL);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            zonelist_name.c_str(), NULL, DB_DOUBLE, NULL);

        return true;
    }

    template<typename T>
    bool add_mesh(const cartesian_mesh<T,3>& msh, const std::string& name)
    {
        std::vector<double> x_coords, y_coords, z_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());
        z_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            /* explicit cast to support rational<>. silo works in double anyways */
            x_coords.push_back( double(pt.x()) );
            y_coords.push_back( double(pt.y()) );
            z_coords.push_back( double(pt.z()) );
        }

        double *coords[] = {x_coords.data(), y_coords.data(), z_coords.data()};

        std::vector<int> nodelist;
        nodelist.reserve( 8*msh.cells_size() );

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            assert(ptids.size() == 8);

            nodelist.push_back( ptids[0] + 1 );
            nodelist.push_back( ptids[1] + 1 );
            nodelist.push_back( ptids[3] + 1 );
            nodelist.push_back( ptids[2] + 1 );
            nodelist.push_back( ptids[4] + 1 );
            nodelist.push_back( ptids[5] + 1 );
            nodelist.push_back( ptids[7] + 1 );
            nodelist.push_back( ptids[6] + 1 );
        }

        int lnodelist = nodelist.size();

        int shapetype[] = { DB_ZONETYPE_HEX };
        int shapesize[] = {8};
        int shapecounts[] = { static_cast<int>(msh.cells_size()) };
        int nshapetypes = 1;
        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 3;

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        DBPutZonelist2(m_siloDb, zonelist_name.c_str(), nzones, ndims,
            nodelist.data(), lnodelist, 1, 0, 0, shapetype, shapesize,
            shapecounts, nshapetypes, NULL);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            zonelist_name.c_str(), NULL, DB_DOUBLE, NULL);

        return true;
    }

    template<typename T>
    bool add_mesh(const cartesian_mesh<T,2>& msh, const std::string& name)
    {
        std::vector<T> x_coords, y_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
        }

        T *coords[] = {x_coords.data(), y_coords.data()};

        std::vector<int> nodelist;
        nodelist.reserve( 4*msh.cells_size() );

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            assert(ptids.size() == 4);

            nodelist.push_back( ptids[0]+1 );
            nodelist.push_back( ptids[1]+1 );
            nodelist.push_back( ptids[3]+1 );
            nodelist.push_back( ptids[2]+1 );
        }

        int lnodelist = nodelist.size();

        int shapesize[] = {4};
        int shapecounts[] = { static_cast<int>(msh.cells_size()) };
        int nshapetypes = 1;
        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 2;

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        DBPutZonelist(m_siloDb, zonelist_name.c_str(), nzones, ndims,
                      nodelist.data(), lnodelist, 1, shapesize, shapecounts, nshapetypes);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            zonelist_name.c_str(), NULL, DB_DOUBLE, NULL);

        return true;
    }

    template<typename T>
    bool add_mesh(generic_mesh<T,3>& msh, const std::string& name)
    {
        using mesh_type = generic_mesh<T,3>;

        static_assert(std::is_same<T,double>::value, "Only double for now");

        std::vector<T> x_coords, y_coords, z_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());
        z_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
            z_coords.push_back(pt.z());
        }

        T *coords[] = {x_coords.data(), y_coords.data(), z_coords.data()};

        int nfaces = msh.faces_size();

        std::vector<int> nodecnts;
        std::vector<int> nodelist;

        for (auto& fc : faces(msh))
        {
            auto ptids = fc.point_ids();
            nodecnts.push_back( ptids.size() );

            for (auto& ptid : ptids)
                nodelist.push_back(ptid);
        }
        int lnodelist = nodelist.size();
        int nzones = msh.cells_size();

        std::vector<int> facecnts, facelist;

        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);
            facecnts.push_back( fcs.size() );

            for (auto& fc : fcs)
            {
                int ofs = offset(msh,fc);
                if (silo_face_orientation(msh, cl, fc) < 0)
                    ofs = ~ofs;
                facelist.push_back( ofs );
            }
        };

        int lfacelist = facelist.size();

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        int err = DBPutPHZonelist(m_siloDb, zonelist_name.c_str(), nfaces, nodecnts.data(),
            lnodelist, nodelist.data(), nullptr, nzones, facecnts.data(), lfacelist,
            facelist.data(), 0, 0, nzones-1, nullptr);

        if (err != 0)
            std::cout << "DBPutPHZonelist() call failed" << std::endl;

        int ndims = 3;
        int nnodes = msh.points_size();

        DBoptlist *optlist = DBMakeOptlist(1);
        DBAddOption(optlist, DBOPT_PHZONELIST, (void*)zonelist_name.c_str());

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            NULL, NULL, DB_DOUBLE, optlist);

        DBFreeOptlist(optlist);

        return true;
    }

    template<typename T>
    bool add_mesh(const generic_mesh<T,2>& msh, const std::string& name)
    {
        static_assert(std::is_same<T,double>::value, "Only double for now");

        std::vector<T> x_coords, y_coords;
        x_coords.reserve(msh.points_size());
        y_coords.reserve(msh.points_size());

        for (auto itor = msh.points_begin(); itor != msh.points_end(); itor++)
        {
            auto pt = *itor;
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
        }

        T *coords[] = {x_coords.data(), y_coords.data()};

        std::vector<int>    nodelist;

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            auto size = ptids.size();
            assert(size > 2);
            for (auto& ptid : ptids)
                nodelist.push_back(ptid+1); // Silo wants 1-based indices
        }

        int lnodelist       = nodelist.size();
        int nshapetypes     = msh.cells_size();

        std::vector<int> shapesize, shapecount;
        for (auto& cl : msh)
        {
            auto fcs = faces(msh, cl);
            shapesize.push_back( fcs.size() );
            shapecount.push_back(1);
        };

        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 2;

        std::stringstream zlname;
        zlname << "zonelist_" << name;
        std::string zonelist_name = zlname.str();

        DBPutZonelist(m_siloDb, zonelist_name.c_str(), nzones, ndims,
            nodelist.data(), lnodelist, 1, shapesize.data(), shapecount.data(), nshapetypes);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            zonelist_name.c_str(), NULL, DB_DOUBLE, NULL);

        return true;
    }


    /* Scalar variable, REAL case */
    template<typename T, variable_centering_t centering>
    bool add_variable(const std::string& mesh_name,
                      silo_variable<T, centering>& var)
    {
        static_assert(std::is_same<T, double>::value, "Sorry, only double for now");

        if (!m_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }

        if (centering == zonal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var.name().c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_ZONECENT, NULL);
        }
        else if (centering == nodal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var.name().c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_NODECENT, NULL);
        }
        else
            return false;

        return true;
    }

    /* Scalar variable, COMPLEX case */
    template<typename T, variable_centering_t centering>
    bool add_variable(const std::string& mesh_name,
                      silo_variable<std::complex<T>, centering>& var)
    {
        static_assert(std::is_same<T, double>::value, "Sorry, only double for now");

        std::vector<T> var_rep(var.size());
        std::vector<T> var_imp(var.size());

        auto r = [](const std::complex<T>& z) -> T { return real(z); };
        auto i = [](const std::complex<T>& z) -> T { return imag(z); };

        std::transform(var.begin(), var.end(), var_rep.begin(), r);
        std::transform(var.begin(), var.end(), var_imp.begin(), i);

        if (!m_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }

        if (centering == zonal_variable_t)
        {
            std::stringstream ss;
            ss << var.name() << "_real";
            DBPutUcdvar1(m_siloDb, ss.str().c_str(), mesh_name.c_str(),
                         var_rep.data(),
                         var_rep.size(), NULL, 0, DB_DOUBLE,
                         DB_ZONECENT, NULL);
            ss.str(std::string());

            ss << var.name() << "_imag";
            DBPutUcdvar1(m_siloDb, ss.str().c_str(), mesh_name.c_str(),
                         var_imp.data(),
                         var_imp.size(), NULL, 0, DB_DOUBLE,
                         DB_ZONECENT, NULL);

        }
        else if (centering == nodal_variable_t)
        {
            std::stringstream ss;
            ss << var.name() << "_real";
            DBPutUcdvar1(m_siloDb, ss.str().c_str(), mesh_name.c_str(),
                         var_rep.data(),
                         var_rep.size(), NULL, 0, DB_DOUBLE,
                         DB_NODECENT, NULL);
            ss.str(std::string());

            ss << var.name() << "_imag";
            DBPutUcdvar1(m_siloDb, ss.str().c_str(), mesh_name.c_str(),
                         var_imp.data(),
                         var_imp.size(), NULL, 0, DB_DOUBLE,
                         DB_NODECENT, NULL);
        }
        else
            return false;

        return true;
    }

    /* Scalar variable, REAL case */
    template<typename T>
    bool add_variable(const std::string& mesh_name,
                      const std::string& var_name,
                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& var,
                      variable_centering_t centering)
    {
        static_assert(std::is_same<T, double>::value, "Sorry, only double for now");

        if (!m_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }

        if (centering == zonal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_ZONECENT, NULL);
        }
        else if (centering == nodal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_NODECENT, NULL);
        }
        else
            return false;

        return true;
    }

    template<typename T>
    bool add_variable(const std::string& mesh_name,
                      const std::string& var_name,
                      const Eigen::Matrix<T, Eigen::Dynamic, 2>& var,
                      variable_centering_t centering)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> c = var.col(0);
        add_variable(mesh_name, var_name+"_x", c, centering);
        c = var.col(1);
        add_variable(mesh_name, var_name+"_y", c, centering);
        std::stringstream exprss;
        exprss << "{" << var_name << "_x" << ", " << var_name << "_y" << "}";
        add_expression(var_name, exprss.str().c_str(), DB_VARTYPE_VECTOR);

        return true;
    }

    template<typename T>
    bool add_variable(const std::string& mesh_name,
                      const std::string& var_name,
                      const Eigen::Matrix<T, Eigen::Dynamic, 3>& var,
                      variable_centering_t centering)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> c = var.col(0);
        add_variable(mesh_name, var_name+"_x", c, centering);
        c = var.col(1);
        add_variable(mesh_name, var_name+"_y", c, centering);
        c = var.col(2);
        add_variable(mesh_name, var_name+"_z", c, centering);
        std::stringstream exprss;
        exprss << "{" << var_name << "_x" << ", " << var_name << "_y";
        exprss << ", " << var_name << "_z" << "}";
        add_expression(var_name, exprss.str().c_str(), DB_VARTYPE_VECTOR);

        return true;
    }

    /* Scalar variable, REAL case */
    template<typename T>
    bool add_variable(const std::string& mesh_name,
                      const std::string& var_name,
                      const std::vector<T>& var,
                      variable_centering_t centering)
    {
        static_assert(std::is_same<T, double>::value, "Sorry, only double for now");

        if (!m_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }

        if (centering == zonal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_ZONECENT, NULL);
        }
        else if (centering == nodal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         var.data(),
                         var.size(), NULL, 0, DB_DOUBLE,
                         DB_NODECENT, NULL);
        }
        else
            return false;

        return true;
    }

    bool add_variable(const std::string& mesh_name, const std::string& vars_name,
        std::vector<std::vector<double>>& vars)
    {
        std::vector<const double*> pvars;
        pvars.resize(vars.size());
        auto tr = [](const std::vector<double>& v) { return v.data(); };
        std::transform(vars.begin(), vars.end(), pvars.begin(), tr);

        int nvars = vars.size();
        int nels = vars[0].size();

        auto err = DBPutPointvar(m_siloDb, vars_name.c_str(), mesh_name.c_str(),
            nvars, pvars.data(), nels, DB_DOUBLE, nullptr);

        return err == 0;
    }

    bool add_variable(const std::string& mesh_name, const std::string& var_name,
        std::vector<double>& var)
    {
        int nels = var.size();

        auto err = DBPutPointvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
            var.data(), nels, DB_DOUBLE, nullptr);

        return err == 0;
    }

    bool add_expression(const std::string& expr_name,
                        const std::string& expr_definition,
                        int expr_type)
    {
        std::stringstream ss;
        ss << "def_" << expr_name;
        const char *name[] = { expr_name.c_str() };
        const char *def[] = { expr_definition.c_str() };
        DBPutDefvars(m_siloDb, ss.str().c_str(), 1, name, &expr_type, def, NULL);
        return true;
    }
};




} //namespaces
