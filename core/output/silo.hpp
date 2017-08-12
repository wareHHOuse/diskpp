/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016, 2017 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
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


#include <silo.h>


namespace disk {
    
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
    bool add_mesh(const simplicial_mesh<T,2>& msh, const std::string& name)
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
        nodelist.reserve( 3*msh.cells_size() );

        for (auto& cl : msh)
        {
            auto ptids = cl.point_ids();
            assert(ptids.size() == 3);

            for (auto& ptid : ptids)
                nodelist.push_back( ptid + 1 );
        }

        int lnodelist = nodelist.size();

        int shapesize[] = {3};
        int shapecounts[] = { static_cast<int>(msh.cells_size()) };
        int nshapetypes = 1;
        int nnodes = msh.points_size();
        int nzones = msh.cells_size();
        int ndims = 2;

        DBPutZonelist(m_siloDb, "zonelist", nzones, ndims, nodelist.data(), lnodelist,
            1, shapesize, shapecounts, nshapetypes);
        
        DBPutUcdmesh(m_siloDb, "mesh", ndims, NULL, coords, nnodes, nzones,
            "zonelist", NULL, DB_DOUBLE, NULL);

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

        DBPutZonelist(m_siloDb, "zonelist", nzones, ndims, nodelist.data(), lnodelist,
            1, shapesize, shapecounts, nshapetypes);
        
        DBPutUcdmesh(m_siloDb, "mesh", ndims, NULL, coords, nnodes, nzones,
            "zonelist", NULL, DB_DOUBLE, NULL);
        
        return true;
    }

    template<typename T>
    bool add_mesh(const generic_mesh<T,2>& msh, const std::string& name)
    {
        std::cout << "Export of generic meshes not yet supported." << std::endl;
        return false;
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
};
    
    


} //namespaces













