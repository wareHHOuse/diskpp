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

#include <cassert>

#include <hdf5.h>
#include <unistd.h>

#include "diskpp/common/eigen.hpp"

namespace disk{

class hdf5_context
{
    std::string     hdf5_filename;
    hid_t           file_id;
    herr_t          status;
    bool            is_open;

public:
    hdf5_context()
    {}

    hdf5_context(const std::string& filename)
        : hdf5_filename(filename), is_open(false)
    {
        file_id = H5Fcreate(hdf5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        is_open = true;
    }

    ~hdf5_context()
    {
        if (is_open)
            status = H5Fclose(file_id);
    }

    hid_t get_file_descriptor() const
    {
        assert(is_open);
        return file_id;
    }
};

hid_t
hdf5_open_or_create(const char *filename)
{
    hid_t file_id;

    if ( access(filename, R_OK | W_OK) == 0 )
        file_id = H5Fopen(filename,  H5F_ACC_RDWR, H5P_DEFAULT);
    else
        file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    return file_id;
}

} // end disk
