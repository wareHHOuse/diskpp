#
#       /\        Matteo Cicuttin (C) 2016-2021
#      /__\       matteo.cicuttin@enpc.fr
#     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
#    /\    /\
#   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
#  /_\/_\/_\/_\   methods.
#
# This file is copyright of the following authors:
# Matteo Cicuttin (C) 2016, 2017, 2018         matteo.cicuttin@enpc.fr
# Nicolas Pignet  (C) 2021                     nicolas.pignet@enpc.fr
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# If you use this code or parts of it for scientific publications, you
# are required to cite it as following:
#
# Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
# polytopal meshes using generic programming.
# M. Cicuttin, D. A. Di Pietro, A. Ern.
# Journal of Computational and Applied Mathematics.
# DOI: 10.1016/j.cam.2017.09.017
#

include(FindPackageHandleStandardArgs)

find_path(MGIS_INCLUDE_DIRS
    NAMES MGIS/Behaviour/Integrate.hxx
    HINTS ENV MGIS_ROOT ${MGIS_ROOT}
    PATH_SUFFIXES include
)

find_library(MGIS_MFRONT_LIBRARIES
    NAMES MFrontGenericInterface
    HINTS ENV MGIS_ROOT ${MGIS_ROOT}
    PATH_SUFFIXES lib
)

find_package_handle_standard_args(MGIS DEFAULT_MSG MGIS_MFRONT_LIBRARIES MGIS_INCLUDE_DIRS)

if (MGIS_FOUND)
    add_library(MGIS INTERFACE)
    target_link_libraries(MGIS INTERFACE "${MGIS_MFRONT_LIBRARIES}")
    target_include_directories(MGIS INTERFACE "${MGIS_INCLUDE_DIRS}")
endif()
