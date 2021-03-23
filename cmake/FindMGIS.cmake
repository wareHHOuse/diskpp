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

find_library(MGIS_MFRONT_LIBRARY
    NAMES   libMFrontGenericInterface.so libMFrontGenericInterface.dylib
    PATHS	/usr/local/lib /usr/lib
)

#if(MGIS_MFRONT_LIBRARY)
#    message("MGIS_MFRONT_LIBRARY found")
#else()
#    message("MGIS_MFRONT_LIBRARY not found")
#endif()


find_path(MGIS_INCLUDE_DIRS
    NAMES   Integrate.hxx
    PATHS   /usr/local/include/MGIS/Behaviour /usr/include/MGIS/Behaviour
)

#if(MGIS_INCLUDE_DIRS)
#    message("MGIS_INCLUDE_DIRS found")
#else()
#    message("MGIS_INCLUDE_DIRS not found")
#endif()



if (MGIS_MFRONT_LIBRARY AND MGIS_INCLUDE_DIRS)
	set(MGIS_FOUND TRUE)
	set(MGIS_LIBRARIES ${MGIS_MFRONT_LIBRARY} )
endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MGIS DEFAULT_MSG MGIS_LIBRARIES MGIS_INCLUDE_DIRS)


