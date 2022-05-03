#
# This source file is part of EMT, the ElectroMagneticTool.
#
# Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
# Department of Electrical Engineering, University of Udine
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Udine nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

include(FindPackageHandleStandardArgs)

#set(MKL_HINTS /usr /opt/intel /opt/intel/mkl /opt/intel/oneapi/mkl)

find_path(MKL_INCLUDE_DIR
    NAMES mkl.h
    HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
    PATH_SUFFIXES mkl mkl/include)

###########################################################
## iomp5
find_library(MKL_iomp5_LIBRARY
    NAMES iomp5
    HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
    PATH_SUFFIXES lib lib/intel64 lib/intel64_lin)

if (MKL_iomp5_LIBRARY)
    message(STATUS "MKL: found iomp5")
    list(APPEND MKL_LIBRARIES ${MKL_iomp5_LIBRARY})
endif()

###########################################################
## mkl_core 
find_library(MKL_mkl_core_LIBRARY
    NAMES   mkl_core
    HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
    PATH_SUFFIXES lib lib/intel64 lib/intel64_lin
                  mkl/lib mkl/lib/intel64 mkl/lib/intel64_lin)

if (MKL_mkl_core_LIBRARY)
    message(STATUS "MKL: found mkl_core")
    list(APPEND MKL_LIBRARIES ${MKL_mkl_core_LIBRARY})
endif()

###########################################################
## mkl_intel_thread 
find_library(MKL_mkl_intel_thread_LIBRARY
    NAMES   mkl_intel_thread
    HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
    PATH_SUFFIXES lib lib/intel64 lib/intel64_lin
                  mkl/lib mkl/lib/intel64 mkl/lib/intel64_lin)

if (MKL_mkl_intel_thread_LIBRARY)
    message(STATUS "MKL: found mkl_intel_thread")
    list(APPEND MKL_LIBRARIES ${MKL_mkl_intel_thread_LIBRARY})
endif()

###########################################################
## mkl_runtime 
if (APPLE)
    find_library(MKL_runtime_LIBRARY
        NAMES   mkl_intel
        HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
        PATH_SUFFIXES lib lib/intel64 lib/intel64_lin
                      mkl/lib mkl/lib/intel64 mkl/lib/intel64_lin)
else ()
    find_library(MKL_runtime_LIBRARY
        NAMES   mkl_rt
        HINTS ENV MKL_ROOT ${MKL_ROOT} ${MKL_HINTS}
        PATH_SUFFIXES lib lib/intel64 lib/intel64_lin
                      mkl/lib mkl/lib/intel64 mkl/lib/intel64_lin)
endif ()

if (MKL_runtime_LIBRARY)
    message(STATUS "MKL: found mkl_runtime")
    list(APPEND MKL_LIBRARIES ${MKL_runtime_LIBRARY})
endif()


find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIR)

if (MKL_FOUND)
    add_library(MKL INTERFACE)
    target_link_libraries(MKL INTERFACE "${MKL_LIBRARIES}")
    target_include_directories(MKL INTERFACE "${MKL_INCLUDE_DIR}")
    target_compile_definitions(MKL INTERFACE -DHAVE_INTEL_MKL)
    target_compile_definitions(MKL INTERFACE -DHAVE_PARDISO)
endif()

