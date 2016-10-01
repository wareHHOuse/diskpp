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

if (CYGWIN)
	SET(_CFLP "${CMAKE_FIND_LIBRARY_PREFIXES}")
	SET(_CFLS "${CMAKE_FIND_LIBRARY_SUFFIXES}")
	SET(CMAKE_FIND_LIBRARY_PREFIXES "")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
endif ()


if ( CYGWIN )
find_library(MKL_iomp5_LIBRARY
    NAMES   libiomp5md
    PATHS	"/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/redist/intel64/compiler/"
            "${INTEL_MKL_LIB_SEARCH_DIRS}"
)
else ()
find_library(MKL_iomp5_LIBRARY
    NAMES   iomp5
    PATHS	/opt/intel/lib
    			/opt/intel/lib/intel64
    		    /opt/intel/lib/intel64_lin
            "${INTEL_MKL_LIB_SEARCH_DIRS}"
)
endif()

find_library(MKL_mkl_core_LIBRARY
    NAMES   mkl_core
    PATHS   /opt/intel/mkl/lib
	    /opt/intel/mkl/lib/intel64
	    "/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/redist/intel64/mkl/"
            "${INTEL_MKL_LIB_SEARCH_DIRS}"
)
find_library(MKL_mkl_intel_thread_LIBRARY
    NAMES   mkl_intel_thread
    PATHS   /opt/intel/mkl/lib
            /opt/intel/mkl/lib/intel64
	    "/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/redist/intel64/mkl/"
            "${INTEL_MKL_LIB_SEARCH_DIRS}"
)

if (APPLE)
	find_library(MKL_runtime_LIBRARY
		NAMES   mkl_intel
		PATHS   /opt/intel/mkl/lib
			"/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/redist/intel64/mkl/"
            "${INTEL_MKL_LIB_SEARCH_DIRS}"
)
else ()
	find_library(MKL_runtime_LIBRARY
		NAMES   mkl_rt
		PATHS   /opt/intel/mkl/lib			
	    		/opt/intel/mkl/lib/intel64
			"/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/redist/intel64/mkl/"
			"${INTEL_MKL_LIB_SEARCH_DIRS}"
	)
endif ()

find_path(MKL_INCLUDE_DIRS
    NAMES   mkl.h
    PATHS   /opt/intel/mkl/include
	    "/cygdrive/C/Program\ Files\ \(x86\)/Intel/Composer\ XE/mkl/include"
            "${INTEL_MKL_INCLUDE_SEARCH_DIRS}/include"
)


if (CYGWIN)
	SET(CMAKE_FIND_LIBRARY_PREFIXES "${_CFLP}")
	SET(CMAKE_FIND_LIBRARY_SUFFIXES "${_CFLS}")
endif()


if (MKL_mkl_core_LIBRARY AND
	MKL_runtime_LIBRARY AND
	MKL_mkl_intel_thread_LIBRARY AND
	MKL_iomp5_LIBRARY AND
	MKL_INCLUDE_DIRS)
	set(MKL_FOUND TRUE)
	set(MKL_LIBRARIES ${MKL_mkl_core_LIBRARY} 
					  ${MKL_runtime_LIBRARY}
					  ${MKL_mkl_intel_thread_LIBRARY} ${MKL_iomp5_LIBRARY})
	set (HAVE_PARDISO TRUE)
endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS)


