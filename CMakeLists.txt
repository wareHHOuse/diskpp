cmake_minimum_required (VERSION 3.0)
project (diskpp)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fpermissive -Wall -Wshadow -Wno-unused-variable -Wno-unused-local-typedef")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -O3 -mavx -g -fpermissive -DNDEBUG -Wall -Wshadow")
set(CMAKE_CXX_FLAGS_RELEASEASSERT "-std=c++14 -O3 -mavx -g -fpermissive -Wall -Wshadow")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# Additional modules path for cmake
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(MKL)
if (MKL_FOUND)
	include_directories(${MKL_INCLUDE_DIRS})
	set(LINK_LIBS "${LINK_LIBS}" "${MKL_LIBRARIES}")
	set(HAVE_INTEL_MKL TRUE)
endif()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/_config_/config.h"
  )

include_directories(core contrib "${PROJECT_BINARY_DIR}/_config_/")

add_subdirectory(contrib)
add_subdirectory(core)

add_library(diskpp SHARED ${diskpp_source})
set(LINK_LIBS "${LINK_LIBS}" diskpp)

add_subdirectory(apps)
