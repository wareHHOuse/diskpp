cmake_minimum_required (VERSION 2.8)
project (diskpp)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

#set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fpermissive -Wall -Wshadow -Wno-unused-variable -Wno-unused-local-typedef")
#set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -O3 -mavx -g -fpermissive -DNDEBUG -Wall -Wshadow")
#set(CMAKE_CXX_FLAGS_RELEASEASSERT "-std=c++14 -O3 -mavx -g -fpermissive -Wall -Wshadow")

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fpermissive")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -O3 -mavx -g -fpermissive -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASEASSERT "-std=c++14 -O3 -mavx -g -fpermissive")

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

set (SOLVER_WRAPPERS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/contrib/solver_wrappers/")

if ( EXISTS "${SOLVER_WRAPPERS_SOURCE_DIR}/CMakeLists.txt" )
    add_subdirectory(${SOLVER_WRAPPERS_SOURCE_DIR})
    include_directories("${SOLVER_WRAPPERS_SOURCE_DIR}")
    set(LINK_LIBS "${LINK_LIBS}" ${SOLVER_WRAPPERS_LIBS})
    set (HAVE_SOLVER_WRAPPERS TRUE)
endif ()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/config.h"
  )
include_directories(core contrib "${PROJECT_SOURCE_DIR}")

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(LINK_LIBS "${LINK_LIBS}" pthread)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
    set(LINK_LIBS "${LINK_LIBS}" pthread)
endif()

add_subdirectory(contrib)
add_subdirectory(core)

add_library(diskpp SHARED ${diskpp_source})
set(LINK_LIBS "${LINK_LIBS}" diskpp)

add_subdirectory(apps)
