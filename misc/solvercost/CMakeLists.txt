cmake_minimum_required (VERSION 3.0)
project (diskpp)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

enable_language(Fortran)

set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fpermissive")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -O3 -mavx -g -fpermissive -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASEASSERT "-std=c++14 -O3 -mavx -g -fpermissive")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(MUMPS REQUIRED)
include_directories(${MUMPS_INCLUDE_DIR})
set(LINK_LIBS ${LINK_LIBS} ${MUMPS_LIBRARIES})

find_package(BLAS REQUIRED)
set(LINK_LIBS ${LINK_LIBS} ${BLAS_LIBRARIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES})

add_executable(solvercost solvercost.cpp)
target_link_libraries(solvercost ${LINK_LIBS})
