cmake_minimum_required (VERSION 3.5)
project (diskpp)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

# project version
set(diskpp_VERSION_MAJOR 0)
set(diskpp_VERSION_MINOR 0)
set(diskpp_VERSION_PATCH 0)

# Additional modules path for cmake
set (CMAKE_MODULE_PATH
     ${CMAKE_MODULE_PATH}
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(safeguards)
include(agmg)

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -g -fpermissive")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -O3 -mavx -g -fpermissive -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASEASSERT "-std=c++14 -O3 -mavx -g -fpermissive")

find_package(MKL)
if (MKL_FOUND)
	include_directories(${MKL_INCLUDE_DIRS})
	set(LINK_LIBS ${LINK_LIBS} ${MKL_LIBRARIES})
	set(HAVE_INTEL_MKL TRUE)
endif()

find_package(SILO REQUIRED)
#if (SILO_FOUND)
    include_directories("${SILO_INCLUDE_DIRS}")
    set(LINK_LIBS ${LINK_LIBS} ${SILO_LIBRARIES})
    set(HAVE_SILO TRUE)
#endif()

find_package(HDF5 REQUIRED)
    include_directories(${HDF5_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
    set(HAVE_HDF5 TRUE)


find_package(Lua REQUIRED)
include_directories("${LUA_INCLUDE_DIR}")
set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})

set(UP_EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen")
if ( EXISTS ${UP_EIGEN_DIR} AND IS_DIRECTORY ${UP_EIGEN_DIR} )
    include_directories(${UP_EIGEN_DIR})
else ( EXISTS ${UP_EIGEN_DIR} AND IS_DIRECTORY ${UP_EIGEN_DIR} )
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ( EXISTS ${UP_EIGEN_DIR} AND IS_DIRECTORY ${UP_EIGEN_DIR} )


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")

include_directories(core contrib)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(LINK_LIBS "${LINK_LIBS}" pthread)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
    set(LINK_LIBS "${LINK_LIBS}" pthread)
endif()

add_subdirectory(contrib)
add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(doc)

add_library(diskpp SHARED ${diskpp_source})
set(LINK_LIBS "${LINK_LIBS}" diskpp)
install(TARGETS diskpp LIBRARY DESTINATION lib)

install(DIRECTORY meshes DESTINATION share)

add_subdirectory(apps)
