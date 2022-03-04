cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(diskpp)

include(CheckLanguage)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

######################################################################
## project version
set(diskpp_VERSION_MAJOR 0)
set(diskpp_VERSION_MINOR 0)
set(diskpp_VERSION_PATCH 0)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

######################################################################
## RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

######################################################################
## Build Release by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

######################################################################
## Disallow in-tree builds
if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please create a build directory and run CMake from there.")
endif()

######################################################################
## Additional modules path for cmake
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set (SITE_CONFIGS "${PROJECT_SOURCE_DIR}/cmake/sites")

set(SITE_PRESET "default" CACHE STRING "Site-specific environment to load")
if ( NOT (${SITE_PRESET} STREQUAL "default") )
    include("${SITE_CONFIGS}/Site${SITE_PRESET}.cmake")
endif()

######################################################################
## Update submodules if necessary
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

######################################################################
## Eigen3
add_subdirectory(contrib/eigen)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

######################################################################
## Solver wrappers
add_subdirectory(contrib/solver_wrappers)
set(LINK_LIBS ${LINK_LIBS} SOLVER_WRAPPERS)

######################################################################
## Helpers for compiler identification
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(COMPILER_IS_CLANG TRUE)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(COMPILER_IS_GNU TRUE)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(COMPILER_IS_INTEL)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(COMPILER_IS_PGI)
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    cmake_host_system_information(RESULT OSVER QUERY OS_RELEASE)
    if(OSVER VERSION_LESS "10.14")
        message(FATAL_ERROR "DiSk++ depends on libraries that use C++17 \
                 features available starting with Mac OS X 10.14. You should \
                 either upgrade your system or use a compiler different than \
                 AppleClang.")
    endif()
endif()

######################################################################
## Package: Threading support
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if (COMPILER_IS_CLANG OR COMPILER_IS_INTEL OR COMPILER_IS_GNU)
    set(LINK_LIBS ${LINK_LIBS} Threads::Threads)
endif()

######################################################################
## Package: Intel MKL support - To be made optional
find_package(MKL REQUIRED)
if (MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
  set(LINK_LIBS ${LINK_LIBS} ${MKL_LIBRARIES})
  set(HAVE_INTEL_MKL TRUE)
  add_definitions(-DHAVE_INTEL_MKL)
endif()

######################################################################
## Package: LLNL Silo library support - To be made optional
find_package(SILO REQUIRED)
if (SILO_FOUND)
    set(LINK_LIBS ${LINK_LIBS} SILO)
endif()

######################################################################
## Package: HDF5 Support - To be made optional
find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS})
set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
set(HAVE_HDF5 TRUE)

######################################################################
## Package: Lua language support
find_package(Lua REQUIRED)
include_directories("${LUA_INCLUDE_DIR}")
set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})

######################################################################
## Package: MGIS support
find_package(MGIS)
if (MGIS_FOUND)
    set(LINK_LIBS ${LINK_LIBS} MGIS)
    set(HAVE_MGIS TRUE)
    add_definitions(-DHAVE_MGIS)
endif()

######################################################################
## Package: GMSH
find_package(GMSH)
if (GMSH_FOUND)
    set(LINK_LIBS ${LINK_LIBS} GMSH)
endif()

######################################################################
## Optimization: Prefer 512 bit vectors on certain architectures.
## See https://reviews.llvm.org/D67259 for more info
option(OPT_PREFER_512bit "Prefer 512 bit vectors with AVX512 (Clang > 10 & GCC)" OFF)
if (OPT_PREFER_512bit)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mprefer-vector-width=512")
endif()

######################################################################
## Optimization: Disallow denormals and activate -ffast-math (violates IEEE754)
option(OPT_AGGRESSIVE_FP "Enable DAZ, FTZ and -ffast-math" ON)
if (OPT_AGGRESSIVE_FP)
    add_definitions(-DDISALLOW_DENORMALS)
    if (COMPILER_IS_CLANG OR COMPILER_IS_GNU)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
    endif()
endif()

######################################################################
## Optimization: Enable vectorizer output
option(OPT_VECTORIZER_REMARKS "Enable vectorizer remarks" OFF)
if (OPT_VECTORIZER_REMARKS)
    if (COMPILER_IS_CLANG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-missed=loop-vectorize")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize")
    endif()

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report-phase=vec -qopt-report=2")
    endif()

    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized")
    endif()
endif()

######################################################################
## Compiler settings: Settings for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address -fpermissive")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -g -fpermissive -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASEASSERT "-O3 -march=native -g -fpermissive")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-anon-enum-enum-conversion")

include_directories("${PROJECT_SOURCE_DIR}")

include_directories(core contrib)

SET(SOL2_BUILD_LUA OFF CACHE BOOL "Override sol2 build lua option")
add_subdirectory(contrib/sol2)
include_directories(contrib/sol2/include)
include_directories(contrib/sol2/single)

add_subdirectory(contrib)
add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(doc)

add_library(diskpp SHARED ${diskpp_source})
set(LINK_LIBS "${LINK_LIBS}" diskpp)
install(TARGETS diskpp LIBRARY DESTINATION lib)

install(DIRECTORY meshes DESTINATION share)

add_subdirectory(apps)
