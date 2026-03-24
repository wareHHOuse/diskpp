cmake_minimum_required(VERSION 3.18)
project(libdiskpp VERSION 0.0.0)

include(CheckLanguage)

if ( WIN32 )
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

######################################################################
## Additional modules path for cmake
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_link_options("-Wl,--no-as-needed")
endif()

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

######################################################################
## RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

######################################################################
## Find Eigen 3
find_package(Eigen3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

######################################################################
## MUMPS solver
find_package(MUMPS)
set(LINK_LIBS ${LINK_LIBS} MUMPS::MUMPS)

######################################################################
## Package: Intel MKL support
set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL CONFIG)
if (MKL_FOUND)
  set(LINK_LIBS ${LINK_LIBS} MKL::MKL)
  set(HAVE_INTEL_MKL TRUE)
  add_compile_definitions(-DHAVE_MKL)
  add_compile_definitions(-DHAVE_INTEL_MKL)
  add_compile_definitions(-DHAVE_PARDISO)
endif()

######################################################################
## Package: LLNL Silo library support
find_package(SILO REQUIRED)
if (SILO_FOUND)
    set(LINK_LIBS ${LINK_LIBS} SILO::SILO)
endif()

######################################################################
## Package: Lua language support
find_package(Lua REQUIRED)
include_directories("${LUA_INCLUDE_DIR}")
set(LINK_LIBS ${LINK_LIBS} ${LUA_LIBRARIES})

######################################################################
## Package: GMSH
find_package(GMSH)
if (GMSH_FOUND)
    set(HAVE_GMSH TRUE)
    set(LINK_LIBS ${LINK_LIBS} GMSH::GMSH)
endif()

######################################################################
## Now DiSk++ sources
set(LIBDISKPP_SOURCES
    src/common/util.cpp src/common/mapped_file.cpp src/common/colormanip.cpp
)

######################################################################
## And the static libraries
add_library(private_stuff INTERFACE)

add_subdirectory(contrib/jburkardt)
target_link_libraries(private_stuff PRIVATE INTERFACE JBURKARDT)

add_subdirectory(contrib/triangle)
target_link_libraries(private_stuff PRIVATE INTERFACE TRIANGLE)

######################################################################
##


set(SOL2_BUILD_LUA OFF CACHE BOOL "Override sol2 build lua option")
add_subdirectory(contrib/sol2)
set(LINK_LIBS ${LINK_LIBS} sol2::sol2)

add_library(diskpp SHARED ${LIBDISKPP_SOURCES})
target_link_libraries(diskpp PUBLIC ${LINK_LIBS})
target_link_libraries(diskpp PRIVATE private_stuff)
set_target_properties(diskpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(diskpp PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

add_library(diskpp::diskpp ALIAS diskpp)




install(TARGETS diskpp
    EXPORT DiskppTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/cmake/"
        DESTINATION lib/cmake/diskpp)

install(DIRECTORY include/ DESTINATION include)

install(EXPORT DiskppTargets
    FILE DiskppTargets.cmake
    NAMESPACE diskpp::
    DESTINATION lib/cmake/diskpp
)

export(TARGETS diskpp
       NAMESPACE diskpp::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/DiskppTargets.cmake")

export(PACKAGE sol2)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DiskppConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DiskppConfigVersion.cmake"
    DESTINATION lib/cmake/diskpp
)


include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DiskppConfigVersion.cmake"
    VERSION 0.0.0
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DiskppConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DiskppConfig.cmake"
    INSTALL_DESTINATION lib/cmake/diskpp
)
