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
## Now the deps
add_library(diskpp_deps INTERFACE)

######################################################################
## Now DiSk++ sources
set(LIBDISKPP_SOURCES
    src/common/util.cpp src/common/mapped_file.cpp src/common/colormanip.cpp
)
add_library(diskpp SHARED ${LIBDISKPP_SOURCES})

######################################################################
## Find Eigen 3
find_package(Eigen3 REQUIRED)
target_link_libraries(diskpp_deps INTERFACE Eigen3::Eigen)

######################################################################
## MUMPS solver
find_package(MUMPS)
if (MUMPS_FOUND)
    target_link_libraries(diskpp_deps INTERFACE MUMPS::MUMPS)
endif()

######################################################################
## Package: Intel MKL support
set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL CONFIG)
if (MKL_FOUND)
    target_link_libraries(diskpp_deps INTERFACE MKL::MKL)
    target_compile_definitions(diskpp_deps INTERFACE
        HAVE_MKL
        HAVE_INTEL_MKL
        HAVE_PARDISO
    )
endif()

######################################################################
## Package: LLNL Silo library support
find_package(SILO REQUIRED)
target_link_libraries(diskpp_deps INTERFACE SILO::SILO)


######################################################################
## Package: Lua language support
find_package(Lua REQUIRED)
if(Lua_FOUND AND NOT TARGET Lua::Lua)
  add_library(Lua::Lua INTERFACE IMPORTED)
  set_target_properties(
    Lua::Lua
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${LUA_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${LUA_LIBRARIES}"
  )
  target_link_libraries(diskpp_deps INTERFACE Lua::Lua)
endif()

######################################################################
## Package: GMSH
find_package(GMSH)
if (GMSH_FOUND)
    target_link_libraries(diskpp_deps INTERFACE GMSH::GMSH)
    target_compile_definitions(diskpp_deps INTERFACE
        HAVE_GMSH
    )
endif()

######################################################################
## Sol2

set(SOL2_BUILD_LUA OFF CACHE BOOL "Override sol2 build lua option")
add_subdirectory(contrib/sol2)
target_link_libraries(diskpp_deps INTERFACE sol2::sol2)

######################################################################
## HDF5 Support
find_package(HDF5)
if (HDF5_FOUND)
    add_subdirectory(contrib/highfive)
    target_link_libraries(diskpp_deps INTERFACE HighFive)
    target_compile_definitions(diskpp_deps INTERFACE 
        HAVE_HDF5
    )
endif()

######################################################################
## Triangle
add_subdirectory(contrib/triangle)
target_link_libraries(diskpp PUBLIC triangle)
target_include_directories(diskpp
    PUBLIC
        $<TARGET_PROPERTY:triangle,INTERFACE_INCLUDE_DIRECTORIES>
)

######################################################################
## JBurkardt stuff
add_subdirectory(contrib/jburkardt)
target_link_libraries(diskpp PUBLIC jburkardt)
target_include_directories(diskpp
    PUBLIC
        $<TARGET_PROPERTY:jburkardt,INTERFACE_INCLUDE_DIRECTORIES>
)

######################################################################
## gmsh_tools
add_subdirectory(contrib/gmsh_tools)
target_link_libraries(diskpp_deps INTERFACE gmshtools)
target_include_directories(diskpp
    PUBLIC
        $<TARGET_PROPERTY:gmshtools,INTERFACE_INCLUDE_DIRECTORIES>
)

######################################################################
## sgr
add_subdirectory(contrib/sgr)
target_link_libraries(diskpp INTERFACE sgr)
target_include_directories(diskpp
    PUBLIC
        $<TARGET_PROPERTY:sgr,INTERFACE_INCLUDE_DIRECTORIES>
)


target_link_libraries(diskpp PUBLIC diskpp_deps)

set_target_properties(diskpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(diskpp PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

add_library(diskpp::diskpp ALIAS diskpp)


install(TARGETS diskpp triangle jburkardt diskpp_deps gmshtools sgr sol2
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

export(TARGETS diskpp triangle jburkardt diskpp_deps gmshtools sgr sol2
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
