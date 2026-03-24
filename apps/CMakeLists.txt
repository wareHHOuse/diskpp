cmake_minimum_required(VERSION 3.30)
project(diskpp_apps)

set(Diskpp_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../libdiskpp/build")

option(BUILD_APP_ACOUSTIC_EIGS "Build acoustic_eigs application" ON)
if (BUILD_APP_ACOUSTIC_EIGS)
    add_subdirectory(acoustic_eigs)
endif()