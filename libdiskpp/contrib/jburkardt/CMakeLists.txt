# Static library to collect code taken from
#   https://people.sc.fsu.edu/~jburkardt/

project(jburkardt)
cmake_minimum_required(VERSION 3.18)

set(jburkardt_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/src/triangle_dunavant_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/jb_common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simplex_gm_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tetrahedron_arbq_rule.cpp)

include_directories(include)
add_library(jburkardt_archive STATIC ${jburkardt_sources})

add_library(JBURKARDT INTERFACE)
target_link_libraries(JBURKARDT INTERFACE jburkardt_archive)
target_include_directories(JBURKARDT INTERFACE include)
