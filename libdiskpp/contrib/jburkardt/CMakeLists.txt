# Static library to collect code taken from
#   https://people.sc.fsu.edu/~jburkardt/

project(jburkardt)
cmake_minimum_required(VERSION 3.20)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/triangle_dunavant_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/jb_common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simplex_gm_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tetrahedron_arbq_rule.cpp)

add_library(jburkardt STATIC ${SOURCES})
set_target_properties(jburkardt PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(jburkardt
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

install(DIRECTORY include DESTINATION include)