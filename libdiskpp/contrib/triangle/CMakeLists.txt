project(triangle)
cmake_minimum_required(VERSION 3.20)


add_library(triangle STATIC triangle_mesher.c)
set_target_properties(triangle PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(triangle PUBLIC -DTRILIBRARY -DREDUCED -DCDT_ONLY)


target_include_directories(triangle
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
        $<INSTALL_INTERFACE:include>
)

install(FILES triangle_mesher.h DESTINATION include)