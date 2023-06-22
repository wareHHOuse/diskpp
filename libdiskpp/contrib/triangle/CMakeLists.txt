

add_library(triangle_archive STATIC ${CMAKE_CURRENT_SOURCE_DIR}/triangle_mesher.c)
target_compile_definitions(triangle_archive PUBLIC -DTRILIBRARY -DREDUCED -DCDT_ONLY)

add_library(TRIANGLE INTERFACE)
target_link_libraries(TRIANGLE INTERFACE triangle_archive)
target_include_directories(TRIANGLE INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

