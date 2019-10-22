option(BUILD_MESH_BASIC_TESTS "Build mesh basic tests" OFF)
if (BUILD_MESH_BASIC_TESTS)
add_executable(mesh_basic_tests mesh_basic_tests.cpp)
target_link_libraries(mesh_basic_tests ${LINK_LIBS})

add_executable(mesh_cube_diff_test mesh_cube_diffusive_test.cpp)
target_link_libraries(mesh_cube_diff_test ${LINK_LIBS})

add_executable(netgen_2d_test netgen_2d_test.cpp)
target_link_libraries(netgen_2d_test ${LINK_LIBS})

add_executable(submesher_test submesher_test.cpp)
target_link_libraries(submesher_test ${LINK_LIBS})

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt4)
if (Qt4_FOUND)
    add_executable(fvca5_meshinfo fvca5_meshinfo.cpp MeshWidget.cpp)
    target_link_libraries(fvca5_meshinfo ${LINK_LIBS} Qt4::QtGui)
endif()
endif(BUILD_MESH_BASIC_TESTS)

