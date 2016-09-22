add_executable(mesh_basic_tests mesh_basic_tests.cpp)
target_link_libraries(mesh_basic_tests ${LINK_LIBS})

add_executable(netgen_2d_test netgen_2d_test.cpp)
target_link_libraries(netgen_2d_test ${LINK_LIBS})

add_executable(submesher_test submesher_test.cpp)
target_link_libraries(submesher_test ${LINK_LIBS})
