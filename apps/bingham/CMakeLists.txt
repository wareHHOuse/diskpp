option(BUILD_BINGHAM_APP "Build Bingham fluids application" ON)
if (BUILD_BINGHAM_APP)
add_executable(bingham_vector bingham_vector.cpp)
target_link_libraries(bingham_vector ${LINK_LIBS})

add_executable(bingham_vector_higher bingham_vector_higher.cpp)
target_link_libraries(bingham_vector_higher ${LINK_LIBS})

add_executable(test_norm test_norm.cpp)
target_link_libraries(test_norm ${LINK_LIBS})
endif(BUILD_BINGHAM_APP)

