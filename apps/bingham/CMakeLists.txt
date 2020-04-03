option(BUILD_APP_BINGHAM "Build Bingham fluids applications" ON)
if (BUILD_APP_BINGHAM)
	add_executable(bingham_vector bingham_vector.cpp)
	target_link_libraries(bingham_vector ${LINK_LIBS})

	add_executable(bingham_vector_higher bingham_vector_higher.cpp)
	target_link_libraries(bingham_vector_higher ${LINK_LIBS})

	add_executable(test_norm test_norm.cpp)
	target_link_libraries(test_norm ${LINK_LIBS})
endif()

