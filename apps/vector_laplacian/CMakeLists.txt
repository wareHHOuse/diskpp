option(BUILD_APP_VECTOR_LAPLACIAN "Build vector laplacian application" ON)
if (BUILD_APP_VECTOR_LAPLACIAN)
	add_executable(vector_laplacian vector_laplacian.cpp)
	target_link_libraries(vector_laplacian ${LINK_LIBS})

	add_executable(vector_laplacian_test vector_laplacian_test.cpp)
	target_link_libraries(vector_laplacian_test ${LINK_LIBS})
	install(TARGETS vector_laplacian_test RUNTIME DESTINATION bin)
endif()