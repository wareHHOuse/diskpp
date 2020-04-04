option(BUILD_APP_NEWCODE "Build newcode application" ON)
if (BUILD_APP_NEWCODE)
	add_executable(new_code_test new_code_test.cpp)
	target_link_libraries(new_code_test ${LINK_LIBS})

	add_executable(test_bases test_bases.cpp)
	target_link_libraries(test_bases ${LINK_LIBS})
endif()