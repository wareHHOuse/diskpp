option(BUILD_APP_LINEAR_ELASTICITY "Build linear elasticity application" ON)
if(BUILD_APP_LINEAR_ELASTICITY)
	add_executable(linear_elasticity linear_elasticity.cpp)
	target_link_libraries(linear_elasticity ${LINK_LIBS})

	add_executable(linear_elasticity_test linear_elasticity_test.cpp)
	target_link_libraries(linear_elasticity_test ${LINK_LIBS})
	install(TARGETS linear_elasticity_test RUNTIME DESTINATION bin)
endif()
