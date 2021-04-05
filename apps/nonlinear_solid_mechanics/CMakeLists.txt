option(BUILD_APP_NSM "Build Non-linear Solid Mechanics applications" ON)
if (BUILD_APP_NSM)
	add_executable(nonlinear_solid_mechanics nonlinear_solid_mechanics.cpp)
	target_link_libraries(nonlinear_solid_mechanics ${LINK_LIBS})
	install(TARGETS nonlinear_solid_mechanics RUNTIME DESTINATION bin)

	add_executable(contact_test contact_test.cpp)
	target_link_libraries(contact_test ${LINK_LIBS})
	install(TARGETS contact_test RUNTIME DESTINATION bin)

	add_executable(elasticity_test elasticity_test.cpp)
	target_link_libraries(elasticity_test ${LINK_LIBS})
	install(TARGETS elasticity_test RUNTIME DESTINATION bin)

	add_executable(hyperelasticity_test hyperelasticity_test.cpp)
	target_link_libraries(hyperelasticity_test ${LINK_LIBS})
	install(TARGETS hyperelasticity_test RUNTIME DESTINATION bin)

endif ()