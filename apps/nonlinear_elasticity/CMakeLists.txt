add_executable(nonlinear_elasticity nonlinear_elasticity.cpp)
target_link_libraries(nonlinear_elasticity ${LINK_LIBS})
install(TARGETS nonlinear_elasticity RUNTIME DESTINATION bin)


add_executable(nonlinear_elasticity_tests nonlinear_elasticity_tests.cpp)
target_link_libraries(nonlinear_elasticity_tests ${LINK_LIBS})
install(TARGETS nonlinear_elasticity_tests RUNTIME DESTINATION bin)
