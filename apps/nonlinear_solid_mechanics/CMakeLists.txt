set(LINK_LIBS diskpp)

add_executable(nonlinear_solid_mechanics src/nonlinear_solid_mechanics.cpp)
target_link_libraries(nonlinear_solid_mechanics ${LINK_LIBS})
install(TARGETS nonlinear_solid_mechanics RUNTIME DESTINATION bin)

add_executable(contact_test src/contact_test.cpp)
target_link_libraries(contact_test ${LINK_LIBS})
install(TARGETS contact_test RUNTIME DESTINATION bin)

add_executable(elasticity_test src/elasticity_test.cpp)
target_link_libraries(elasticity_test ${LINK_LIBS})
install(TARGETS elasticity_test RUNTIME DESTINATION bin)

add_executable(hyperelasticity_test src/hyperelasticity_test.cpp)
target_link_libraries(hyperelasticity_test ${LINK_LIBS})
install(TARGETS hyperelasticity_test RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
