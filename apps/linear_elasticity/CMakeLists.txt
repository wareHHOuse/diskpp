set(LINK_LIBS diskpp)

add_executable(linear_elasticity src/linear_elasticity.cpp)
target_link_libraries(linear_elasticity ${LINK_LIBS})
install(TARGETS linear_elasticity RUNTIME DESTINATION bin)

add_executable(linear_elasticity_test src/linear_elasticity_test.cpp)
target_link_libraries(linear_elasticity_test ${LINK_LIBS})
install(TARGETS linear_elasticity_test RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
