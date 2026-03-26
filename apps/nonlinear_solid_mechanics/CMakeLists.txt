cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(nonlinear_solid_mechanics src/nonlinear_solid_mechanics.cpp)
target_link_libraries(nonlinear_solid_mechanics diskpp::diskpp)
install(TARGETS nonlinear_solid_mechanics RUNTIME DESTINATION bin)

add_executable(contact_test src/contact_test.cpp)
target_link_libraries(contact_test diskpp::diskpp)
install(TARGETS contact_test RUNTIME DESTINATION bin)

add_executable(elasticity_test src/elasticity_test.cpp)
target_link_libraries(elasticity_test diskpp::diskpp)
install(TARGETS elasticity_test RUNTIME DESTINATION bin)

add_executable(elastodynamic_test src/elastodynamic_test.cpp)
target_link_libraries(elastodynamic_test diskpp::diskpp)
install(TARGETS elastodynamic_test RUNTIME DESTINATION bin)

add_executable(hyperelasticity_test src/hyperelasticity_test.cpp)
target_link_libraries(hyperelasticity_test diskpp::diskpp)
install(TARGETS hyperelasticity_test RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
