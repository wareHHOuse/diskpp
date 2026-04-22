cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(linear_transmission linear_transmission.cpp)
target_link_libraries(linear_transmission diskpp::diskpp)
install(TARGETS linear_transmission RUNTIME DESTINATION bin)

add_executable(nonlinear_transmission nonlinear_transmission.cpp)
target_link_libraries(nonlinear_transmission diskpp::diskpp)

add_executable(test_projections test_projections.cpp)
target_link_libraries(test_projections diskpp::diskpp)