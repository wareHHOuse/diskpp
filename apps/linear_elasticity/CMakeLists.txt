cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(linear_elasticity src/linear_elasticity.cpp)
target_link_libraries(linear_elasticity diskpp::diskpp)
install(TARGETS linear_elasticity RUNTIME DESTINATION bin)

add_executable(linear_elasticity_test src/linear_elasticity_test.cpp)
target_link_libraries(linear_elasticity_test diskpp::diskpp)
install(TARGETS linear_elasticity_test RUNTIME DESTINATION bin)

#install(DIRECTORY share/ DESTINATION share)
