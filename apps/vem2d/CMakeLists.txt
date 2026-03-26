cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(test_local_matrices src/test_local_matrices.cpp)
target_link_libraries(test_local_matrices diskpp::diskpp)
install(TARGETS test_local_matrices RUNTIME DESTINATION bin)

