cmake_minimum_required(VERSION 3.30)
project(barycenter_test)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(barycenter_test barycenter_test.cpp)
target_link_libraries(barycenter_test diskpp::diskpp)
install(TARGETS barycenter_test RUNTIME DESTINATION bin)