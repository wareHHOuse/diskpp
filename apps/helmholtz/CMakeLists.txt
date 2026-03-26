cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(helmholtz helmholtz.cpp)
target_link_libraries(helmholtz diskpp::diskpp)
install(TARGETS helmholtz RUNTIME DESTINATION bin)