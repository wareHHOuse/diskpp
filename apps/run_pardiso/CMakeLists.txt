cmake_minimum_required(VERSION 3.30)
project(run_pardiso)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(run_pardiso src/run_pardiso.cpp)
target_link_libraries(run_pardiso diskpp::diskpp)
install(TARGETS run_pardiso RUNTIME DESTINATION bin)
