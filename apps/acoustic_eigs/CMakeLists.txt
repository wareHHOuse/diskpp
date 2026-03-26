cmake_minimum_required(VERSION 3.30)
project(acoustic_eigs)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(acoustic_eigs src/acoustic_eigs.cpp)
target_link_libraries(acoustic_eigs diskpp::diskpp)
install(TARGETS acoustic_eigs RUNTIME DESTINATION bin)

