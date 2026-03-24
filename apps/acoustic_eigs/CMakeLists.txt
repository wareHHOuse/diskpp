cmake_minimum_required(VERSION 3.30)
project(acoustic_eigs)

find_package(Diskpp REQUIRED)

add_executable(acoustic_eigs src/acoustic_eigs.cpp)
target_link_libraries(acoustic_eigs diskpp::diskpp)
install(TARGETS acoustic_eigs RUNTIME DESTINATION bin)

