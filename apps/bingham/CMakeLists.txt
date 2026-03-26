cmake_minimum_required(VERSION 3.30)
project(bingham)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(bingham_vector src/bingham_vector.cpp)
target_link_libraries(bingham_vector diskpp::diskpp)
install(TARGETS bingham_vector RUNTIME DESTINATION bin)

add_executable(bingham_vector_higher src/bingham_vector_higher.cpp)
target_link_libraries(bingham_vector_higher diskpp::diskpp)
install(TARGETS bingham_vector_higher RUNTIME DESTINATION bin)

add_executable(test_norm src/test_norm.cpp)
target_link_libraries(test_norm diskpp::diskpp)
install(TARGETS test_norm RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
