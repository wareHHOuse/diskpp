cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(vector_laplacian src/vector_laplacian.cpp)
target_link_libraries(vector_laplacian diskpp::diskpp)
install(TARGETS vector_laplacian RUNTIME DESTINATION bin)

add_executable(vector_laplacian_test src/vector_laplacian_test.cpp)
target_link_libraries(vector_laplacian_test diskpp::diskpp)
install(TARGETS vector_laplacian_test RUNTIME DESTINATION bin)
