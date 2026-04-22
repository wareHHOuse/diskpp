cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(unsteady_laplacian src/unsteady_laplacian.cpp)
target_link_libraries(unsteady_laplacian diskpp::diskpp)

add_executable(heat_dG src/heat_dG.cpp)
target_link_libraries(heat_dG diskpp::diskpp)
