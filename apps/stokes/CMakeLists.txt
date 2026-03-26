cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(stokes stokes.cpp)
target_link_libraries(stokes diskpp::diskpp)
install(TARGETS stokes RUNTIME DESTINATION bin)

add_executable(stokes_driven stokes_driven_cavity.cpp)
target_link_libraries(stokes_driven diskpp::diskpp)
install(TARGETS stokes_driven RUNTIME DESTINATION bin)

