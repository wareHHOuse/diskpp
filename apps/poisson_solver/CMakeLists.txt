cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(poisson_solver
    src/poisson_solver.cpp
    src/lua_interface.cpp
)
target_link_libraries(poisson_solver diskpp::diskpp)
install(TARGETS poisson_solver RUNTIME DESTINATION bin)

add_executable(poisson_jump
    src/poisson_jump.cpp
)
target_link_libraries(poisson_jump diskpp::diskpp)
install(TARGETS poisson_jump RUNTIME DESTINATION bin)