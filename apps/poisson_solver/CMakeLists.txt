cmake_minimum_required(VERSION 3.18)
project(poisson_solver)

set(LINK_LIBS diskpp)

add_executable(poisson_solver
    src/poisson_solver.cpp
    src/lua_interface.cpp
)
target_link_libraries(poisson_solver ${LINK_LIBS})
install(TARGETS poisson_solver RUNTIME DESTINATION bin)

add_executable(poisson_jump
    src/poisson_jump.cpp
)
target_link_libraries(poisson_jump ${LINK_LIBS})
install(TARGETS poisson_jump RUNTIME DESTINATION bin)