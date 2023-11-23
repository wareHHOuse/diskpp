cmake_minimum_required(VERSION 3.18)
project(diskpp_steklov_eigsolver)

set(LINK_LIBS diskpp)

add_executable(steklov src/main.cpp)
target_link_libraries(steklov ${LINK_LIBS})
install(TARGETS steklov RUNTIME DESTINATION bin)