cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(vemasm src/vemasm.cpp)
target_link_libraries(vemasm diskpp::diskpp)
install(TARGETS vemasm RUNTIME DESTINATION bin)
