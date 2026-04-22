cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle diskpp::diskpp)
install(TARGETS obstacle RUNTIME DESTINATION bin)

install(DIRECTORY gnuplot_scripts/ DESTINATION share)

