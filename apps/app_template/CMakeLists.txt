cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(app_template app_template.cpp)
target_link_libraries(app_template diskpp::diskpp)
install(TARGETS app_template RUNTIME DESTINATION bin)
