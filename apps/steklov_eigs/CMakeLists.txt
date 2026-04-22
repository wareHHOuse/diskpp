cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(steklov src/main.cpp)
target_link_libraries(steklov diskpp::diskpp)
install(TARGETS steklov RUNTIME DESTINATION bin)