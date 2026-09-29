cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(sc_ad src/sc_ad.cpp)
target_link_libraries(sc_ad diskpp::diskpp)
install(TARGETS sc_ad RUNTIME DESTINATION bin)