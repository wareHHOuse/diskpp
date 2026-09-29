cmake_minimum_required(VERSION 3.30)
project(static_coulomb)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wall_slide_v2 wall_slide_v2.cpp)
target_link_libraries(wall_slide_v2 diskpp::diskpp)
install(TARGETS wall_slide_v2 RUNTIME DESTINATION bin)

add_executable(sinum62_v2 sinum62_v2.cpp)
target_link_libraries(sinum62_v2 diskpp::diskpp)
install(TARGETS sinum62_v2 RUNTIME DESTINATION bin)

add_executable(coulomb3d_bostanhan coulomb3d_bostanhan.cpp)
target_link_libraries(coulomb3d_bostanhan diskpp::diskpp)
install(TARGETS coulomb3d_bostanhan RUNTIME DESTINATION bin)