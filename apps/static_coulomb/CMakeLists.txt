cmake_minimum_required(VERSION 3.30)
project(static_coulomb)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(coulomb_v2 coulomb_v2.cpp)
target_link_libraries(coulomb_v2 diskpp::diskpp)
install(TARGETS coulomb_v2 RUNTIME DESTINATION bin)
 
add_executable(araya25 araya25.cpp)
target_link_libraries(araya25 diskpp::diskpp)
install(TARGETS araya25 RUNTIME DESTINATION bin)

add_executable(wall_slide wall_slide.cpp)
target_link_libraries(wall_slide diskpp::diskpp)
install(TARGETS wall_slide RUNTIME DESTINATION bin)

add_executable(AMAEPNTN_522_2D AMAEPNTN_522_2D.cpp)
target_link_libraries(AMAEPNTN_522_2D diskpp::diskpp)
install(TARGETS AMAEPNTN_522_2D RUNTIME DESTINATION bin)

add_executable(AMAEPNTN_522_3D AMAEPNTN_522_3D.cpp)
target_link_libraries(AMAEPNTN_522_3D diskpp::diskpp)
install(TARGETS AMAEPNTN_522_3D RUNTIME DESTINATION bin)

add_executable(rahmi_422 rahmi_422.cpp)
target_link_libraries(rahmi_422 diskpp::diskpp)
install(TARGETS rahmi_422 RUNTIME DESTINATION bin)