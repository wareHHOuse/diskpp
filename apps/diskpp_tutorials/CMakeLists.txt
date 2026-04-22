cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(diskpp_tutorial_1 diskpp_tutorial_1.cpp)
target_link_libraries(diskpp_tutorial_1 diskpp::diskpp)
install(TARGETS diskpp_tutorial_1 DESTINATION bin)

add_executable(diskpp_tutorial_2 diskpp_tutorial_2.cpp)
target_link_libraries(diskpp_tutorial_2 diskpp::diskpp)
install(TARGETS diskpp_tutorial_2 DESTINATION bin)

add_executable(diskpp_tutorial_3 diskpp_tutorial_3.cpp)
target_link_libraries(diskpp_tutorial_3 diskpp::diskpp)
install(TARGETS diskpp_tutorial_3 DESTINATION bin)

#if (HAVE_GMSH)
	add_executable(diskpp_tutorial_4 diskpp_tutorial_4.cpp)
	target_link_libraries(diskpp_tutorial_4 diskpp::diskpp)
	install(TARGETS diskpp_tutorial_4 DESTINATION bin)
#endif()