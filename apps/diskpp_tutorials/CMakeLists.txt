cmake_minimum_required(VERSION 3.18)
project(diskpp_tutorials)

set(LINK_LIBS diskpp)

add_executable(diskpp_tutorial_1 diskpp_tutorial_1.cpp)
target_link_libraries(diskpp_tutorial_1 ${LINK_LIBS})
install(TARGETS diskpp_tutorial_1 DESTINATION bin)

add_executable(diskpp_tutorial_2 diskpp_tutorial_2.cpp)
target_link_libraries(diskpp_tutorial_2 ${LINK_LIBS})
install(TARGETS diskpp_tutorial_2 DESTINATION bin)

add_executable(diskpp_tutorial_3 diskpp_tutorial_3.cpp)
target_link_libraries(diskpp_tutorial_3 ${LINK_LIBS})
install(TARGETS diskpp_tutorial_3 DESTINATION bin)


