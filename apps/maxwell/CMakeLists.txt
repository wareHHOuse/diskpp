cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(maxwell_dga src/maxwell_dga.cpp)
target_link_libraries(maxwell_dga diskpp::diskpp)
install(TARGETS maxwell_dga RUNTIME DESTINATION bin)

add_executable(maxwell src/maxwell.cpp)
target_link_libraries(maxwell diskpp::diskpp)
install(TARGETS maxwell RUNTIME DESTINATION bin)

add_executable(maxwell_sip_dg src/maxwell_sip_dg.cpp)
target_link_libraries(maxwell_sip_dg diskpp::diskpp)
install(TARGETS maxwell_sip_dg RUNTIME DESTINATION bin)

add_executable(rectest src/rectest.cpp)
target_link_libraries(rectest diskpp::diskpp)
install(TARGETS rectest RUNTIME DESTINATION bin)

add_executable(maxwell_solver_hho src/maxwell_solver_hho.cpp)
target_link_libraries(maxwell_solver_hho diskpp::diskpp)
install(TARGETS maxwell_solver_hho RUNTIME DESTINATION bin)

add_executable(maxwell_td src/maxwell_td.cpp)
target_link_libraries(maxwell_td diskpp::diskpp)
install(TARGETS maxwell_td RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
