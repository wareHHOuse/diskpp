cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(nitsche src/nitsche.cpp)
target_link_libraries(nitsche diskpp::diskpp)
install(TARGETS nitsche RUNTIME DESTINATION bin)

add_executable(linelast src/linelast.cpp)
target_link_libraries(linelast diskpp::diskpp)
install(TARGETS linelast RUNTIME DESTINATION bin)

#add_executable(test_divrec src/test_divrec.cpp)
#target_link_libraries(test_divrec diskpp::diskpp)


add_executable(bridge src/bridge.cpp)
target_link_libraries(bridge diskpp::diskpp)