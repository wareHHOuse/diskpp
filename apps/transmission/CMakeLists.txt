cmake_minimum_required(VERSION 3.18)
project(diskpp_tutorials)

set(LINK_LIBS diskpp)

add_executable(linear_transmission linear_transmission.cpp)
target_link_libraries(linear_transmission ${LINK_LIBS})
install(TARGETS linear_transmission RUNTIME DESTINATION bin)

add_executable(nonlinear_transmission nonlinear_transmission.cpp)
target_link_libraries(nonlinear_transmission ${LINK_LIBS})

add_executable(test_projections test_projections.cpp)
target_link_libraries(test_projections ${LINK_LIBS})

add_executable(pardisotest pardisotest.cpp)
target_link_libraries(pardisotest ${LINK_LIBS})