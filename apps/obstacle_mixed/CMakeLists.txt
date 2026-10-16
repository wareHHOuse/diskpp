set(LINK_LIBS diskpp)

## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(membranes_Pk_Pk src/membranes_Pk_Pk.cpp)
target_link_libraries(membranes_Pk_Pk ${LINK_LIBS})
install(TARGETS membranes_Pk_Pk RUNTIME DESTINATION bin)

add_executable(obstacle_Pk_Pk src/obstacle_Pk_Pk.cpp)
target_link_libraries(obstacle_Pk_Pk ${LINK_LIBS})
install(TARGETS obstacle_Pk_Pk RUNTIME DESTINATION bin)

add_executable(obstacle_P2_P0 src/obstacle_P2_P0.cpp)
target_link_libraries(obstacle_P2_P0 ${LINK_LIBS})
install(TARGETS obstacle_P2_P0 RUNTIME DESTINATION bin)
