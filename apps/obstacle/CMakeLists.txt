option(BUILD_OBSTACLE "Build Obstacle Problem application" ON)
if (BUILD_OBSTACLE)
add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle ${LINK_LIBS})
endif (BUILD_OBSTACLE)

