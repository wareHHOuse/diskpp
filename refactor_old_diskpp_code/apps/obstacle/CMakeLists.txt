option(BUILD_APP_OBSTACLE "Build Obstacle Problem application" ON)
if (BUILD_APP_OBSTACLE)
	add_executable(obstacle obstacle.cpp)
	target_link_libraries(obstacle ${LINK_LIBS})
endif()

