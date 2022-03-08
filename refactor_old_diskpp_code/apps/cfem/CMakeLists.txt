option(BUILD_APP_CFEM "Build Continuous FEM application" ON)
if (BUILD_APP_CFEM)
	add_executable(cfem cfem.cpp)
	target_link_libraries(cfem ${LINK_LIBS})
endif()

