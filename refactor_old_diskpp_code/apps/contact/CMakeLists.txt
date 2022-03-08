option(BUILD_APP_CONTACT "Build Contact applications" ON)
if (BUILD_APP_CONTACT)
	add_executable(signorini_whho signorini_newton.cpp)
	target_link_libraries(signorini_whho ${LINK_LIBS})

	add_executable(signorini_fhho signorini_fix.cpp)
	target_link_libraries(signorini_fhho ${LINK_LIBS})

	add_executable(signorini_ffem signorini_cfem_fix.cpp)
	target_link_libraries(signorini_ffem ${LINK_LIBS})

	add_executable(signorini_wfem signorini_cfem_newton.cpp)
	target_link_libraries(signorini_wfem ${LINK_LIBS})

	add_executable(signorini_whier signorini_hierarchical_newton.cpp)
	target_link_libraries(signorini_whier ${LINK_LIBS})

	add_executable(diffusion_nitsche diffusion_nitsche.cpp)
	target_link_libraries(diffusion_nitsche ${LINK_LIBS})

	add_executable(diffusion_nitsche_test diffusion_nitsche_test.cpp)
	target_link_libraries(diffusion_nitsche_test ${LINK_LIBS})

	#add_executable(signorini_whho_par signorini_newton_parameter.cpp)
	#target_link_libraries(signorini_whho_par ${LINK_LIBS})

	add_executable(signorini_whho_test signorini_newton_test.cpp)
	target_link_libraries(signorini_whho_test ${LINK_LIBS})
endif()
