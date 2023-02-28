set(LINK_LIBS diskpp)

add_executable(signorini_whho src/signorini_newton.cpp)
target_link_libraries(signorini_whho ${LINK_LIBS})

add_executable(signorini_fhho src/signorini_fix.cpp)
target_link_libraries(signorini_fhho ${LINK_LIBS})

add_executable(signorini_ffem src/signorini_cfem_fix.cpp)
target_link_libraries(signorini_ffem ${LINK_LIBS})

add_executable(signorini_wfem src/signorini_cfem_newton.cpp)
target_link_libraries(signorini_wfem ${LINK_LIBS})

add_executable(signorini_whier src/signorini_hierarchical_newton.cpp)
target_link_libraries(signorini_whier ${LINK_LIBS})

add_executable(diffusion_nitsche src/diffusion_nitsche.cpp)
target_link_libraries(diffusion_nitsche ${LINK_LIBS})

add_executable(diffusion_nitsche_test src/diffusion_nitsche_test.cpp)
target_link_libraries(diffusion_nitsche_test ${LINK_LIBS})

#add_executable(signorini_whho_par src/signorini_newton_parameter.cpp)
#target_link_libraries(signorini_whho_par ${LINK_LIBS})

add_executable(signorini_whho_test src/signorini_newton_test.cpp)
target_link_libraries(signorini_whho_test ${LINK_LIBS})
