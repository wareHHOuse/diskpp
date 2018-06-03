add_executable(stokes stokes.cpp)
target_link_libraries(stokes ${LINK_LIBS})

add_executable(stokes_low_grad stokes_low_grad.cpp)
target_link_libraries(stokes_low_grad ${LINK_LIBS})

add_executable(stokes_driven stokes_alg_driven_cavity.cpp)
target_link_libraries(stokes_driven ${LINK_LIBS})

add_executable(stokes_alg stokes_alg.cpp)
target_link_libraries(stokes_alg ${LINK_LIBS})

add_executable(stokes_alg_sym stokes_alg_sym.cpp)
target_link_libraries(stokes_alg_sym ${LINK_LIBS})

#add_executable(stokes_alg_no_div stokes_alg_no_div.cpp)
#target_link_libraries(stokes_alg_no_div ${LINK_LIBS})

add_executable(stokes_cylinder stokes_over_cylinder.cpp)
target_link_libraries(stokes_cylinder ${LINK_LIBS})
