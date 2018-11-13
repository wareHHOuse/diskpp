add_executable(diffusion_hho diffusion_hho.cpp)
target_link_libraries(diffusion_hho ${LINK_LIBS})
install(TARGETS diffusion_hho RUNTIME DESTINATION bin)

add_executable(diffusion_hho_test diffusion_hho_test.cpp)
target_link_libraries(diffusion_hho_test ${LINK_LIBS})
install(TARGETS diffusion_hho_test RUNTIME DESTINATION bin)

add_executable(diffusion_dga diffusion_dga.cpp)
target_link_libraries(diffusion_dga ${LINK_LIBS})
install(TARGETS diffusion_dga RUNTIME DESTINATION bin)

add_executable(diffusion_alg diffusion_alg.cpp)
target_link_libraries(diffusion_alg ${LINK_LIBS})
install(TARGETS diffusion_alg RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad diffusion_low_grad.cpp)
target_link_libraries(diffusion_low_grad ${LINK_LIBS})
install(TARGETS diffusion_low_grad RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad_vector diffusion_low_grad_vector.cpp)
target_link_libraries(diffusion_low_grad_vector ${LINK_LIBS})
install(TARGETS diffusion_low_grad_vector RUNTIME DESTINATION bin)
