
add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion ${LINK_LIBS})
install(TARGETS diffusion RUNTIME DESTINATION bin)

add_executable(diffusion_test diffusion_test.cpp)
target_link_libraries(diffusion_test ${LINK_LIBS})
install(TARGETS diffusion_test RUNTIME DESTINATION bin)

add_executable(diffusion_hho diffusion_hho.cpp)
target_link_libraries(diffusion_hho ${LINK_LIBS})
install(TARGETS diffusion_hho RUNTIME DESTINATION bin)

add_executable(diffusion_dga diffusion_dga.cpp)
target_link_libraries(diffusion_dga ${LINK_LIBS})
install(TARGETS diffusion_dga RUNTIME DESTINATION bin)
