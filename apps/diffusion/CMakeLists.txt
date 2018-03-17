
add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion ${LINK_LIBS})
install(TARGETS diffusion RUNTIME DESTINATION bin)

add_executable(diffusion_test diffusion_test.cpp)
target_link_libraries(diffusion_test ${LINK_LIBS})
install(TARGETS diffusion_test RUNTIME DESTINATION bin)

add_executable(diffusion_new diffusion_new.cpp)
target_link_libraries(diffusion_new ${LINK_LIBS})
install(TARGETS diffusion_new RUNTIME DESTINATION bin)
