
add_executable(diffusion2 diffusion2.cpp)
target_link_libraries(diffusion2 ${LINK_LIBS})
install(TARGETS diffusion2 RUNTIME DESTINATION bin)

add_executable(diffusion2_test diffusion2_test.cpp)
target_link_libraries(diffusion2_test ${LINK_LIBS})
install(TARGETS diffusion2_test RUNTIME DESTINATION bin)
