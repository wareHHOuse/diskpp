add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion ${LINK_LIBS})
install(TARGETS diffusion RUNTIME DESTINATION bin)
