set(LINK_LIBS diskpp)

add_executable(wave_propagation src/wave_propagation.cpp)
target_link_libraries(wave_propagation ${LINK_LIBS})
install(TARGETS wave_propagation RUNTIME DESTINATION bin)
