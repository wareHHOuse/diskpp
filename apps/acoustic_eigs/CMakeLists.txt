set(LINK_LIBS diskpp)

add_executable(acoustic_eigs src/acoustic_eigs.cpp)
target_link_libraries(acoustic_eigs ${LINK_LIBS})
install(TARGETS acoustic_eigs RUNTIME DESTINATION bin)
