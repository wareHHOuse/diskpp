set(LINK_LIBS diskpp)

add_executable(stokes stokes.cpp)
target_link_libraries(stokes ${LINK_LIBS})
install(TARGETS stokes RUNTIME DESTINATION bin)

add_executable(stokes_driven stokes_driven_cavity.cpp)
target_link_libraries(stokes_driven ${LINK_LIBS})
install(TARGETS stokes_driven RUNTIME DESTINATION bin)

