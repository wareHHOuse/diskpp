add_executable(multiscale multiscale.cpp)
target_link_libraries(multiscale ${LINK_LIBS})
install(TARGETS multiscale RUNTIME DESTINATION bin)

add_executable(reference_solver reference_solver.cpp)
target_link_libraries(reference_solver ${LINK_LIBS})
install(TARGETS reference_solver RUNTIME DESTINATION bin)
