set(LINK_LIBS diskpp)

add_executable(run_pardiso src/run_pardiso.cpp)
target_link_libraries(run_pardiso ${LINK_LIBS})
install(TARGETS run_pardiso RUNTIME DESTINATION bin)
