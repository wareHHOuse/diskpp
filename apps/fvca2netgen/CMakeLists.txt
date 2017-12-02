add_executable(fvca2netgen fvca2netgen.cpp)
target_link_libraries(fvca2netgen ${LINK_LIBS})
install(TARGETS fvca2netgen RUNTIME DESTINATION bin)
