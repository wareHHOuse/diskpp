set(LINK_LIBS diskpp)

add_executable(mesh2silo mesh2silo.cpp)
target_link_libraries(mesh2silo ${LINK_LIBS})
install(TARGETS mesh2silo RUNTIME DESTINATION bin)