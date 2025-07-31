set(LINK_LIBS diskpp)

add_executable(nitsche src/nitsche.cpp)
target_link_libraries(nitsche ${LINK_LIBS})
install(TARGETS nitsche RUNTIME DESTINATION bin)