set(LINK_LIBS diskpp)

add_executable(cfem src/cfem.cpp)
target_link_libraries(cfem ${LINK_LIBS})
install(TARGETS cfem RUNTIME DESTINATION bin)