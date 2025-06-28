set(LINK_LIBS diskpp)

add_executable(unfitted_HHO src/unfitted_HHO.cpp)
target_link_libraries(unfitted_HHO ${LINK_LIBS})
install(TARGETS unfitted_HHO RUNTIME DESTINATION bin)
