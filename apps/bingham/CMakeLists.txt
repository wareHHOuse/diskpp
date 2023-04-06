set(LINK_LIBS diskpp)

add_executable(bingham_vector src/bingham_vector.cpp)
target_link_libraries(bingham_vector ${LINK_LIBS})
install(TARGETS bingham_vector RUNTIME DESTINATION bin)

add_executable(bingham_vector_higher src/bingham_vector_higher.cpp)
target_link_libraries(bingham_vector_higher ${LINK_LIBS})
install(TARGETS bingham_vector_higher RUNTIME DESTINATION bin)

add_executable(test_norm src/test_norm.cpp)
target_link_libraries(test_norm ${LINK_LIBS})
install(TARGETS test_norm RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)


