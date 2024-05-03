set(LINK_LIBS diskpp)

add_executable(test_local_matrices src/test_local_matrices.cpp)
target_link_libraries(test_local_matrices ${LINK_LIBS})
install(TARGETS test_local_matrices RUNTIME DESTINATION bin)


install(DIRECTORY share/ DESTINATION share)
