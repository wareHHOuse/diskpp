set(LINK_LIBS diskpp)

add_executable(vector_laplacian src/vector_laplacian.cpp)
target_link_libraries(vector_laplacian ${LINK_LIBS})
install(TARGETS vector_laplacian RUNTIME DESTINATION bin)

add_executable(vector_laplacian_test src/vector_laplacian_test.cpp)
target_link_libraries(vector_laplacian_test ${LINK_LIBS})
install(TARGETS vector_laplacian_test RUNTIME DESTINATION bin)
