set(LINK_LIBS diskpp)

add_executable(barycenter_test barycenter_test.cpp)
target_link_libraries(barycenter_test ${LINK_LIBS})
install(TARGETS barycenter_test RUNTIME DESTINATION bin)
