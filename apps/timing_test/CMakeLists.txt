add_executable(timing_test timing_test.cpp)
target_link_libraries(timing_test ${LINK_LIBS})
install(TARGETS timing_test RUNTIME DESTINATION bin)
