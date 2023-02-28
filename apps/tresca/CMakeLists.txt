set(LINK_LIBS diskpp)

add_executable(tresca src/tresca.cpp)
target_link_libraries(tresca ${LINK_LIBS})
install(TARGETS tresca RUNTIME DESTINATION bin)

add_executable(tresca_test src/tresca_test.cpp)
target_link_libraries(tresca_test ${LINK_LIBS})
install(TARGETS tresca_test RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)


