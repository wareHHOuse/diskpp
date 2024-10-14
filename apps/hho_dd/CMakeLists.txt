set(LINK_LIBS diskpp)

add_executable(hhodd src/hhodd.cpp)
target_link_libraries(hhodd ${LINK_LIBS})
install(TARGETS hhodd RUNTIME DESTINATION bin)


install(DIRECTORY share/ DESTINATION share)
