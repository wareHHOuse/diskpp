set(LINK_LIBS diskpp)

add_executable(conditioning src/conditioning.cpp)
target_link_libraries(conditioning ${LINK_LIBS})
install(TARGETS conditioning RUNTIME DESTINATION bin)