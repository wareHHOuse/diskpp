set(LINK_LIBS diskpp)

add_executable(nitsche src/nitsche.cpp)
target_link_libraries(nitsche ${LINK_LIBS})
install(TARGETS nitsche RUNTIME DESTINATION bin)

add_executable(linelast src/linelast.cpp)
target_link_libraries(linelast ${LINK_LIBS})
install(TARGETS linelast RUNTIME DESTINATION bin)
