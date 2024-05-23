set(LINK_LIBS diskpp)

add_executable(vemasm src/vemasm.cpp)
target_link_libraries(vemasm ${LINK_LIBS})
install(TARGETS vemasm RUNTIME DESTINATION bin)
