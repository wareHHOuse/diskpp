set(LINK_LIBS diskpp)

add_executable(obstacle obstacle.cpp)
target_link_libraries(obstacle ${LINK_LIBS})
install(TARGETS obstacle RUNTIME DESTINATION bin)

install(DIRECTORY gnuplot_scripts/ DESTINATION share)

