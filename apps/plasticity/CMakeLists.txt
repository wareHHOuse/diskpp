set(LINK_LIBS diskpp)

add_executable(plasticity src/plasticity.cpp)
target_link_libraries(plasticity ${LINK_LIBS})
install(TARGETS plasticity RUNTIME DESTINATION bin)

install(DIRECTORY gnuplot_scripts/ DESTINATION share)