add_executable(plasticity plasticity.cpp)
target_link_libraries(plasticity ${LINK_LIBS})
install(TARGETS plasticity RUNTIME DESTINATION bin)
