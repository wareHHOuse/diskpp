add_executable(LerayLions LerayLions.cpp)
target_link_libraries(LerayLions ${LINK_LIBS})
install(TARGETS LerayLions RUNTIME DESTINATION bin)