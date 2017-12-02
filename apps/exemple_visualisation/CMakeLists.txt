add_executable(exemple_visualisation  exemple_visualisation.cpp)
target_link_libraries(exemple_visualisation ${LINK_LIBS})
install(TARGETS exemple_visualisation RUNTIME DESTINATION bin)
