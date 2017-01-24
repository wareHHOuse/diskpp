add_executable(multiscale multiscale.cpp)
target_link_libraries(multiscale ${LINK_LIBS})
install(TARGETS multiscale RUNTIME DESTINATION bin)
