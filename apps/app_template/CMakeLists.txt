set(LINK_LIBS diskpp)

add_executable(app_template app_template.cpp)
target_link_libraries(app_template ${LINK_LIBS})
install(TARGETS app_template RUNTIME DESTINATION bin)
