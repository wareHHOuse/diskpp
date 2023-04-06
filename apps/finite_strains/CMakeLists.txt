set(LINK_LIBS diskpp)

add_executable(finite_strains src/finite_strains.cpp)
target_link_libraries(finite_strains ${LINK_LIBS})
install(TARGETS finite_strains RUNTIME DESTINATION bin)
