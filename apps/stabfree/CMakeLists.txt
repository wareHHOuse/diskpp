set(LINK_LIBS diskpp)

add_executable(stabfree stabfree.cpp)
target_link_libraries(stabfree ${LINK_LIBS})

