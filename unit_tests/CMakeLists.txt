set(LINK_LIBS diskpp)

add_executable(monomial_bases_1D monomial_bases_1D.cpp)
target_link_libraries(monomial_bases_1D ${LINK_LIBS})

add_executable(monomial_bases_2D monomial_bases_2D.cpp)
target_link_libraries(monomial_bases_2D ${LINK_LIBS})