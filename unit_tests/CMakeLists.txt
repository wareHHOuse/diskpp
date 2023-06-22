set(LINK_LIBS diskpp)

add_executable(monomial_bases_1D monomial_bases_1D.cpp)
target_link_libraries(monomial_bases_1D ${LINK_LIBS})

add_executable(monomial_bases_2D monomial_bases_2D.cpp)
target_link_libraries(monomial_bases_2D ${LINK_LIBS})

add_executable(quadratures quadratures.cpp)
target_link_libraries(quadratures ${LINK_LIBS})

add_executable(quad_raw_triangle_gauss quad_raw_triangle_gauss.cpp)
target_link_libraries(quad_raw_triangle_gauss ${LINK_LIBS})

add_executable(quad_raw_dunavant quad_raw_dunavant.cpp)
target_link_libraries(quad_raw_dunavant ${LINK_LIBS})

add_executable(quad_raw_tensorized_gl quad_raw_tensorized_gl.cpp)
target_link_libraries(quad_raw_tensorized_gl ${LINK_LIBS})