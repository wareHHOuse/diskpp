set(LINK_LIBS diskpp)

add_executable(primitives primitives.cpp)
target_link_libraries(primitives ${LINK_LIBS})
add_test(NAME primitives COMMAND primitives)

add_executable(quad_raw_triangle_gauss quad_raw_triangle_gauss.cpp)
target_link_libraries(quad_raw_triangle_gauss ${LINK_LIBS})
add_test(NAME quad_raw_triangle_gauss COMMAND quad_raw_triangle_gauss)

add_executable(quad_raw_dunavant quad_raw_dunavant.cpp)
target_link_libraries(quad_raw_dunavant ${LINK_LIBS})
add_test(NAME quad_raw_dunavant COMMAND quad_raw_dunavant)

add_executable(quad_raw_tensorized_gl quad_raw_tensorized_gl.cpp)
target_link_libraries(quad_raw_tensorized_gl ${LINK_LIBS})
add_test(NAME quad_raw_tensorized_gl COMMAND quad_raw_tensorized_gl)

add_executable(quad_raw_tetra quad_raw_tetra.cpp)
target_link_libraries(quad_raw_tetra ${LINK_LIBS})
add_test(NAME quad_raw_tetra COMMAND quad_raw_tetra)

add_executable(quadratures quadratures.cpp)
target_link_libraries(quadratures ${LINK_LIBS})
add_test(NAME quadratures COMMAND quadratures)

add_executable(quadratures_divthm quadratures_divthm.cpp)
target_link_libraries(quadratures_divthm ${LINK_LIBS})
add_test(NAME quadratures_divthm COMMAND quadratures_divthm)

add_executable(monomial_bases_1D monomial_bases_1D.cpp)
target_link_libraries(monomial_bases_1D ${LINK_LIBS})
add_test(NAME monomial_bases_1D COMMAND monomial_bases_1D)

add_executable(monomial_bases_2D monomial_bases_2D.cpp)
target_link_libraries(monomial_bases_2D ${LINK_LIBS})
add_test(NAME monomial_bases_2D COMMAND monomial_bases_2D)
