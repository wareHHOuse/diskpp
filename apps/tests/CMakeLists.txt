add_executable(scalar_reconstruction scalar_reconstruction.cpp)
target_link_libraries(scalar_reconstruction ${LINK_LIBS})

add_executable(scalar_stabilization scalar_stabilization.cpp)
target_link_libraries(scalar_stabilization ${LINK_LIBS})

add_executable(vector_reconstruction vector_reconstruction.cpp)
target_link_libraries(vector_reconstruction ${LINK_LIBS})


add_executable(divergence_reconstruction divergence_reconstruction.cpp)
target_link_libraries(divergence_reconstruction ${LINK_LIBS})

add_executable(vector_symmetric_reconstruction vector_symmetric_reconstruction.cpp)
target_link_libraries(vector_symmetric_reconstruction ${LINK_LIBS})

add_executable(vector_stabilization vector_stabilization.cpp)
target_link_libraries(vector_stabilization ${LINK_LIBS})

add_executable(perf perf.cpp)
target_link_libraries(perf ${LINK_LIBS})

add_executable(stokes_t stokes.cpp)
target_link_libraries(stokes_t ${LINK_LIBS})

add_executable(basis_quadratures basis_quadratures.cpp)
target_link_libraries(basis_quadratures ${LINK_LIBS})

add_executable(scalar_reconstruction_contact scalar_reconstruction_contact.cpp)
target_link_libraries(scalar_reconstruction_contact ${LINK_LIBS})

add_executable(raw_quadratures raw_quadratures.cpp)
target_link_libraries(raw_quadratures ${LINK_LIBS})

add_executable(quadratures_on_mesh quadratures_on_mesh.cpp)
target_link_libraries(quadratures_on_mesh ${LINK_LIBS})

add_executable(primitives primitives.cpp)
target_link_libraries(primitives ${LINK_LIBS})

add_executable(netgen_primitives netgen_primitives.cpp)
target_link_libraries(netgen_primitives ${LINK_LIBS})

