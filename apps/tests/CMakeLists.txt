
set(LINK_LIBS diskpp)

add_executable(monomial_bases_1D monomial_bases_1D.cpp)
target_link_libraries(monomial_bases_1D ${LINK_LIBS})

add_executable(monomial_bases_2D monomial_bases_2D.cpp)
target_link_libraries(monomial_bases_2D ${LINK_LIBS})

add_executable(bases_new bases_new.cpp)
target_link_libraries(bases_new ${LINK_LIBS})

add_executable(basis_function_plot basis_function_plot.cpp)
target_link_libraries(basis_function_plot ${LINK_LIBS})

add_executable(test_scalar_reconstruction scalar_reconstruction.cpp)
target_link_libraries(test_scalar_reconstruction ${LINK_LIBS})
install(TARGETS test_scalar_reconstruction RUNTIME DESTINATION bin)

    #add_executable(scalar_reconstruction_fullspace_test scalar_reconstruction_fullspace_test.cpp)
    #target_link_libraries(scalar_reconstruction_fullspace_test ${LINK_LIBS})

add_executable(test_scalar_stabilization scalar_stabilization.cpp)
target_link_libraries(test_scalar_stabilization ${LINK_LIBS})
install(TARGETS test_scalar_stabilization RUNTIME DESTINATION bin)

add_executable(test_vector_reconstruction vector_reconstruction.cpp)
target_link_libraries(test_vector_reconstruction ${LINK_LIBS})
install(TARGETS test_vector_reconstruction RUNTIME DESTINATION bin)

add_executable(test_divergence_reconstruction divergence_reconstruction.cpp)
target_link_libraries(test_divergence_reconstruction ${LINK_LIBS})

add_executable(test_vector_symmetric_reconstruction vector_symmetric_reconstruction.cpp)
target_link_libraries(test_vector_symmetric_reconstruction ${LINK_LIBS})

add_executable(test_vector_symmetric_stabilization vector_symmetric_stabilization.cpp)
target_link_libraries(test_vector_symmetric_stabilization ${LINK_LIBS})

add_executable(test_vector_stabilization vector_stabilization.cpp)
target_link_libraries(test_vector_stabilization ${LINK_LIBS})

add_executable(test_vector_gradient vector_gradient.cpp)
target_link_libraries(test_vector_gradient ${LINK_LIBS})

add_executable(test_perf perf.cpp)
target_link_libraries(test_perf ${LINK_LIBS})

add_executable(test_stokes stokes.cpp)
target_link_libraries(test_stokes ${LINK_LIBS})

add_executable(test_basis_quadratures basis_quadratures.cpp)
target_link_libraries(test_basis_quadratures ${LINK_LIBS})

    #add_executable(scalar_reconstruction_contact scalar_reconstruction_contact.cpp)
    #target_link_libraries(scalar_reconstruction_contact ${LINK_LIBS})

add_executable(test_raw_quadratures raw_quadratures.cpp)
target_link_libraries(test_raw_quadratures ${LINK_LIBS})

add_executable(test_quadratures_on_mesh quadratures_on_mesh.cpp)
target_link_libraries(test_quadratures_on_mesh ${LINK_LIBS})

add_executable(test_primitives primitives.cpp)
target_link_libraries(test_primitives ${LINK_LIBS})

add_executable(test_netgen_primitives netgen_primitives.cpp)
target_link_libraries(test_netgen_primitives ${LINK_LIBS})

add_executable(test_behaviors behaviors_test.cpp)
target_link_libraries(test_behaviors ${LINK_LIBS})

add_executable(test_integration integration.cpp)
target_link_libraries(test_integration ${LINK_LIBS})

add_executable(test_timing timing_test.cpp)
target_link_libraries(test_timing ${LINK_LIBS})

#add_executable(test_divide_polygons divide_polygons.cpp)
#target_link_libraries(test_divide_polygons ${LINK_LIBS})

add_executable(test_curl_reconstruction curl_reconstruction.cpp)
target_link_libraries(test_curl_reconstruction ${LINK_LIBS})

add_executable(test_curl_reconstruction_pk curl_reconstruction_pk.cpp)
target_link_libraries(test_curl_reconstruction_pk ${LINK_LIBS})

add_executable(test_curl_stabilization curl_stabilization.cpp)
target_link_libraries(test_curl_stabilization ${LINK_LIBS})

add_executable(test_curl_interpolation curl_interpolation.cpp)
target_link_libraries(test_curl_interpolation ${LINK_LIBS})

add_executable(test_curl_reconstruction_nedelec curl_reconstruction_nedelec.cpp)
target_link_libraries(test_curl_reconstruction_nedelec ${LINK_LIBS})

add_executable(test_curl_stabilization_nedelec curl_stabilization_nedelec.cpp)
target_link_libraries(test_curl_stabilization_nedelec ${LINK_LIBS})

add_executable(test_conditionning_basis conditionning_basis.cpp)
target_link_libraries(test_conditionning_basis ${LINK_LIBS})

#add_executable(test_conditionning test_conditionning.cpp)
#target_link_libraries(test_conditionning ${LINK_LIBS})

add_executable(meshgen meshgen.cpp)
target_link_libraries(meshgen ${LINK_LIBS})

if(HAVE_GMSH)
    add_executable(test_gmsh_geometry_loader gmsh_geometry_loader.cpp)
    target_link_libraries(test_gmsh_geometry_loader ${LINK_LIBS})
endif(HAVE_GMSH)
