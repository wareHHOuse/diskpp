cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(elemtriangulation elemtriangulation.cpp)
target_link_libraries(elemtriangulation diskpp::diskpp)

add_executable(bases_new bases_new.cpp)
target_link_libraries(bases_new diskpp::diskpp)

add_executable(basis_function_plot basis_function_plot.cpp)
target_link_libraries(basis_function_plot diskpp::diskpp)

add_executable(test_scalar_reconstruction scalar_reconstruction.cpp)
target_link_libraries(test_scalar_reconstruction diskpp::diskpp)
install(TARGETS test_scalar_reconstruction RUNTIME DESTINATION bin)

    #add_executable(scalar_reconstruction_fullspace_test scalar_reconstruction_fullspace_test.cpp)
    #target_link_libraries(scalar_reconstruction_fullspace_test diskpp::diskpp)

add_executable(test_scalar_stabilization scalar_stabilization.cpp)
target_link_libraries(test_scalar_stabilization diskpp::diskpp)
install(TARGETS test_scalar_stabilization RUNTIME DESTINATION bin)

add_executable(test_vector_reconstruction vector_reconstruction.cpp)
target_link_libraries(test_vector_reconstruction diskpp::diskpp)
install(TARGETS test_vector_reconstruction RUNTIME DESTINATION bin)

add_executable(test_divergence_reconstruction divergence_reconstruction.cpp)
target_link_libraries(test_divergence_reconstruction diskpp::diskpp)

add_executable(test_vector_symmetric_reconstruction vector_symmetric_reconstruction.cpp)
target_link_libraries(test_vector_symmetric_reconstruction diskpp::diskpp)

add_executable(test_vector_symmetric_stabilization vector_symmetric_stabilization.cpp)
target_link_libraries(test_vector_symmetric_stabilization diskpp::diskpp)

add_executable(test_vector_stabilization vector_stabilization.cpp)
target_link_libraries(test_vector_stabilization diskpp::diskpp)

add_executable(test_vector_gradient vector_gradient.cpp)
target_link_libraries(test_vector_gradient diskpp::diskpp)

add_executable(test_perf perf.cpp)
target_link_libraries(test_perf diskpp::diskpp)

add_executable(test_stokes stokes.cpp)
target_link_libraries(test_stokes diskpp::diskpp)

add_executable(test_basis_quadratures basis_quadratures.cpp)
target_link_libraries(test_basis_quadratures diskpp::diskpp)

    #add_executable(scalar_reconstruction_contact scalar_reconstruction_contact.cpp)
    #target_link_libraries(scalar_reconstruction_contact diskpp::diskpp)

add_executable(test_raw_quadratures raw_quadratures.cpp)
target_link_libraries(test_raw_quadratures diskpp::diskpp)

add_executable(test_quadratures_on_mesh quadratures_on_mesh.cpp)
target_link_libraries(test_quadratures_on_mesh diskpp::diskpp)

add_executable(test_primitives primitives.cpp)
target_link_libraries(test_primitives diskpp::diskpp)

add_executable(test_netgen_primitives netgen_primitives.cpp)
target_link_libraries(test_netgen_primitives diskpp::diskpp)

add_executable(test_behaviors behaviors_test.cpp)
target_link_libraries(test_behaviors diskpp::diskpp)

add_executable(test_integration integration.cpp)
target_link_libraries(test_integration diskpp::diskpp)

add_executable(test_timing timing_test.cpp)
target_link_libraries(test_timing diskpp::diskpp)

#add_executable(test_divide_polygons divide_polygons.cpp)
#target_link_libraries(test_divide_polygons diskpp::diskpp)

add_executable(test_curl_reconstruction curl_reconstruction.cpp)
target_link_libraries(test_curl_reconstruction diskpp::diskpp)

add_executable(test_curl_reconstruction_pk curl_reconstruction_pk.cpp)
target_link_libraries(test_curl_reconstruction_pk diskpp::diskpp)

add_executable(test_curl_stabilization curl_stabilization.cpp)
target_link_libraries(test_curl_stabilization diskpp::diskpp)

add_executable(test_curl_interpolation curl_interpolation.cpp)
target_link_libraries(test_curl_interpolation diskpp::diskpp)

add_executable(test_curl_reconstruction_nedelec curl_reconstruction_nedelec.cpp)
target_link_libraries(test_curl_reconstruction_nedelec diskpp::diskpp)

add_executable(test_curl_stabilization_nedelec curl_stabilization_nedelec.cpp)
target_link_libraries(test_curl_stabilization_nedelec diskpp::diskpp)

add_executable(test_conditionning_basis conditionning_basis.cpp)
target_link_libraries(test_conditionning_basis diskpp::diskpp)

#add_executable(test_conditionning test_conditionning.cpp)
#target_link_libraries(test_conditionning diskpp::diskpp)

add_executable(meshgen meshgen.cpp)
target_link_libraries(meshgen diskpp::diskpp)

if(HAVE_GMSH)
    add_executable(test_gmsh_geometry_loader gmsh_geometry_loader.cpp)
    target_link_libraries(test_gmsh_geometry_loader diskpp::diskpp)
endif(HAVE_GMSH)
