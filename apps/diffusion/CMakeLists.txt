cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(sfexp SHARED src/stabfree_poly_explorer.cpp)
set_target_properties(sfexp PROPERTIES PREFIX "")
set_target_properties(sfexp PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(sfexp ${LUA_LIBRARIES} sol2::sol2 diskpp::diskpp)

#add_executable(stabfree_find_bad_poly src/stabfree_find_bad_poly.cpp)
#target_link_libraries(stabfree_find_bad_poly diskpp::diskpp)
#install(TARGETS stabfree_find_bad_poly RUNTIME DESTINATION bin)

add_library(rkernel SHARED src/reconstruction_kernel.cpp)
set_target_properties(rkernel PROPERTIES PREFIX "")
set_target_properties(rkernel PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(rkernel ${LUA_LIBRARIES} sol2::sol2 diskpp::diskpp)

add_executable(hpoly_test src/hpoly_test.cpp)
target_link_libraries(hpoly_test diskpp::diskpp)
install(TARGETS hpoly_test RUNTIME DESTINATION bin)

add_executable(diffusion_vem src/diffusion_vem.cpp)
target_link_libraries(diffusion_vem diskpp::diskpp)
install(TARGETS diffusion_vem RUNTIME DESTINATION bin)

add_executable(consistency_test src/consistency_test.cpp)
target_link_libraries(consistency_test diskpp::diskpp)
install(TARGETS consistency_test RUNTIME DESTINATION bin)

add_executable(consistency src/consistency.cpp)
target_link_libraries(consistency diskpp::diskpp)
install(TARGETS consistency RUNTIME DESTINATION bin)

add_executable(new_hho src/new_hho.cpp)
target_link_libraries(new_hho diskpp::diskpp)
install(TARGETS new_hho RUNTIME DESTINATION bin)

add_executable(diffusion_hho src/diffusion_hho.cpp)
target_link_libraries(diffusion_hho diskpp::diskpp)
install(TARGETS diffusion_hho RUNTIME DESTINATION bin)

add_executable(diffusion_hho_stabfree_test src/diffusion_hho_stabfree_test.cpp)
target_link_libraries(diffusion_hho_stabfree_test diskpp::diskpp)
install(TARGETS diffusion_hho_stabfree_test RUNTIME DESTINATION bin)

add_executable(hho_stabfree_eigs src/hho_stabfree_eigs.cpp)
target_link_libraries(hho_stabfree_eigs diskpp::diskpp)
install(TARGETS hho_stabfree_eigs RUNTIME DESTINATION bin)

add_executable(diffusion_dg src/diffusion_dg.cpp)
target_link_libraries(diffusion_dg diskpp::diskpp)
install(TARGETS diffusion_dg RUNTIME DESTINATION bin)

#add_executable(diffusion_hho_test src/diffusion_hho_test.cpp)
#target_link_libraries(diffusion_hho_test diskpp::diskpp)
#install(TARGETS diffusion_hho_test RUNTIME DESTINATION bin)

add_executable(diffusion_dga src/diffusion_dga.cpp)
target_link_libraries(diffusion_dga diskpp::diskpp)
install(TARGETS diffusion_dga RUNTIME DESTINATION bin)

add_executable(diffusion_alg src/diffusion_alg.cpp)
target_link_libraries(diffusion_alg diskpp::diskpp)
install(TARGETS diffusion_alg RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad src/diffusion_low_grad.cpp)
target_link_libraries(diffusion_low_grad diskpp::diskpp)
install(TARGETS diffusion_low_grad RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad_vector src/diffusion_low_grad_vector.cpp)
target_link_libraries(diffusion_low_grad_vector diskpp::diskpp)
install(TARGETS diffusion_low_grad_vector RUNTIME DESTINATION bin)

add_executable(variable_diffusion_hho src/variable_diffusion_hho.cpp)
target_link_libraries(variable_diffusion_hho diskpp::diskpp)
install(TARGETS variable_diffusion_hho RUNTIME DESTINATION bin)

