cmake_minimum_required(VERSION 3.18)
project(diskpp_diffusion)

set(LINK_LIBS diskpp)

add_executable(diffusion_vem src/diffusion_vem.cpp)
target_link_libraries(diffusion_vem ${LINK_LIBS})
install(TARGETS diffusion_vem RUNTIME DESTINATION bin)

add_executable(consistency src/consistency.cpp)
target_link_libraries(consistency ${LINK_LIBS})
install(TARGETS consistency RUNTIME DESTINATION bin)

add_executable(new_hho src/new_hho.cpp)
target_link_libraries(new_hho ${LINK_LIBS})
install(TARGETS new_hho RUNTIME DESTINATION bin)

add_executable(diffusion_hho src/diffusion_hho.cpp)
target_link_libraries(diffusion_hho ${LINK_LIBS})
install(TARGETS diffusion_hho RUNTIME DESTINATION bin)

add_executable(diffusion_hho_stabfree_test src/diffusion_hho_stabfree_test.cpp)
target_link_libraries(diffusion_hho_stabfree_test ${LINK_LIBS})
install(TARGETS diffusion_hho_stabfree_test RUNTIME DESTINATION bin)

add_executable(hho_stabfree_eigs src/hho_stabfree_eigs.cpp)
target_link_libraries(hho_stabfree_eigs ${LINK_LIBS})
install(TARGETS hho_stabfree_eigs RUNTIME DESTINATION bin)

add_executable(diffusion_dg src/diffusion_dg.cpp)
target_link_libraries(diffusion_dg ${LINK_LIBS})
install(TARGETS diffusion_dg RUNTIME DESTINATION bin)

#add_executable(diffusion_hho_test src/diffusion_hho_test.cpp)
#target_link_libraries(diffusion_hho_test ${LINK_LIBS})
#install(TARGETS diffusion_hho_test RUNTIME DESTINATION bin)

add_executable(diffusion_dga src/diffusion_dga.cpp)
target_link_libraries(diffusion_dga ${LINK_LIBS})
install(TARGETS diffusion_dga RUNTIME DESTINATION bin)

add_executable(diffusion_alg src/diffusion_alg.cpp)
target_link_libraries(diffusion_alg ${LINK_LIBS})
install(TARGETS diffusion_alg RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad src/diffusion_low_grad.cpp)
target_link_libraries(diffusion_low_grad ${LINK_LIBS})
install(TARGETS diffusion_low_grad RUNTIME DESTINATION bin)

add_executable(diffusion_low_grad_vector src/diffusion_low_grad_vector.cpp)
target_link_libraries(diffusion_low_grad_vector ${LINK_LIBS})
install(TARGETS diffusion_low_grad_vector RUNTIME DESTINATION bin)

add_executable(variable_diffusion_hho src/variable_diffusion_hho.cpp)
target_link_libraries(variable_diffusion_hho ${LINK_LIBS})
install(TARGETS variable_diffusion_hho RUNTIME DESTINATION bin)

