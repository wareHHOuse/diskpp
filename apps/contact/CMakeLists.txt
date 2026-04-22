cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(signorini_whho src/signorini_newton.cpp)
target_link_libraries(signorini_whho diskpp::diskpp)
install(TARGETS signorini_whho RUNTIME DESTINATION bin)

add_executable(signorini_fhho src/signorini_fix.cpp)
target_link_libraries(signorini_fhho diskpp::diskpp)
install(TARGETS signorini_fhho RUNTIME DESTINATION bin)

add_executable(signorini_ffem src/signorini_cfem_fix.cpp)
target_link_libraries(signorini_ffem diskpp::diskpp)
install(TARGETS signorini_ffem RUNTIME DESTINATION bin)

add_executable(signorini_wfem src/signorini_cfem_newton.cpp)
target_link_libraries(signorini_wfem diskpp::diskpp)
install(TARGETS signorini_wfem RUNTIME DESTINATION bin)

add_executable(signorini_whier src/signorini_hierarchical_newton.cpp)
target_link_libraries(signorini_whier diskpp::diskpp)
install(TARGETS signorini_whier RUNTIME DESTINATION bin)

add_executable(diffusion_nitsche src/diffusion_nitsche.cpp)
target_link_libraries(diffusion_nitsche diskpp::diskpp)
install(TARGETS diffusion_nitsche RUNTIME DESTINATION bin)

add_executable(diffusion_nitsche_test src/diffusion_nitsche_test.cpp)
target_link_libraries(diffusion_nitsche_test diskpp::diskpp)
install(TARGETS diffusion_nitsche_test RUNTIME DESTINATION bin)

#add_executable(signorini_whho_par src/signorini_newton_parameter.cpp)
#target_link_libraries(signorini_whho_par ${LINK_LIBS})

add_executable(signorini_whho_test src/signorini_newton_test.cpp)
target_link_libraries(signorini_whho_test diskpp::diskpp)
install(TARGETS signorini_whho_test RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)
