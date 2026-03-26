cmake_minimum_required(VERSION 3.30)
project(app_template)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(mesh2silo mesh2silo.cpp)
target_link_libraries(mesh2silo diskpp::diskpp)
install(TARGETS mesh2silo RUNTIME DESTINATION bin)

add_executable(make_fvca5_hexmesh make_fvca5_hexmesh.cpp)
target_link_libraries(make_fvca5_hexmesh diskpp::diskpp)
install(TARGETS make_fvca5_hexmesh RUNTIME DESTINATION bin)

add_executable(inertia_transform_plot inertia_transform_plot.cpp)
target_link_libraries(inertia_transform_plot diskpp::diskpp)
install(TARGETS inertia_transform_plot RUNTIME DESTINATION bin)

add_executable(inertia_eval_conditioning inertia_eval_conditioning.cpp)
target_link_libraries(inertia_eval_conditioning diskpp::diskpp)
install(TARGETS inertia_eval_conditioning RUNTIME DESTINATION bin)

add_executable(mkgauss mkgauss.cpp)
target_link_libraries(mkgauss diskpp::diskpp)
install(TARGETS mkgauss RUNTIME DESTINATION bin)

add_executable(check_triangle_output check_triangle_output.cpp)
target_link_libraries(check_triangle_output diskpp::diskpp)
install(TARGETS check_triangle_output RUNTIME DESTINATION bin)

add_executable(h5io h5io.cpp)
target_link_libraries(h5io diskpp::diskpp)
install(TARGETS h5io RUNTIME DESTINATION bin)