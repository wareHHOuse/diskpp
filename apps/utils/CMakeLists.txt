set(LINK_LIBS diskpp)

add_executable(mesh2silo mesh2silo.cpp)
target_link_libraries(mesh2silo ${LINK_LIBS})
install(TARGETS mesh2silo RUNTIME DESTINATION bin)

add_executable(make_fvca5_hexmesh make_fvca5_hexmesh.cpp)
target_link_libraries(make_fvca5_hexmesh ${LINK_LIBS})
install(TARGETS make_fvca5_hexmesh RUNTIME DESTINATION bin)

add_executable(inertia_transform_plot inertia_transform_plot.cpp)
target_link_libraries(inertia_transform_plot ${LINK_LIBS})
install(TARGETS inertia_transform_plot RUNTIME DESTINATION bin)

add_executable(inertia_eval_conditioning inertia_eval_conditioning.cpp)
target_link_libraries(inertia_eval_conditioning ${LINK_LIBS})
install(TARGETS inertia_eval_conditioning RUNTIME DESTINATION bin)