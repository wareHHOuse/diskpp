set(LINK_LIBS diskpp)

add_executable(unsteady_laplacian src/unsteady_laplacian.cpp)
target_link_libraries(unsteady_laplacian ${LINK_LIBS})

add_executable(heat_dG src/heat_dG.cpp)
target_link_libraries(heat_dG ${LINK_LIBS})
