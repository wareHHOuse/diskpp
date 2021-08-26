option(BUILD_APP_ULAPL "Build unsteady laplacian application" ON)
if (BUILD_APP_ULAPL)
    add_executable(unsteady_laplacian unsteady_laplacian.cpp)
    target_link_libraries(unsteady_laplacian ${LINK_LIBS})

    add_executable(heat_dG heat_dG.cpp)
    target_link_libraries(heat_dG ${LINK_LIBS})
endif ()
