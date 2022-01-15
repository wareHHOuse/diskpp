option(BUILD_APP_MAXWELL "Build Maxwell application" ON)
if (BUILD_APP_MAXWELL)
    add_executable(maxwell maxwell.cpp)
    target_link_libraries(maxwell ${LINK_LIBS})
    set_target_properties(maxwell PROPERTIES
                            CXX_STANDARD 20
                            CXX_STANDARD_REQUIRED YES
                            CXX_EXTENSIONS NO)

    add_executable(maxwell_sip_dg maxwell_sip_dg.cpp)
    target_link_libraries(maxwell_sip_dg ${LINK_LIBS})

    add_executable(rectest rectest.cpp)
    target_link_libraries(rectest ${LINK_LIBS})

    add_executable(meshio meshio.cpp)
    target_link_libraries(meshio ${LINK_LIBS})

    add_executable(maxwell_solver_hho maxwell_solver_hho.cpp)
    target_link_libraries(maxwell_solver_hho ${LINK_LIBS})
    set_target_properties(maxwell_solver_hho PROPERTIES
                            CXX_STANDARD 20
                            CXX_STANDARD_REQUIRED YES
                            CXX_EXTENSIONS NO)

    add_executable(maxwell_td maxwell_td.cpp)
    target_link_libraries(maxwell_td ${LINK_LIBS})
    set_target_properties(maxwell_td PROPERTIES
                            CXX_STANDARD 20
                            CXX_STANDARD_REQUIRED YES
                            CXX_EXTENSIONS NO)
endif()

