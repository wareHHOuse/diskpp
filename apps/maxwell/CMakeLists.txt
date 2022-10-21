set(LINK_LIBS diskpp)

add_executable(maxwell src/maxwell.cpp)
target_link_libraries(maxwell ${LINK_LIBS})
set_target_properties(maxwell PROPERTIES
                        CXX_STANDARD 20
                        CXX_STANDARD_REQUIRED YES
                        CXX_EXTENSIONS NO)
install(TARGETS maxwell RUNTIME DESTINATION bin)

add_executable(maxwell_sip_dg src/maxwell_sip_dg.cpp)
target_link_libraries(maxwell_sip_dg ${LINK_LIBS})
install(TARGETS maxwell RUNTIME DESTINATION bin)

add_executable(rectest src/rectest.cpp)
target_link_libraries(rectest ${LINK_LIBS})
install(TARGETS rectest RUNTIME DESTINATION bin)

add_executable(meshio src/meshio.cpp)
target_link_libraries(meshio ${LINK_LIBS})
install(TARGETS meshio RUNTIME DESTINATION bin)

add_executable(maxwell_solver_hho src/maxwell_solver_hho.cpp)
target_link_libraries(maxwell_solver_hho ${LINK_LIBS})
set_target_properties(maxwell_solver_hho PROPERTIES
                        CXX_STANDARD 20
                        CXX_STANDARD_REQUIRED YES
                        CXX_EXTENSIONS NO)
install(TARGETS maxwell_solver_hho RUNTIME DESTINATION bin)

add_executable(maxwell_td src/maxwell_td.cpp)
target_link_libraries(maxwell_td ${LINK_LIBS})
set_target_properties(maxwell_td PROPERTIES
                        CXX_STANDARD 20
                        CXX_STANDARD_REQUIRED YES
                        CXX_EXTENSIONS NO)
install(TARGETS maxwell_td RUNTIME DESTINATION bin)

install(DIRECTORY share/ DESTINATION share)