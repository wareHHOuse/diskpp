add_executable(signorini_whho signorini_newton.cpp)
target_link_libraries(signorini_whho ${LINK_LIBS})

add_executable(signorini_fhho signorini_fix.cpp)
target_link_libraries(signorini_fhho ${LINK_LIBS})

add_executable(signorini_ffem signorini_cfem_fix.cpp)
target_link_libraries(signorini_ffem ${LINK_LIBS})

add_executable(signorini_wfem signorini_cfem_newton.cpp)
target_link_libraries(signorini_wfem ${LINK_LIBS})

add_executable(signorini_whier signorini_hierarchical_newton.cpp)
target_link_libraries(signorini_whier ${LINK_LIBS})
