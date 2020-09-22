option(BUILD_APP_MAXWELL "Build Maxwell application" ON)
if (BUILD_APP_MAXWELL)
    add_executable(maxwell maxwell.cpp)
    target_link_libraries(maxwell ${LINK_LIBS})
    set_target_properties(maxwell PROPERTIES
                            CXX_STANDARD 20
                            CXX_STANDARD_REQUIRED YES
                            CXX_EXTENSIONS NO)
endif()

