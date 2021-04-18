option(BUILD_DISKPP_TUTORIALS "Build DiSk++ tutorials" ON)
if (BUILD_DISKPP_TUTORIALS)
    add_executable(diskpp_tutorial_1 diskpp_tutorial_1.cpp)
    target_link_libraries(diskpp_tutorial_1 ${LINK_LIBS})

    add_executable(diskpp_tutorial_2 diskpp_tutorial_2.cpp)
    target_link_libraries(diskpp_tutorial_2 ${LINK_LIBS})

    add_executable(diskpp_tutorial_3 diskpp_tutorial_3.cpp)
    target_link_libraries(diskpp_tutorial_3 ${LINK_LIBS})
endif()
