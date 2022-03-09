option(BUILD_DISKPP_TUTORIALS "Build DiSk++ tutorials" ON)
if (BUILD_DISKPP_TUTORIALS)
    add_subdirectory(diskpp_tutorials)
endif()

option(BUILD_APP_DIFFUSION "Build Diffusion applications" ON)
if (BUILD_APP_DIFFUSION)
    add_subdirectory(diffusion)
endif()

option(BUILD_APP_DISKPP_TESTS "Build DiSk++ tests" ON)
if (BUILD_APP_DISKPP_TESTS)
    add_subdirectory(tests)
endif()