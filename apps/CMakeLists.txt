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

option(BUILD_APP_MAXWELL "Build Maxwell application" ON)
if (BUILD_APP_MAXWELL)
    add_subdirectory(maxwell)
endif()

option(BUILD_APP_TRANSMISSION "Build transmission applications" ON)
if (BUILD_APP_TRANSMISSION)
    add_subdirectory(transmission)
endif()

option(BUILD_APP_OBSTACLE "Build obstacle application" ON)
if (BUILD_APP_OBSTACLE)
    add_subdirectory(obstacle)
endif()

option(BUILD_APP_BINGHAM "Build Bingham applications" ON)
if (BUILD_APP_BINGHAM)
    add_subdirectory(bingham)
endif()

option(BUILD_APP_HELMHOLTZ "Build Helmholtz applications" ON)
if (BUILD_APP_HELMHOLTZ)
    add_subdirectory(helmholtz)
endif()

option(BUILD_APP_CONTACT "Build contact applications" ON)
if (BUILD_APP_CONTACT)
    add_subdirectory(contact)
endif()

