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

option(BUILD_APP_LINEAR_ELASTICITY "Build linear elasticity applications" ON)
if (BUILD_APP_LINEAR_ELASTICITY)
    add_subdirectory(linear_elasticity)
endif()

option(BUILD_APP_CONTACT "Build contact applications" ON)
if (BUILD_APP_CONTACT)
    add_subdirectory(contact)
endif()

option(BUILD_APP_STOKES "Build Stokes applications" ON)
if (BUILD_APP_STOKES)
    add_subdirectory(stokes)
endif()

option(BUILD_APP_TEMPLATE "Build template applications" ON)
if (BUILD_APP_TEMPLATE)
    add_subdirectory(app_template)
endif()

option(BUILD_APP_CFEM "Build Continuous FEM application" ON)
if (BUILD_APP_CFEM)
    add_subdirectory(cfem)
endif()

option(BUILD_APP_FINITE_STRAINS "Build Finite Strains application" ON)
if (BUILD_APP_FINITE_STRAINS)
    add_subdirectory(finite_strains)
endif()

option(BUILD_UTILS "Build utilities" ON)
if (BUILD_UTILS)
    add_subdirectory(utils)
endif()

option(BUILD_APP_BARYTEST "Build barycenter test application" ON)
if (BUILD_APP_BARYTEST)
    add_subdirectory(barycenter_test)
endif()

