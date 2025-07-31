option(BUILD_CONDITIONING_TESTS "Build conditioning tests" ON)
if (BUILD_CONDITIONING_TESTS)
    add_subdirectory(conditioning)
endif()

option(BUILD_APP_VEMASM "Build VEM assembly testing app" ON)
if (BUILD_APP_VEMASM)
    add_subdirectory(vemasm)
endif()

option(BUILD_DISKPP_TUTORIALS "Build DiSk++ tutorials" ON)
if (BUILD_DISKPP_TUTORIALS)
    add_subdirectory(diskpp_tutorials)
endif()

option(BUILD_APP_STEKLOV_EIGS "Build Steklov eigsolver" ON)
if (BUILD_APP_STEKLOV_EIGS)
    add_subdirectory(steklov_eigs)
endif()

option(BUILD_APP_DIFFUSION "Build Diffusion applications" ON)
if (BUILD_APP_DIFFUSION)
    add_subdirectory(diffusion)
endif()

option(BUILD_APP_POISSON "Build Poisson equation solver" ON)
if (BUILD_APP_POISSON)
    add_subdirectory(poisson_solver)
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

option(BUILD_APP_NSM "Build Non-linear Solid Mechanics applications" ON)
if (BUILD_APP_NSM)
    add_subdirectory(nonlinear_solid_mechanics)
endif()

option(BUILD_APP_UNST_LAPLACIAN "Build unsteady laplacian application" ON)
if (BUILD_APP_UNST_LAPLACIAN)
    add_subdirectory(unsteady_laplacian)
endif()

option(BUILD_APP_VECTOR_LAPLACIAN "Build vector laplacian application" ON)
if (BUILD_APP_VECTOR_LAPLACIAN)
    add_subdirectory(vector_laplacian)
endif()

option(BUILD_UTILS "Build utilities" ON)
if (BUILD_UTILS)
    add_subdirectory(utils)
endif()

option(BUILD_APP_BARYTEST "Build barycenter test application" ON)
if (BUILD_APP_BARYTEST)
    add_subdirectory(barycenter_test)
endif()

option(BUILD_APP_VEM_2D "Build VEM 2d application" ON)
if (BUILD_APP_VEM_2D)
    add_subdirectory(vem2d)
endif()

option(BUILD_APP_NITSCHE "Build Nitsche application" ON)
if (BUILD_APP_NITSCHE)
    add_subdirectory(nitsche)
endif()

