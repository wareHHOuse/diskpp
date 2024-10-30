set(LINK_LIBS diskpp)

add_executable(hhodd src/hhodd.cpp)
target_link_libraries(hhodd ${LINK_LIBS})
install(TARGETS hhodd RUNTIME DESTINATION bin)

set(HHODD_LOCAL_SOLVER "Eigen::SparseLU" CACHE STRING "HHO-DD local solver")
set_property(CACHE HHODD_LOCAL_SOLVER PROPERTY STRINGS "Eigen::SparseLU" "PARDISO" "MUMPS")

if (HHODD_LOCAL_SOLVER STREQUAL "Eigen::SparseLU")
    target_compile_definitions(hhodd PUBLIC -DHHODD_USE_EIGEN_SPARSELU)
endif()

if (HHODD_LOCAL_SOLVER STREQUAL "PARDISO")
    target_compile_definitions(hhodd PUBLIC -DHHODD_USE_PARDISO)
endif()

if (HHODD_LOCAL_SOLVER STREQUAL "MUMPS")
    target_compile_definitions(hhodd PUBLIC -DHHODD_USE_MUMPS)
endif()

install(DIRECTORY share/ DESTINATION share)
