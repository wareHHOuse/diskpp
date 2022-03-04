include(FindPackageHandleStandardArgs)

find_path(GMSH_INCLUDE_DIR
    NAMES gmsh.h
    HINTS ENV GMSH_ROOT ${GMSH_ROOT}
    PATH_SUFFIXES include)

find_library(GMSH_LIBRARIES
    NAMES gmsh
    HINTS ENV GMSH_ROOT ${GMSH_ROOT}
    PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH_LIBRARIES GMSH_INCLUDE_DIR)

if (GMSH_FOUND)
    add_library(GMSH INTERFACE)
    target_link_libraries(GMSH INTERFACE "${GMSH_LIBRARIES}")
    target_include_directories(GMSH INTERFACE "${GMSH_INCLUDE_DIR}")
    target_compile_definitions(GMSH INTERFACE -DHAVE_GMSH)
endif()
