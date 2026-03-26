include(FindPackageHandleStandardArgs)

find_path(GMSH_INCLUDE_DIR
    NAMES gmsh.h
    HINTS ENV GMSH_ROOT ${GMSH_ROOT}
    PATH_SUFFIXES include)

find_library(GMSH_LIBRARY
    NAMES gmsh
    HINTS ENV GMSH_ROOT ${GMSH_ROOT}
    PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH_LIBRARY GMSH_INCLUDE_DIR)

if(GMSH_FOUND AND NOT TARGET GMSH::GMSH)
    add_library(GMSH::GMSH INTERFACE IMPORTED)
    set_target_properties(GMSH::GMSH PROPERTIES
        IMPORTED_LOCATION "${GMSH_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
    )
endif()
