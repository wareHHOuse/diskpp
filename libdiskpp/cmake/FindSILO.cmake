include(FindPackageHandleStandardArgs)

find_path(SILO_INCLUDE_DIR
    NAMES silo.h
    HINTS ENV SILO_ROOT ${SILO_ROOT}
    PATH_SUFFIXES include)

find_library(SILO_LIBRARY
    NAMES silo siloh5
    HINTS ENV SILO_ROOT ${SILO_ROOT}
    PATH_SUFFIXES lib)

find_package_handle_standard_args(SILO DEFAULT_MSG SILO_LIBRARY SILO_INCLUDE_DIR)

if(SILO_FOUND AND NOT TARGET SILO::SILO)
    add_library(SILO::SILO UNKNOWN IMPORTED)
    set_target_properties(SILO::SILO PROPERTIES
        IMPORTED_LOCATION "${SILO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SILO_INCLUDE_DIR}"
    )
endif()