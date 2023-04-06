include(FindPackageHandleStandardArgs)

find_path(SILO_INCLUDE_DIR
    NAMES silo.h
    HINTS ENV SILO_ROOT ${SILO_ROOT}
    PATH_SUFFIXES include)

find_library(SILO_LIBRARIES
    NAMES silo siloh5
    HINTS ENV SILO_ROOT ${SILO_ROOT}
    PATH_SUFFIXES lib)

find_package_handle_standard_args(SILO DEFAULT_MSG SILO_LIBRARIES SILO_INCLUDE_DIR)

if (SILO_FOUND)
    add_library(SILO INTERFACE)
    target_link_libraries(SILO INTERFACE ${SILO_LIBRARIES})
    target_include_directories(SILO INTERFACE "${SILO_INCLUDE_DIR}")
    target_compile_definitions(SILO INTERFACE -DHAVE_SILO)
endif()