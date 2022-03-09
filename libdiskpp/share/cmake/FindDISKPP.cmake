include(FindPackageHandleStandardArgs)

find_path(DISKPP_INCLUDE_DIR
    NAMES diskpp_signature.hpp
    HINTS ENV DISKPP_ROOT ${DISKPP_ROOT}
    PATH_SUFFIXES include)

find_library(DISKPP_LIBRARIES
    NAMES diskpp
    HINTS ENV DISKPP_ROOT ${DISKPP_ROOT}
    PATH_SUFFIXES lib)

find_package_handle_standard_args(DISKPP DEFAULT_MSG DISKPP_LIBRARIES DISKPP_INCLUDE_DIR)

if (DISKPP_FOUND)
    add_library(DISKPP INTERFACE)
    target_link_libraries(DISKPP INTERFACE ${DISKPP_LIBRARIES})
    target_include_directories(DISKPP INTERFACE "${DISKPP_INCLUDE_DIR}")
endif()
