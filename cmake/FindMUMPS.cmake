include(FindPackageHandleStandardArgs)

if(NOT MUMPS_SEQUENTIAL_OR_PARALLEL)
  set(MUMPS_FIND_COMPONENTS SEQ)
endif()

set(MUMPS_SEQUENTIAL_OR_PARALLEL "Sequential"
    CACHE STRING "Find sequential or parallel MUMPS.")

set_property(CACHE MUMPS_SEQUENTIAL_OR_PARALLEL
    PROPERTY STRINGS "Sequential" "Parallel")

if(${MUMPS_SEQUENTIAL_OR_PARALLEL} STREQUAL "Sequential")
    set(MUMPSS_NAMES smumps_seq smumps)
    set(MUMPSD_NAMES dmumps_seq dmumps)
    set(MUMPSC_NAMES cmumps_seq cmumps)
    set(MUMPSZ_NAMES zmumps_seq zmumps)
    set(MUMPS_COMMON_NAMES mumps_common_seq mumps_common)
elseif(${MUMPS_SEQUENTIAL_OR_PARALLEL} STREQUAL "Parallel")
  set(MUMPSS_NAMES smumps_ptscotch smumps_scotch smumps)
  set(MUMPSD_NAMES dmumps_ptscotch dmumps_scotch dmumps)
  set(MUMPSC_NAMES cmumps_ptscotch cmumps_scotch cmumps)
  set(MUMPSZ_NAMES zmumps_ptscotch zmumps_scotch zmumps)
  set(MUMPS_COMMON_NAMES mumps_common_ptscotch mumps_common_scotch mumps_common)
else()
  message(FATAL_ERROR "Invalid value \"${MUMPS_FIND_COMPONENTS}\" for MUMPS_SEQUENTIAL_OR_PARALLEL.")
endif()

######################################################################
## Libraries common to all MUMPS versions
find_library(MUMPS_COMMON_LIBRARY
    NAMES ${MUMPS_COMMON_NAMES}
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES lib lib64)

if(MUMPS_COMMON_LIBRARY)
    add_library(MUMPS::COMMON UNKNOWN IMPORTED)
    set_target_properties(MUMPS::COMMON PROPERTIES
        IMPORTED_LOCATION "${MUMPS_COMMON_LIBRARY}"
    )
    set(MUMPS_COMMON_FOUND TRUE)
else()
    set(MUMPS_COMMON_FOUND FALSE)
endif()



######################################################################
## mpiseq
if(${MUMPS_SEQUENTIAL_OR_PARALLEL} STREQUAL "Sequential")
    find_library(MUMPS_MPISEQ_LIBRARY
        NAMES mpiseq
        HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
        PATH_SUFFIXES lib lib64)
    
    if(MUMPS_MPISEQ_LIBRARY)
        add_library(MUMPS::MPISEQ UNKNOWN IMPORTED)
        set_target_properties(MUMPS::MPISEQ PROPERTIES
            IMPORTED_LOCATION "${MUMPS_MPISEQ_LIBRARY}"
        )
        set(MUMPS_MPISEQ_FOUND TRUE)
        message(STATUS "Found MUMPS mpiseq")
    else()
        set(MUMPS_MPISEQ_FOUND FALSE)
    endif()
endif()

######################################################################
## Extras
set(extras_libs ptscotch scotch metis pord pord_seq)
foreach(extra IN LISTS extras_libs)
    find_library(MUMPS_${extra}_LIBRARY
        NAMES ${extra}
        HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
        PATH_SUFFIXES lib lib64)
    if(MUMPS_${extra}_LIBRARY)
        message(STATUS "MUMPS extra: " ${extra})
        list(APPEND ${MUMPS_LIBRARIES} ${MUMPS_${extra}_LIBRARY})
    endif()
endforeach()

######################################################################
## Real single precision MUMPS
find_path(MUMPSS_INCLUDE_DIR
    NAMES smumps_c.h
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES include
  )

find_library(MUMPSS_LIBRARY
    NAMES ${MUMPSS_NAMES}
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES lib lib64)

if (MUMPSS_INCLUDE_DIR AND MUMPSS_LIBRARY)
    add_library(MUMPS::S UNKNOWN IMPORTED)
    set_target_properties(MUMPS::S PROPERTIES
        IMPORTED_LOCATION "${MUMPSS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPSS_INCLUDE_DIR}"
    )
    set(MUMPSS_FOUND TRUE)
    set(FOUND_SOME_MUMPS TRUE)
    message(STATUS "Found MUMPS real single precision")
else()
    set(MUMPSS_FOUND FALSE)
endif()

######################################################################
## Real double precision MUMPS
find_path(MUMPSD_INCLUDE_DIR
    NAMES dmumps_c.h
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES include
  )

find_library(MUMPSD_LIBRARY
    NAMES ${MUMPSD_NAMES}
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES lib lib64)

if (MUMPSD_INCLUDE_DIR AND MUMPSD_LIBRARY)
    add_library(MUMPS::D UNKNOWN IMPORTED)
    set_target_properties(MUMPS::D PROPERTIES
        IMPORTED_LOCATION "${MUMPSD_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPSD_INCLUDE_DIR}"
    )
    set(MUMPSD_FOUND TRUE)
    set(FOUND_SOME_MUMPS TRUE)
    message(STATUS "Found MUMPS real double precision")
else()
    set(MUMPSD_FOUND FALSE)
endif()

######################################################################
## Complex single precision MUMPS
find_path(MUMPSC_INCLUDE_DIR
    NAMES cmumps_c.h
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES include
  )

find_library(MUMPSC_LIBRARY
    NAMES ${MUMPSC_NAMES}
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES lib lib64)

if (MUMPSC_INCLUDE_DIR AND MUMPSC_LIBRARY)
    add_library(MUMPS::C UNKNOWN IMPORTED)
    set_target_properties(MUMPS::C PROPERTIES
        IMPORTED_LOCATION "${MUMPSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPSC_INCLUDE_DIR}"
    )
    set(MUMPSC_FOUND TRUE)
    set(FOUND_SOME_MUMPS TRUE)
    message(STATUS "Found MUMPS complex single precision")
else()
    set(MUMPSC_FOUND FALSE)
endif()

######################################################################
## Complex double precision MUMPS
find_path(MUMPSZ_INCLUDE_DIR
    NAMES zmumps_c.h
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES include
  )

find_library(MUMPSZ_LIBRARY
    NAMES ${MUMPSZ_NAMES}
    HINTS ENV MUMPS_ROOT ${MUMPS_ROOT}
    PATH_SUFFIXES lib lib64)

if (MUMPSZ_INCLUDE_DIR AND MUMPSZ_LIBRARY)
    add_library(MUMPS::Z UNKNOWN IMPORTED)
    set_target_properties(MUMPS::Z PROPERTIES
        IMPORTED_LOCATION "${MUMPSZ_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${MUMPSZ_INCLUDE_DIR}"
    )
    set(MUMPS_Z_FOUND TRUE)
    set(FOUND_SOME_MUMPS TRUE)
    message(STATUS "Found MUMPS complex double precision")
else()
    set(MUMPS_Z_FOUND FALSE)
endif()

find_package_handle_standard_args(MUMPS
    REQUIRED_VARS
        MUMPS_COMMON_LIBRARY
    HANDLE_COMPONENTS
)

if (MUMPS_FOUND)
    message(STATUS "xxx")
    add_library(MUMPS::MUMPS INTERFACE IMPORTED)

    target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::COMMON)

    if (MUMPS_MPISEQ_FOUND)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::MPISEQ)
    endif()

    if (MUMPSS_FOUND)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::S)
    endif()

    if (MUMPSD_FOUND)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::D)
    endif()

    if (MUMPSC_FOUND)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::C)
    endif()

    if (MUMPS_Z_FOUND)
        target_link_libraries(MUMPS::MUMPS INTERFACE MUMPS::Z)
    endif()

endif()

######################################################################
## Interface
if (FOUND_SOME_MUMPS)
    add_library(MUMPS INTERFACE)
    target_link_libraries(MUMPS INTERFACE ${MUMPS_VERSIONS})
    target_include_directories(MUMPS INTERFACE include)

    list(LENGTH MUMPS_VERSIONS num_vers)
    if(num_vers GREATER 0)
        target_compile_definitions(MUMPS INTERFACE -DHAVE_MUMPS)
    endif()

    
endif()

