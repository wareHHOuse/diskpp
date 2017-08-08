
enable_language (Fortran)

set (AGMG_DIR "${PROJECT_SOURCE_DIR}/contrib/agmg/")
set (AGMG_SOURCE_DIR "${AGMG_DIR}/agmg_source/")

if ( EXISTS "${AGMG_SOURCE_DIR}" AND IS_DIRECTORY "${AGMG_SOURCE_DIR}" )
    set (HAVE_AGMG TRUE)

    get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

    if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
        set (CMAKE_Fortran_FLAGS " -w")
    endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

    set (COMMON_SOURCES "${AGMG_SOURCE_DIR}/blas_agmg.f" "${AGMG_SOURCE_DIR}/lapack_agmg.f")
    #set (COMMON_SOURCES "")
    set (SAGMG_SOURCES  "${AGMG_SOURCE_DIR}/sagmg.f90" "${AGMG_SOURCE_DIR}/sagmg_mumps.f90")
    set (DAGMG_SOURCES  "${AGMG_SOURCE_DIR}/dagmg.f90" "${AGMG_SOURCE_DIR}/dagmg_mumps.f90")
    set (CAGMG_SOURCES  "${AGMG_SOURCE_DIR}/cagmg.f90" "${AGMG_SOURCE_DIR}/cagmg_mumps.f90")
    set (ZAGMG_SOURCES  "${AGMG_SOURCE_DIR}/zagmg.f90" "${AGMG_SOURCE_DIR}/zagmg_mumps.f90")

    add_library (sagmg SHARED ${SAGMG_SOURCES} ${COMMON_SOURCES})
    add_library (dagmg SHARED ${DAGMG_SOURCES} ${COMMON_SOURCES})
    add_library (cagmg SHARED ${CAGMG_SOURCES} ${COMMON_SOURCES})
    add_library (zagmg SHARED ${ZAGMG_SOURCES} ${COMMON_SOURCES})

    if ( CYGWIN )
        set_target_properties (sagmg PROPERTIES LINK_FLAGS
            "-Wl,--version-script=\"${AGMG_DIR}/libsagmg.version\"")

        set_target_properties (dagmg PROPERTIES LINK_FLAGS
            "-Wl,--version-script=\"${AGMG_DIR}/libdagmg.version\"")

        set_target_properties (cagmg PROPERTIES LINK_FLAGS
            "-Wl,--version-script=\"${AGMG_DIR}/libcagmg.version\"")

        set_target_properties (zagmg PROPERTIES LINK_FLAGS
            "-Wl,--version-script=\"${AGMG_DIR}/libzagmg.version\"")
    endif (CYGWIN)

    if ( APPLE )
        set_target_properties (sagmg PROPERTIES LINK_FLAGS "-Wl,-exported_symbols_list \
                               \"${AGMG_DIR}/sagmg_export.map\"")
                               #target_link_libraries(sagmg "-lblas -llapack")

        set_target_properties (dagmg PROPERTIES LINK_FLAGS "-Wl,-exported_symbols_list \
                               \"${AGMG_DIR}/dagmg_export.map\"")
                               #target_link_libraries(dagmg "-lblas -llapack")

        set_target_properties (cagmg PROPERTIES LINK_FLAGS "-Wl,-exported_symbols_list \
                               \"${AGMG_DIR}/cagmg_export.map\"")
                               #target_link_libraries(cagmg "-lblas -llapack")

        set_target_properties (zagmg PROPERTIES LINK_FLAGS "-Wl,-exported_symbols_list \
                               \"${AGMG_DIR}/zagmg_export.map\"")
                               #target_link_libraries(zagmg "-lblas -llapack")

    endif ( APPLE )

    set(LINK_LIBS ${LINK_LIBS} sagmg dagmg cagmg zagmg)

endif ()
