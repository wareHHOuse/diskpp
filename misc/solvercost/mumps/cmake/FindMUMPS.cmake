

#find_path (MUMPS_DIR include/mumps_compat.h HINTS ENV MUMPS_DIR PATHS $ENV{HOME}/mumps DOC "Mumps Directory")

include(FindPackageHandleStandardArgs)

find_path(MUMPS_INCLUDE_DIRS
          NAMES mumps_compat.h
          PATHS /usr/include
		/usr/local/include)

find_library(MUMPS_smumps_LIBRARY
             NAMES smumps
             PATHS /usr/lib
		   /usr/local/lib)

find_library(MUMPS_dmumps_LIBRARY
             NAMES dmumps
             PATHS /usr/lib
		   /usr/local/lib)

find_library(MUMPS_cmumps_LIBRARY
             NAMES cmumps
             PATHS /usr/lib
		   /usr/local/lib)

find_library(MUMPS_zmumps_LIBRARY
             NAMES zmumps
             PATHS /usr/lib
		   /usr/local/lib)

find_library(MUMPS_mumps_common_LIBRARY
             NAMES mumps_common
             PATHS /usr/lib
                   /usr/local/lib)

find_library(MUMPS_pord_LIBRARY
             NAMES pord
             PATHS /usr/lib
             /usr/local/lib)

if ( MUMPS_smumps_LIBRARY AND
     MUMPS_dmumps_LIBRARY AND 
     MUMPS_cmumps_LIBRARY AND
     MUMPS_zmumps_LIBRARY AND
     MUMPS_mumps_common_LIBRARY AND
     MUMPS_pord_LIBRARY)
 set(MUMPS_LIBRARIES ${MUMPS_smumps_LIBRARY} ${MUMPS_dmumps_LIBRARY}  ${MUMPS_cmumps_LIBRARY} ${MUMPS_zmumps_LIBRARY} ${MUMPS_mumps_common_LIBRARY} ${MUMPS_pord_LIBRARY})
endif ()
 
#set (MUMPS_DIR "/usr/local/Cellar/mumps/4.10.0/")

#IF(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)
#  SET(MUMPS_FOUND YES)
#  SET(MUMPS_INCLUDES ${MUMPS_DIR})
#  find_path (MUMPS_INCLUDE_DIR mumps_compat.h HINTS "${MUMPS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
#  list(APPEND MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})
#  FILE(GLOB MUMPS_LIBRARIES RELATIVE "${MUMPS_DIR}/lib" "${MUMPS_DIR}/lib/lib*.dylib")
#ELSE(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)
#  SET(MUMPS_FOUND NO)
#ENDIF(EXISTS ${MUMPS_DIR}/include/mumps_compat.h)

find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES MUMPS_INCLUDE_DIRS)

