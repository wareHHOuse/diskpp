include(FindPackageHandleStandardArgs)


find_path(SILO_INCLUDE_DIRS
          NAMES silo.h
          PATHS /usr/include /usr/local/include)

find_library(SILO_LIBRARIES
             NAMES silo siloh5
             PATHS /usr/lib /usr/local/include)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SILO DEFAULT_MSG SILO_LIBRARIES SILO_INCLUDE_DIRS)
