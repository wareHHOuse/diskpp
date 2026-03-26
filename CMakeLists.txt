cmake_minimum_required(VERSION 3.30)
project(diskpp)

######################################################################
## Compiler reporting: disable deprecation warnings
option(OPT_DISABLE_DEPRECATION_WARNINGS "Disable deprecation warnings" ON)
if(OPT_DISABLE_DEPRECATION_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

add_subdirectory(libdiskpp)
add_subdirectory(apps)


