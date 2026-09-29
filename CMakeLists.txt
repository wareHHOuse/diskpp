cmake_minimum_required(VERSION 3.30)
project(diskpp)

######################################################################
## Address Sanitizer
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

if (ENABLE_ASAN)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g -O1)

    add_link_options(-fsanitize=address,undefined)
endif()

######################################################################
## Compiler reporting: disable deprecation warnings
option(OPT_DISABLE_DEPRECATION_WARNINGS "Disable deprecation warnings" ON)
if(OPT_DISABLE_DEPRECATION_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

add_subdirectory(libdiskpp)
add_subdirectory(apps)
add_subdirectory(unit_tests)


