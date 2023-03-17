cmake_minimum_required(VERSION 3.18)
project(libdiskpp)

######################################################################
## Required C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

######################################################################
## RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

######################################################################
## Find Eigen 3
#find_package(Eigen3 REQUIRED)
#set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

include_directories(include)

add_subdirectory(contrib/jburkardt)
set(LINK_LIBS ${LINK_LIBS} JBURKARDT)

add_subdirectory(contrib/solver_wrappers)
set(LINK_LIBS ${LINK_LIBS} SOLVER_WRAPPERS)

add_subdirectory(contrib/matplotplusplus)
set(LINK_LIBS ${LINK_LIBS} matplot)

add_library(diskpp SHARED src/common/util.cpp src/common/mapped_file.cpp src/common/colormanip.cpp)
target_link_libraries(diskpp ${LINK_LIBS})
target_include_directories(diskpp PUBLIC include)

install(TARGETS diskpp DESTINATION lib)
install(DIRECTORY include/diskpp DESTINATION include)
