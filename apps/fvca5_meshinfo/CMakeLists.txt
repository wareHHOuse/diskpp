set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt4)
if (Qt4_FOUND)
    add_executable(fvca5_meshinfo fvca5_meshinfo.cpp MeshWidget.cpp)
    target_link_libraries(fvca5_meshinfo ${LINK_LIBS} Qt4::QtGui)
endif()
