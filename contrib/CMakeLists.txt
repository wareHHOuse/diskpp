

add_subdirectory(jburkardt)
set(diskpp_source ${diskpp_source} ${quadrature_source})

add_subdirectory(triangle)
set(diskpp_source ${diskpp_source} ${triangle_source})

set(diskpp_source ${diskpp_source} 
    #${CMAKE_CURRENT_SOURCE_DIR}/simpleps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecounter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapped_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colormanip.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gmsh_tools/gmshElement.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gmsh_tools/gmshMesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gmsh_tools/gmshData.cpp
    PARENT_SCOPE)
