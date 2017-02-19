
set(diskpp_source ${diskpp_source}
    ${CMAKE_CURRENT_SOURCE_DIR}/gmshElement.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gmshMesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gmshData.cpp
    PARENT_SCOPE)


set(LINK_LIBS "${LINK_LIBS}" PARENT_SCOPE)
