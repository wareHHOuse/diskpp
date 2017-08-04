
add_subdirectory(agmg)

set(diskpp_source ${diskpp_source}
    ${CMAKE_CURRENT_SOURCE_DIR}/simpleps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecounter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapped_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colormanip.cpp
    PARENT_SCOPE)

add_subdirectory(jburkardt)
set(LINK_LIBS ${LINK_LIBS} ${SOLVER_WRAPPERS_LIBRARIES} PARENT_SCOPE)
