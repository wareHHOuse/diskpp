

add_subdirectory(jburkardt)
set(diskpp_source ${diskpp_source} ${quadratures_source}
    ${CMAKE_CURRENT_SOURCE_DIR}/simpleps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecounter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapped_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/colormanip.cpp
    PARENT_SCOPE)


