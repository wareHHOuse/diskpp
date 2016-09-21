set(quadratures_source
    ${CMAKE_CURRENT_SOURCE_DIR}/triangle_dunavant_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/jb_common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simplex_gm_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tetrahedron_arbq_rule.cpp
    )

add_library(quadratures SHARED ${quadratures_source})
set(LINK_LIBS "${LINK_LIBS}" quadratures PARENT_SCOPE)
