#add_executable(signorini signorini.cpp)
#target_link_libraries(signorini ${LINK_LIBS})

add_executable(signorini_cfem signorini_fem_fix_point.cpp)
target_link_libraries(signorini_cfem ${LINK_LIBS})
