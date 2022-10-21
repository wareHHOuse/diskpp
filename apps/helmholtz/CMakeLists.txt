set(LINK_LIBS diskpp)

add_executable(helmholtz helmholtz.cpp)
target_link_libraries(helmholtz ${LINK_LIBS})
set_target_properties(helmholtz PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)


