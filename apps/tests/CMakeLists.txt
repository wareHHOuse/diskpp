add_executable(scalar_reconstruction scalar_reconstruction.cpp)
target_link_libraries(scalar_reconstruction ${LINK_LIBS})

add_executable(scalar_stabilization scalar_stabilization.cpp)
target_link_libraries(scalar_stabilization ${LINK_LIBS})

add_executable(vector_reconstruction vector_reconstruction.cpp)
target_link_libraries(vector_reconstruction ${LINK_LIBS})

add_executable(vector_stabilization vector_stabilization.cpp)
target_link_libraries(vector_stabilization ${LINK_LIBS})