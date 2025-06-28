add_executable(elastodynamics elastodynamics.cpp ${fitted_waves_headers} ${fitted_waves_sources})
target_link_libraries(elastodynamics ${LINK_LIBS})
