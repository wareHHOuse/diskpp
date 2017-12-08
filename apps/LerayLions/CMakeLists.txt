add_executable(LerayLions LerayLions.cpp)
target_link_libraries(LerayLions ${LINK_LIBS})
install(TARGETS LerayLions RUNTIME DESTINATION bin)


#add_executable(LerayLions_test_convergence_rates LerayLions_test_convergence_rates.cpp)
#target_link_libraries(LerayLions_test_convergence_rates ${LINK_LIBS})
#install(TARGETS LerayLions_test_convergence_rates RUNTIME DESTINATION bin)