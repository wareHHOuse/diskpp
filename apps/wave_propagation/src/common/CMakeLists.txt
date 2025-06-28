
set(fitted_waves_sources ${fitted_waves_sources} 
	#${CMAKE_CURRENT_SOURCE_DIR}/fitted_geometry_builders.cpp
PARENT_SCOPE)

set(fitted_waves_headers ${fitted_waves_headers} 
	${CMAKE_CURRENT_SOURCE_DIR}/display_settings.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/preprocessor.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/postprocessor.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/acoustic_material_data.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/elastic_material_data.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/scal_analytic_functions.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/vec_analytic_functions.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/scal_vec_analytic_functions.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/fitted_geometry_builders.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/assembly_index.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/acoustic_one_field_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/acoustic_two_fields_assembler.hpp	
	${CMAKE_CURRENT_SOURCE_DIR}/elastodynamic_one_field_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/elastodynamic_two_fields_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/elastodynamic_three_fields_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/elastoacoustic_two_fields_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/elastoacoustic_four_fields_assembler.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/fitted_geometry_builders.hpp	
	${CMAKE_CURRENT_SOURCE_DIR}/linear_solver.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/dirk_hho_scheme.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/erk_hho_scheme.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/ssprk_hho_scheme.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/dirk_butcher_tableau.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/erk_butcher_tableau.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/ssprk_shu_osher_tableau.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/lsrk_butcher_tableau.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/sourceterm.hpp
PARENT_SCOPE)

