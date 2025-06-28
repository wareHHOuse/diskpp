

//  Created by Romain Mottier

// ../../elastoacoustic -k 4 -s 0 -r 0 -c 1 -p 1 -l 9 -n 2000 -f 0 -e 0

void EHHOFirstOrder_conv_tests(int argc, char **argv);

void EHHOFirstOrder_conv_tests(int argc, char **argv){
  
    // ##################################################
    // ################################################## Simulation paramaters 
    // ##################################################
    
    std::cout << std::endl << bold << red << "   EXPLICIT COUPLING CONV TEST" << reset << std::endl;
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // ##################################################
    // ################################################## HHO setting 
    // ##################################################

    for (size_t k = 1; k <= sim_data.m_k_degree; k++) {

        std::cout << std::endl << bold << red << "   Polynomial degree k : " << k << reset << std::endl;

        // Creating HHO approximation spaces and corresponding linear operator
        size_t cell_k_degree = k;
        if (sim_data.m_hdg_stabilization_Q)
            cell_k_degree++;
        disk::hho_degree_info hho_di(cell_k_degree, k);
    
        // ##################################################
        // ################################################## Loop over level of space refinement 
        // ##################################################

        for (size_t l = 0; l <= sim_data.m_n_divs; l++) {

            std::cout << bold << cyan << "      Space refinment level -l : " << l << reset << std::endl;

            // ##################################################
            // ################################################## Mesh generation 
            // ##################################################
                       
            typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
            typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
            typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
            mesh_type msh;
            
            if (sim_data.m_polygonal_mesh_Q) {
                // Mesh availability
                auto validate_l = [](size_t l) -> size_t {
                    if (!((0 <= l) && (l < 15))) {
                        std::cout << std::endl << std::endl;
                        std::cout << bold << red << "Warning:: Only few polygonal meshes available.";
                        std::cout << std::endl << std::endl;
                        return 4;
                    }
                };
                // size_t l = validate_l(sim_data.m_n_divs);
                polygon_2d_mesh_reader<RealType> mesh_builder;
                std::vector<std::string> mesh_files;
                // Simplicial meshes
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l0_conv_test_1.0.txt");    // l = 0
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l1_conv_test_0.35.txt");   // 
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l2_conv_test_0.15.txt");   // l = 2
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l3_conv_test_0.07.txt");   // 
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l4_conv_test_0.035.txt");  // l = 4
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l5_conv_test_0.026.txt");  // 
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l6_conv_test_0.017.txt");  // l = 6
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l7_conv_test_0.0125.txt"); // 
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l8_conv_test_0.0085.txt"); // l = 8
                // mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l9_conv_test_0.005.txt");  // 
                // Polyhedral meshes
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_32.txt");    // -l 0
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_64.txt");    // -l 1
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_128.txt");   // -l 2
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_256.txt");   // -l 3
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_512.txt");   // -l 4
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_1024.txt");  // -l 5 
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_2048.txt");  // -l 6
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_4096.txt");  // -l 7 
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_8192.txt");  // -l 8
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_16384.txt"); // -l 9
                mesh_files.push_back("../../meshes/conv_test/poly_bad_quality/poly_32768.txt"); // -l 10
                // Reading the mesh
                mesh_builder.set_poly_mesh_file(mesh_files[l]);
                mesh_builder.build_mesh();
                mesh_builder.move_to_mesh_storage(msh);
                mesh_builder.remove_duplicate_points();
            }
            else {
                RealType lx = 2.0;  
                RealType ly = 1.0;          
                size_t nx = 4;
                size_t ny = 2;
                cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
                mesh_builder.refine_mesh(sim_data.m_n_divs);
                mesh_builder.set_translation_data(-1.0, 0.0);
                mesh_builder.build_mesh();
                mesh_builder.move_to_mesh_storage(msh);
            }
            RealType h = 10;
            for (auto & cell : msh ) {
                auto cell_ind = msh.lookup(cell);
                mesh_type::point_type bar = barycenter(msh, cell);
                RealType h_l = diameter(msh, cell);
                if (h_l < h) 
                    h = h_l;
            }
            
            // ##################################################
            // ################################################## Time discretization infos 
            // ##################################################
            size_t nt = 10;
            for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) 
                nt  = sim_data.m_nt_divs;  
            RealType ti = 0.0;
            RealType tf = 1.0;
            RealType dt = (tf-ti)/nt;
            RealType t = ti;

            // ERK schemes
            Matrix<RealType, Dynamic, Dynamic> a;
            Matrix<RealType, Dynamic, 1> b;
            Matrix<RealType, Dynamic, 1> c;            
            int s = 4;
            erk_butcher_tableau::erk_tables(s, a, b, c);

            // ##################################################
            // ################################################## Manufactured solution 
            // ##################################################
            
            scal_vec_analytic_functions functions;
            // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
            functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInTime);
            // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
            // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial_paper);
            
            // Elastic analytical functions
            auto u_fun    = functions.Evaluate_u(t);
            auto v_fun    = functions.Evaluate_v(t);
            auto a_fun    = functions.Evaluate_a(t);
            auto f_fun    = functions.Evaluate_f(t);
            auto flux_fun = functions.Evaluate_sigma(t);
            
            // Acoustic analytical functions
            auto s_u_fun    = functions.Evaluate_s_u(t);
            auto s_v_fun    = functions.Evaluate_s_v(t);
            auto s_a_fun    = functions.Evaluate_s_a(t);
            auto s_f_fun    = functions.Evaluate_s_f(t);
            auto s_flux_fun = functions.Evaluate_s_q(t);
            
            // ##################################################
            // ################################################## Material data 
            // ##################################################

            // Classify cells per material data and bc faces
            auto elastic_mat_fun = [](const typename mesh_type::point_type& pt) -> elastic_material_data<RealType> {
                double x,y;
                x = pt.x();
                y = pt.y();
                RealType rho, vp, vs;
                rho = 1.0;            // Solid mass density
                vp  = std::sqrt(3.0); // Seismic compressional velocity vp
                vs  = 1.0;            // Seismic shear velocity vs
                elastic_material_data<RealType> material(rho,vp,vs);
                return material;
            };
            auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
                double x,y;
                x = pt.x();
                y = pt.y();
                RealType rho, vp;
                rho = 1.0; // Fluid mass density
                vp  = 1.0; // Seismic compressional velocity vp
                acoustic_material_data<RealType> material(rho,vp);
                return material;
            };
            
            // ##################################################
            // ################################################## Structure setting
            // ##################################################
            
            std::map<size_t,elastic_material_data<RealType>> e_material;
            std::map<size_t,acoustic_material_data<RealType>> a_material;
            std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
            std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
            RealType eps = 1.0e-10;
            for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
                const auto face = *face_it;
                mesh_type::point_type bar = barycenter(msh, face);
                auto fc_id = msh.lookup(face);
                if (std::fabs(bar.x()) < eps) {
                    interface_face_indexes.insert(fc_id);
                    continue;
                } 
            }
            for (auto & cell : msh ) {
                auto cell_ind = msh.lookup(cell);
                mesh_type::point_type bar = barycenter(msh, cell);
                // Assigning the material properties
                if (bar.x() > 0) {
                    acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
                    a_material.insert(std::make_pair(cell_ind,material));
                }
                else {
                    elastic_material_data<RealType> material = elastic_mat_fun(bar);
                    e_material.insert(std::make_pair(cell_ind,material));
                }
                // Detection of faces on the interfaces
                auto cell_faces = faces(msh,cell);
                for (auto face :cell_faces) {
                    auto fc_id = msh.lookup(face);
                    bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
                    if (is_member_Q) {
                        if (bar.x() > 0) 
                            interface_cell_pair_indexes[fc_id].second = cell_ind;
                        else 
                            interface_cell_pair_indexes[fc_id].first = cell_ind;
                    }
                }
            }            
            // Internal faces structure 
            std::set<size_t> elastic_internal_faces;
            std::set<size_t> acoustic_internal_faces;
            for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
                const auto face = *face_it;
                mesh_type::point_type bar = barycenter(msh, face);
                auto fc_id = msh.lookup(face);      
                bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
                if (!is_member_Q) {
                    if (bar.y() > 0) 
                        acoustic_internal_faces.insert(fc_id);    
                    else 
                        elastic_internal_faces.insert(fc_id);
                }
            }
            size_t bc_elastic_id  = 0;
            size_t bc_acoustic_id = 1;
            for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
                auto face = *face_it;
                mesh_type::point_type bar = barycenter(msh, face);
                auto fc_id = msh.lookup(face);
                if (bar.x() > 0) {
                    disk::boundary_descriptor bi{bc_acoustic_id, true};
                    msh.backend_storage()->boundary_info.at(fc_id) = bi;
                    acoustic_bc_face_indexes.insert(fc_id);
                }
                else {
                    disk::boundary_descriptor bi{bc_elastic_id, true};
                    msh.backend_storage()->boundary_info.at(fc_id) = bi;
                    elastic_bc_face_indexes.insert(fc_id);
                }   
            }
            e_boundary_type e_bnd(msh);
            a_boundary_type a_bnd(msh);
            e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
            a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);
            
            // ##################################################
            // ################################################## Solving a primal HHO mixed problem 
            // ##################################################
            
            auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
            assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
            assembler.set_coupling_stabilization();
            if (sim_data.m_scaled_stabilization_Q)
                assembler.set_scaled_stabilization();
            assembler.assemble_mass(msh);
            assembler.assemble_coupling_terms(msh);
            
            // ##################################################
            // ################################################## Projecting initial data 
            // ##################################################
            
            Matrix<RealType, Dynamic, 1> x_dof;
            assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
            assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);
            
            // ##################################################
            // ################################################## Solving a primal HHO mixed problem 
            // ##################################################
            
            assembler.assemble(msh, f_fun, s_f_fun, true);
            assembler.LHS += assembler.COUPLING; 
            size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
            size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
            size_t e_face_dofs = assembler.get_e_face_dof();
            size_t a_face_dofs = assembler.get_a_face_dof();
            erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
            erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
            erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
            erk_an.refresh_faces_unknowns(x_dof);
                
            // ##################################################
            // ################################################## Preprocessor
            // ##################################################  
            
            std::ostringstream filename;
            filename << "explicit_l_" << l << "_n_" << sim_data.m_nt_divs << "_k_" << k << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
            std::string filename_str = filename.str();
            std::ofstream simulation_log(filename_str);
            sim_data.write_simulation_data(simulation_log);
            simulation_log << "Number of ERK steps =  " << s << std::endl;
            simulation_log << "Number of time steps =  " << nt << std::endl;
            simulation_log << "Step size =  " << dt << std::endl;
            simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
            simulation_log << "Space step = " << h << std::endl;
            simulation_log.flush();
            
            // ##################################################
            // ################################################## Time marching
            // ##################################################
            
            Matrix<RealType, Dynamic, 1> x_dof_n;
            for(size_t it = 1; it <= nt; it++) {
                RealType tn = dt*(it-1)+ti;
                // ERK step
                {
                    size_t n_dof = x_dof.rows();
                    Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
                    Matrix<double, Dynamic, 1> Fg, Fg_c, xd;
                    xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);    
                    Matrix<double, Dynamic, 1> yn, ki;
                    x_dof_n = x_dof;
                    for (int i = 0; i < s; i++) {
                        yn = x_dof;
                        for (int j = 0; j < s - 1; j++) 
                            yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                        t = tn + c(i,0) * dt;
                        {// Manufactured solution
                            auto v_fun      = functions.Evaluate_v(t);
                            auto f_fun      = functions.Evaluate_f(t);
                            auto s_v_fun    = functions.Evaluate_s_v(t);
                            auto s_f_fun    = functions.Evaluate_s_f(t);
                            assembler.get_e_bc_conditions().updateDirichletFunction(v_fun, 0);
                            assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
                            assembler.assemble_rhs(msh, f_fun, s_f_fun, true);
                            erk_an.SetFg(assembler.RHS);
                            erk_an.erk_weight(yn, ki);
                        }
                        // Accumulated solution
                        x_dof_n += dt*b(i,0)*ki;
                        k.block(0, i, n_dof, 1) = ki;
                    }
                }
                x_dof = x_dof_n;
                t = tn + dt;
                auto v_fun      = functions.Evaluate_v(t);
                auto flux_fun   = functions.Evaluate_sigma(t);
                auto s_v_fun    = functions.Evaluate_s_v(t);
                auto s_flux_fun = functions.Evaluate_s_q(t);
                if (it == nt) {
                    postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun, simulation_log);
                    postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic_energy_norm(msh, hho_di, assembler, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun, simulation_log);
                }
            }
            bool mesh_quality = true;
            if (mesh_quality) {
                std::ostringstream filename_mesh_quality;
                filename_mesh_quality << "mesh_quality_l_" << l << ".txt";
                std::string filename_str = filename_mesh_quality.str();
                std::ofstream mesh_file(filename_str);
                postprocessor<mesh_type>::mesh_quality(msh, assembler, mesh_file);
            }
        }
    }    
    std::cout << std::endl << std::endl;
}
        
