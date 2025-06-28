
//  Created by Romain Mottier

void EAcoustic_stability(int argc, char **argv);

void EAcoustic_stability(int argc, char **argv){
    
    // ##################################################
    // ################################################## Simulation paramaters 
    // ##################################################
    
    std::cout << std::endl << bold << red << "   ACOUSTIC STABILITY" << reset << std::endl;
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    
    // ##################################################
    // ################################################## Mesh generation 
    // ##################################################
    
    cpu.tic();
    tc.tic();
    
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;
    
    if (sim_data.m_polygonal_mesh_Q) {
        auto validate_l = [](size_t l) -> size_t {
            if ((0 <= l) && (l < 15) ) {
                return l;
            }
            else {
                std::cout << std::endl << std::endl;
                std::cout << "Warning:: Only few polygonal meshes available.";
                std::cout << std::endl << std::endl;
                return 4;
            }
        };
        
        size_t l = validate_l(sim_data.m_n_divs);
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l0_conv_test_1.0.txt");    // l = 0
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l1_conv_test_0.35.txt");   // 
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l2_conv_test_0.15.txt");   // l = 2
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l3_conv_test_0.07.txt");   // 
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l4_conv_test_0.035.txt");  // l = 4
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l5_conv_test_0.026.txt");  // 
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l6_conv_test_0.017.txt");  // l = 6
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l7_conv_test_0.0125.txt"); // 
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l8_conv_test_0.0085.txt"); // l = 8
        mesh_files.push_back("../../meshes/conv_test/simplices/unstructured/l9_conv_test_0.005.txt");  // 
        
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_32.txt");
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_64.txt");     // -l 1 
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_128.txt");    
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_256.txt");    // -l 3
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_512.txt");
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_1024.txt");   // -l 5 
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_2048.txt");
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_4096.txt");   // -l 7 
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_8192.txt");
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_16384.txt");  // -l 9
        // mesh_files.push_back("../../meshes/conv_test/poly/poly_32768.txt");
        
        // Reading the polygonal mesh
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
        if (h_l < h) {
            h = h_l;
        }
    }
    
    tc.toc();
    std::cout << bold << red << "   MESH GENERATION : ";
    std::cout << tc << " seconds" << reset << std::endl << std::endl << std::endl;
    
    // ##################################################
    // ################################################## Time controls 
    // ##################################################
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        // nt *= 2;
        nt = sim_data.m_nt_divs;
    }
    
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = (tf-ti)/nt;
    RealType t = ti;
     
    // ##################################################
    // ################################################## Manufactured solution 
    // ##################################################
    
    scal_vec_analytic_functions functions;
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInTime);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial_paper);
    
    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_matrix<double,2,2> {
        double x,y;
        x = pt.x();
        y = pt.y();
        disk::static_matrix<double, 2, 2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
        return sigma;
    };
    
    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
        disk::static_vector<double, 2> f{0,0};
        return f;
    };
    
    auto null_s_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> double {
        return 0.0;
    }; 
    
    // Acoustic analytical functions
    auto s_v_fun    = functions.Evaluate_s_v(t);
    auto s_flux_fun = functions.Evaluate_s_q(t);
    
    // ##################################################
    // ################################################## HHO setting 
    // ##################################################
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);
    
    // ##################################################
    // ################################################## Material data 
    // ##################################################
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> acoustic_material_data<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        RealType rho, vp;
        rho = 1.0; // Fluid mass density
        vp  = 1.0; // Vp
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
    
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
        a_material.insert(std::make_pair(cell_ind,material));
    }
    
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        disk::boundary_descriptor bi{bc_acoustic_id, true};
        msh.backend_storage()->boundary_info.at(fc_id) = bi;
        acoustic_bc_face_indexes.insert(fc_id);
    }
    
    // Internal faces structure (for explicit schemes)
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);      
        bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
    }
    
    // boundary conditions
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    
    // ##################################################
    // ################################################## Solving a primal HHO mixed problem 
    // ##################################################
    
    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    
    tc.toc();
    std::cout << bold << red << "   ASSEMBLY 1 : " << std::endl;
    std::cout << bold << cyan << "      Assembler generation : ";
    std::cout << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "      Mass Assembly : ";
    std::cout << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_coupling_terms(msh);
    tc.toc();
    std::cout << bold << cyan << "      Coupling assembly : ";
    std::cout << tc << " seconds" << reset << std::endl << std::endl;
    
    // ##################################################
    // ################################################## Projecting initial data 
    // ##################################################
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun);
    
    // ##################################################
    // ################################################## Solving a first order equation HDG/HHO propagation problem
    // ##################################################
    
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    // ERK(s) schemes
    int s = 4;
    erk_butcher_tableau::erk_tables(s, a, b, c);
    
    std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
    std::cout << bold << cyan << "      First stiffness assembly completed: ";
    tc.tic();
    assembler.assemble(msh, null_fun, null_s_fun, true);
    tc.toc();
    std::cout << tc << " seconds" << reset << std::endl;
    
    size_t elastic_cell_dofs  = 0.0;                            //assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = 0.0;                                   //assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();
    
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);//assembler.get_interfaces());
    
    tc.toc();
    std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;
    tc.tic();
    erk_an.refresh_faces_unknowns(x_dof); 
    tc.toc();
    std::cout << bold << cyan << "      Inverse of Sff + Coupling in: " << tc << " seconds" << reset << std::endl;
    
    // ##################################################
    // ################################################## Preprocessor
    // ##################################################
    
    std::ostringstream filename;
    filename << "A_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
    std::string filename_str = filename.str();
    std::ofstream simulation_log(filename_str);
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Space step = " << h << std::endl;
    // auto block_dimension = assembler.Scc_block_dimension();
    // erk_an.compute_eigenvalues_bis(assembler.LHS_STAB, block_dimension, simulation_log);
    erk_an.compute_eigenvalues(simulation_log);
    simulation_log.flush();
    
    std::ostringstream filename_e;
    filename_e << "energy_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
    std::string filename_e_str = filename_e.str();
    std::ofstream energy_file(filename_e_str);
    auto energy_0 = postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);
    
    // ##################################################
    // ################################################## Time marching
    // ##################################################
    
    Matrix<RealType, Dynamic, 1> x_dof_n;
    bool approx_fail_check_Q = false;
    RealType energy = energy_0;
    
    for(size_t it = 1; it <= nt; it++){    
        
        std::cout << bold << red << "   Time step number " << it << reset << std::endl;
        RealType tn = dt*(it-1)+ti;
        
        // ERK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c, xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
            
            Matrix<double, Dynamic, 1> yn, ki;
            x_dof_n = x_dof;
            
            for (int i = 0; i < s; i++) {
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }
                t = tn + c(i,0) * dt;
                {
                    assembler.assemble_rhs(msh, null_fun, null_s_fun, true);
                    erk_an.SetFg(assembler.RHS);
                    erk_an.erk_weight(yn, ki);
                }
                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        
        tc.toc();
        std::cout << bold << cyan << "      ERK step completed: " << tc << " seconds" << reset;
        x_dof = x_dof_n;
        
        t = tn + dt;
        
        // Energy evaluation
        RealType energy_n = postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);
        RealType relative_energy   = (energy_n - energy)   / energy;
        RealType relative_energy_0 = (energy_n - energy_0) / energy_0;
        bool unstable_check_Q = (relative_energy > 1.0e-2) || (relative_energy_0 >= 1.0e-2);
        if (unstable_check_Q) { 
            std::cout << std::endl << std::endl << bold << red << "   Simulation is unstable" << reset << std::endl;
            break;
        }
        energy = energy_n;
    }
    
    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << std::endl << bold << red << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
    
}

