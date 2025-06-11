
//  Contributions by Omar Durán and Romain Mottier

#ifndef IHHOFirstOrder_hpp
#define IHHOFirstOrder_hpp

void IHHOFirstOrder(char **argv){
    
    simulation_data sim_data = preprocessor::process_acoustics_lua_file(argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();
    
    if (sim_data.m_polygonal_mesh_Q) {
        if (sim_data.m_n_divs > 8) {
            sim_data.m_n_divs = 8;
        }
    }
    
    mesh_type msh;
    if (sim_data.m_polygonal_mesh_Q) {
        size_t l = sim_data.m_n_divs;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        mesh_files.push_back("meshes/unit_square_polymesh_nel_20.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_40.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_80.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_160.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_320.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_640.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_1280.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_2560.txt");
        mesh_files.push_back("meshes/unit_square_polymesh_nel_5120.txt");
        
        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
    }else{
        RealType lx = 1.0;
        RealType ly = 1.0;
        size_t nx = 2;
        size_t ny = 2;
        
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
    }
    
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = (tf-ti)/nt;

    scal_analytic_functions functions;
    switch (sim_data.m_exact_functions) {
        case 1:
            functions.set_function_type(scal_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);
            break;
        case 2:
            functions.set_function_type(scal_analytic_functions::EFunctionType::EFunctionQuadraticInTime);
            break;
        default:
            functions.set_function_type(scal_analytic_functions::EFunctionType::EFunctionNonPolynomial);
            break;
    }

    RealType t = ti;
    auto exact_vel_fun  = functions.Evaluate_v(t);
    auto exact_flux_fun = functions.Evaluate_q(t);
    auto rhs_fun        = functions.Evaluate_f(t);
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_vel_fun);
    tc.tic();
    auto assembler = acoustic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);
    
    // simple material
    RealType rho = 1.0;
    RealType vp = 1.0;
    acoustic_material_data<RealType> material(rho,vp);
    assembler.load_material_data(msh,material);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    // Projecting initial data
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, exact_vel_fun, exact_flux_fun);
    assembler.project_over_faces(msh, x_dof, exact_vel_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "scalar_mixed_";
        postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
    }
    
    std::ofstream simulation_log("acoustic_two_fields.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    // DIRK(s) schemes
    int s = 3;
    bool is_sdirk_Q = true;
    
    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }
    
    tc.tic();
    assembler.assemble(msh, rhs_fun);
    tc.toc();
    std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
    dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
    
    if (sim_data.m_sc_Q) {
        dirk_an.set_static_condensation_data(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()), assembler.get_n_face_dof());
    }
    
    if (is_sdirk_Q) {
        // SDIRK case
        double scale = a(0,0) * dt;
        dirk_an.SetScale(scale);
        tc.tic();
        dirk_an.ComposeMatrix();
        
        if (sim_data.m_iterative_solver_Q) {
            dirk_an.setIterativeSolver();
        }
        
        dirk_an.DecomposeMatrix();
        tc.toc();
        std::cout << bold << cyan << "Matrix decomposed: " << tc << " seconds" << reset << std::endl;
    }else{
        // DIRK case
        if (sim_data.m_iterative_solver_Q) {
            dirk_an.setIterativeSolver();
        }
    }
    
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;
        
        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<RealType, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<RealType, Dynamic, 1>::Zero(n_dof, 1);
            
            Matrix<RealType, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {
                
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }
                
                {
                    RealType t = tn + c(i,0) * dt;
                    auto exact_vel_fun      = functions.Evaluate_v(t);
                    auto rhs_fun            = functions.Evaluate_f(t);
                    assembler.get_bc_conditions().updateDirichletFunction(exact_vel_fun, 0);
                    assembler.assemble_rhs(msh, rhs_fun);
                    assembler.apply_bc(msh);
                    dirk_an.SetFg(assembler.RHS);
                    dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);
                }

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "DIRK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;
        
        t = tn + dt;
        auto exact_vel_fun = functions.Evaluate_v(t);
        auto exact_flux_fun = functions.Evaluate_q(t);
        
        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "scalar_mixed_";
            postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, exact_vel_fun, exact_flux_fun, false);
        }
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }
        
        if(it == nt){
            // Computing errors
            postprocessor<mesh_type>::compute_errors_two_fields(msh, hho_di, assembler, x_dof, exact_vel_fun, exact_flux_fun,simulation_log);
        }

    }
    
    simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
    simulation_log << "Number of DIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
    
}

#endif 
