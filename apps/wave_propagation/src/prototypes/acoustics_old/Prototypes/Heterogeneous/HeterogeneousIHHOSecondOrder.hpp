
//  Contributions by Omar Durán and Romain Mottier

#ifndef HeterogeneousIHHOSecondOrder_hpp
#define HeterogeneousIHHOSecondOrder_hpp

void HeterogeneousIHHOSecondOrder(char **argv){
    
    simulation_data sim_data = preprocessor::process_acoustics_lua_file(argv);
    sim_data.print_simulation_data();
    
    ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////   Building a cartesian mesh   //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 0.2;
    size_t nx = 10;
    size_t ny = 1;
    
    mesh_type msh;
    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh_x_direction(sim_data.m_n_divs);
    mesh_builder.set_translation_data(0.0, 0.0);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);

    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////   Time controls : Final time value 0.5   //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    
    size_t nt = 10; // Number of temp iterations
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.5;
    RealType dt = (tf-ti)/nt;

    
    scal_analytic_functions functions;
    functions.set_function_type(scal_analytic_functions::EFunctionType::EFunctionInhomogeneousInSpace);
    RealType t = ti;
    auto exact_scal_fun     = functions.Evaluate_u(t);
    auto exact_vel_fun      = functions.Evaluate_v(t);
    auto exact_accel_fun    = functions.Evaluate_a(t);
    auto exact_flux_fun     = functions.Evaluate_q(t);
    auto rhs_fun            = functions.Evaluate_f(t);
    
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////    Creating HHO approximation spaces and corresponding linear operator   //////
    ///////////////////////////////////////////////////////////////////////////////////////
    
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////   Solving a primal HHO mixed problem   ////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(exact_scal_fun); // easy because boundary assumes zero every where any time.
    tc.tic();
    auto assembler = acoustic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(2);
        RealType rho, vp;
        if (x < 0.5) {
            vp = 10.0;
        }else{
            vp = 1.0;
        }
        rho = 1.0/(vp*vp); // this is required to make both formulations compatible by keeping kappa = 1
        mat_data[0] = rho; // rho
        mat_data[1] = vp;  // seismic compressional velocity vp
        return mat_data;
    };
    
    assembler.load_material_data(msh,acoustic_mat_fun);
    if(sim_data.m_hdg_stabilization_Q){
        assembler.set_hdg_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    /////////////   Projecting initial scalar, velocity and acceleration   ////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    
    Matrix<RealType, Dynamic, 1> p_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, p_dof_n, exact_scal_fun);
    assembler.project_over_faces(msh, p_dof_n, exact_scal_fun);
    assembler.project_over_cells(msh, v_dof_n, exact_vel_fun);
    assembler.project_over_faces(msh, v_dof_n, exact_vel_fun);
    assembler.project_over_cells(msh, a_dof_n, exact_accel_fun);
    assembler.project_over_faces(msh, a_dof_n, exact_accel_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "scalar_";
        postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, v_dof_n, exact_vel_fun, false);
    }
    
    std::ofstream simulation_log("acoustic_one_field.txt");
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, t, p_dof_n, v_dof_n, simulation_log);
    }
    
    linear_solver<RealType> analysis;
    bool standar_Q = true;
        
    ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////   Newmark process  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.25;
        RealType gamma = 0.5;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 0.6;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, rhs_fun);
        SparseMatrix<RealType> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        linear_solver<RealType> analysis;
        if (sim_data.m_sc_Q) {
            analysis.set_Kg(assembler.LHS, assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
        }else{
            analysis.set_Kg(assembler.LHS);
        }
        
        if (sim_data.m_iterative_solver_Q) {
            analysis.set_iterative_solver(true);
        }else{
            analysis.set_direct_solver(true); // symmetric matrix case
        }
        
        analysis.factorize();
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        ///////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////   Temporal Loop  ///////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////  
             
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            // Manufactured solution
            RealType t = dt*it+ti;
            auto exact_scal_fun  = functions.Evaluate_u(t);
            auto exact_vel_fun   = functions.Evaluate_v(t);
            auto exact_flux_fun  = functions.Evaluate_q(t);

            assembler.get_bc_conditions().updateDirichletFunction(exact_scal_fun, 0);
            assembler.RHS.setZero(); // Optimization: this is a problem with (f = 0)
            assembler.apply_bc(msh);

            // Compute intermediate state for scalar and rate
            p_dof_n = p_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*p_dof_n;

            assembler.RHS -= res;
            tc.toc();
            std::cout << bold << cyan << "Rhs assembly completed: " << tc << " seconds" << reset << std::endl;

            tc.tic();
            a_dof_np = analysis.solve(assembler.RHS); // new acceleration
            tc.toc();
            
            std::cout << bold << cyan << "Solution completed: " << tc << " seconds" << reset << std::endl;

            // update scalar and rate
            p_dof_n += beta*dt*dt*a_dof_np;
            v_dof_n += gamma*dt*a_dof_np;
            a_dof_n  = a_dof_np;
            
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "scalar_";
                postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, v_dof_n, exact_vel_fun, false);
            }
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, t, p_dof_n, v_dof_n, simulation_log);
            }
            
            if(it == nt){
                postprocessor<mesh_type>::compute_errors_one_field(msh, hho_di, assembler, p_dof_n, exact_scal_fun, exact_flux_fun, simulation_log);
            }
            
        }
        
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
}

#endif 
