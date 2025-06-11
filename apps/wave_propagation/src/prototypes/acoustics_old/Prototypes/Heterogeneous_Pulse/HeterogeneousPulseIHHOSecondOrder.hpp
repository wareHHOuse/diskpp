
//  Contributions by Omar Durán and Romain Mottier

#ifndef HeterogeneousPulseIHHOSecondOrder_hpp
#define HeterogeneousPulseIHHOSecondOrder_hpp

void HeterogeneousPulseIHHOSecondOrder(char **argv){
    
    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////////   Préprocessing   //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    
    simulation_data sim_data = preprocessor::process_acoustics_lua_file(argv);
    sim_data.print_simulation_data();
    
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////   Meshes   //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();
    
    mesh_type msh;
    if (sim_data.m_polygonal_mesh_Q) {
        
        auto validate_l = [](size_t l) -> size_t {
            if ((0 < l) && (l < 5) ) {
                return l;
            }else{
                std::cout << "Warning:: Only five meshes available. Running level 4." << std::endl;
                return 4;
            }
        };
        
        size_t l = validate_l(sim_data.m_n_divs);
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        mesh_files.push_back("meshes/mexican_hat_polymesh_nel_4096.txt");
        mesh_files.push_back("meshes/mexican_hat_polymesh_nel_5120.txt");
        mesh_files.push_back("meshes/mexican_hat_polymesh_nel_10240.txt");
        mesh_files.push_back("meshes/mexican_hat_polymesh_nel_16384.txt");
        mesh_files.push_back("meshes/mexican_hat_polymesh_nel_65536.txt");
        
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
    
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////   Time controls   ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = tf/nt;

    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////   Source term ?   //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    auto null_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            return 0.0;
    };
    
    auto vel_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y,xc,yc,r,wave;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 0.5;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = 0.1*(-4*std::sqrt(10.0/3.0)*(-1 + 1600.0*r*r))/(std::exp(800*r*r)*std::pow(M_PI,0.25));
            return wave;
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///   Creating HHO approximation spaces and corresponding linear operator   ///
    ///////////////////////////////////////////////////////////////////////////////
    
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////   Solving a primal HHO mixed problem   ///////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun); // easy because boundary assumes zero every where any time.
    tc.tic();
    auto assembler = acoustic_one_field_assembler<mesh_type>(msh, hho_di, bnd);

    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////   Acoustic properties   /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(2);
        RealType rho, vp;
        rho = 1.0;
        if (y < 0.5) {
            vp = 1.0;
        }else{
            vp = 1.0;
        }
        mat_data[0] = rho; // rho
        mat_data[1] = vp; // seismic compressional velocity vp
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
    
    ///////////////////////////////////////////////////////////////////////////////
    ///////////   Projecting initial scalar, velocity and acceleration   //////////
    ///////////////////////////////////////////////////////////////////////////////
    
    Matrix<RealType, Dynamic, 1> p_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, p_dof_n, null_fun);
    assembler.project_over_faces(msh, p_dof_n, null_fun);
    assembler.project_over_cells(msh, v_dof_n, vel_fun);
    assembler.project_over_faces(msh, v_dof_n, vel_fun);
    assembler.project_over_cells(msh, a_dof_n, null_fun);
    assembler.project_over_faces(msh, a_dof_n, null_fun);
 
    /////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////   Postprocessing  /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    	    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_scalar_";
        postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, v_dof_n, vel_fun, false);
    }
    
    std::ofstream simulation_log("inhomogeneous_acoustic_one_field.txt");
    
    std::ofstream sensor_log1("sensor1.csv");
    std::ofstream sensor_log2("sensor2.csv");
    std::ofstream sensor_log3("sensor3.csv");
    typename mesh_type::point_type pt1(1.0/3, 2.0/3);
    typename mesh_type::point_type pt2(1.0/2, 2.0/3);
    typename mesh_type::point_type pt3(2.0/2, 2.0/3);
    std::pair<typename mesh_type::point_type,size_t> pt1_cell = std::make_pair(pt1, -1);
    std::pair<typename mesh_type::point_type,size_t> pt2_cell = std::make_pair(pt2, -1);
    std::pair<typename mesh_type::point_type,size_t> pt3_cell = std::make_pair(pt3, -1);
    postprocessor<mesh_type>::record_data_acoustic_one_field(0, pt1_cell, msh, hho_di, v_dof_n, sensor_log1);
    postprocessor<mesh_type>::record_data_acoustic_one_field(0, pt2_cell, msh, hho_di, v_dof_n, sensor_log2);
    postprocessor<mesh_type>::record_data_acoustic_one_field(0, pt3_cell, msh, hho_di, v_dof_n, sensor_log3);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, ti, p_dof_n, v_dof_n, simulation_log);
    }
    
    timecounter simulation_tc;
    bool standar_Q = true;
    
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////   Newmark process   //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    {
        Matrix<RealType, Dynamic, 1> a_dof_np = a_dof_n;

        RealType beta = 0.25;
        RealType gamma = 0.5;
        if (!standar_Q) {
            RealType kappa = 0.25;
            gamma = 1.0;
            beta = kappa*(gamma+0.5)*(gamma+0.5);
        }
        
        tc.tic();
        assembler.assemble(msh, null_fun);
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
        
        simulation_tc.tic();
        
        
    	///////////////////////////////////////////////////////////////////////////////
    	////////////////////////////   Temporal loop   ////////////////////////////////
   	///////////////////////////////////////////////////////////////////////////////
        
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            // Manufactured solution
            RealType t = dt*it+ti;

            tc.tic();
            // Compute intermediate state for scalar and rate
            p_dof_n = p_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*p_dof_n;
            
            assembler.RHS.setZero(); // Optimization: this is a problem with (f = 0) and (p_D = 0)
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
            
            /////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////   Postprocessing  //////////////////////////////////
    	    /////////////////////////////////////////////////////////////////////////////////
    	    
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "inhomogeneous_scalar_";
                postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, v_dof_n, vel_fun, false);
            }
            
            postprocessor<mesh_type>::record_data_acoustic_one_field(it, pt1_cell, msh, hho_di, v_dof_n, sensor_log1);
            postprocessor<mesh_type>::record_data_acoustic_one_field(it, pt2_cell, msh, hho_di, v_dof_n, sensor_log2);
            postprocessor<mesh_type>::record_data_acoustic_one_field(it, pt3_cell, msh, hho_di, v_dof_n, sensor_log3);

            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, t, p_dof_n, v_dof_n, simulation_log);
            }
                	    
	    //postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(it, top_pt_cell, msh, hho_di, v_dof_n, sensor_top_log);
    	    //postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(it, bot_pt_cell, msh, hho_di, v_dof_n, sensor_bot_log);    	    
                       
            /////////////////////////////////////////////////////////////////////////////////
            //////////////////////////   Display Terminal  //////////////////////////////////
    	    /////////////////////////////////////////////////////////////////////////////////
        }
        simulation_tc.toc();
        simulation_log << "Simulation time : " << simulation_tc << " seconds" << std::endl;
        simulation_log << "Number of equations : " << analysis.n_equations() << std::endl;
        simulation_log << "Number of time steps =  " << nt << std::endl;
        simulation_log << "Step size =  " << dt << std::endl;
        simulation_log.flush();
    }
}

#endif 
