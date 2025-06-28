
//  Contributions by Omar Durán and Romain Mottier

#ifndef HeterogeneousGar6more2DIHHOSecondOrder_hpp
#define HeterogeneousGar6more2DIHHOSecondOrder_hpp

void HeterogeneousGar6more2DIHHOSecondOrder(int argc, char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 3.0;
    RealType ly = 3.0;
    size_t nx = 3;
    size_t ny = 3;
    
    mesh_type msh;

    cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
    mesh_builder.refine_mesh(sim_data.m_n_divs);
    mesh_builder.set_translation_data(-1.5, -1.5);
    mesh_builder.build_mesh();
    mesh_builder.move_to_mesh_storage(msh);
    
    std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    
    // Time controls : Final time value 1.0
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt     = tf/nt;

    auto null_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            return 0.0;
    };
    
    auto p_fun = [](const mesh_type::point_type& pt) -> RealType {
        RealType x,y,xc,yc,r,wave,vx,vy,v,c,lp,factor;
        x = pt.x();
        y = pt.y();
        xc = 0.0;
        yc = 2.0/3.0;
        c = 10.0;
        lp = 1.0*std::sqrt(9.0)/10.0;
        r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
        factor = (lp*lp/(2.0*M_PI*M_PI));
        return factor*wave;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun); // easy because boundary assumes zero every where any time.
    tc.tic();
    auto assembler = acoustic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
    
    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(2);
        RealType rho, vp;
        rho = 1.0;
        if (y < 0.0) {
            vp = 1.0*std::sqrt(3.0);
        }else{
            vp = 1.0*std::sqrt(9.0);
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
    
    // Projecting initial scalar, velocity and acceleration
    Matrix<RealType, Dynamic, 1> p_dof_n, v_dof_n, a_dof_n;
    assembler.project_over_cells(msh, p_dof_n, p_fun);
    assembler.project_over_faces(msh, p_dof_n, p_fun);
    assembler.project_over_cells(msh, v_dof_n, null_fun);
    assembler.project_over_faces(msh, v_dof_n, null_fun);
    assembler.project_over_cells(msh, a_dof_n, null_fun);
    assembler.project_over_faces(msh, a_dof_n, null_fun);
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "inhomogeneous_scalar_";
        postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, p_dof_n, p_fun, false);
    }
    
    std::ofstream simulation_log("inhomogeneous_acoustic_one_field.txt");
    
    std::ofstream sensor_1_log("s1_acoustic_one_field_h.csv");
    std::ofstream sensor_2_log("s2_acoustic_one_field_h.csv");
    std::ofstream sensor_3_log("s3_acoustic_one_field_h.csv");
    typename mesh_type::point_type s1_pt(-1.0/3.0, -1.0/3.0);
    typename mesh_type::point_type s2_pt( 0.0, -1.0/3.0);
    typename mesh_type::point_type s3_pt(+1.0/3.0, -1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
    
    postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(0, s1_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_1_log);
    postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(0, s2_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_2_log);
    postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(0, s3_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_3_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, ti, p_dof_n, v_dof_n, simulation_log);
    }
        
    timecounter simulation_tc;
    bool standar_Q = true;
    // Newmark process
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
        SparseMatrix<double> Kg = assembler.LHS;
        assembler.LHS *= beta*(dt*dt);
        assembler.LHS += assembler.MASS;
        tc.toc();
        std::cout << bold << cyan << "Stiffness assembly completed: " << tc << " seconds" << reset << std::endl;
        
        linear_solver<RealType> analysis;
        if (sim_data.m_sc_Q) {
            tc.tic();
            analysis.set_Kg(assembler.LHS,assembler.get_n_face_dof());
            analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
            tc.toc();
            std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
            
//            analysis.set_iterative_solver(true, 1.0e-10);
            analysis.set_direct_solver(true);
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
        }else{
            tc.tic();
            analysis.set_Kg(assembler.LHS);
            tc.toc();
            std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
            
            tc.tic();
            analysis.factorize();
            tc.toc();
            std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
        }
        
        simulation_tc.tic();
        for(size_t it = 1; it <= nt; it++){

            std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;

            // Manufactured solution
            RealType t = dt*it+ti;

            tc.tic();
            // Compute intermediate state for scalar and rate
            p_dof_n = p_dof_n + dt*v_dof_n + 0.5*dt*dt*(1-2.0*beta)*a_dof_n;
            v_dof_n = v_dof_n + dt*(1-gamma)*a_dof_n;
            Matrix<RealType, Dynamic, 1> res = Kg*p_dof_n;
            
            assembler.RHS.setZero();
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
                std::string silo_file_name = "inhomogeneous_scalar_";
                postprocessor<mesh_type>::write_silo_one_field(silo_file_name, it, msh, hho_di, p_dof_n, p_fun, false);
            }
            
            postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(it, s1_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_1_log);
            postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(it, s2_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_2_log);
            postprocessor<mesh_type>::record_velocity_data_acoustic_one_field(it, s3_pt_cell, msh, hho_di, assembler, p_dof_n, sensor_3_log);
            
            if (sim_data.m_report_energy_Q) {
                postprocessor<mesh_type>::compute_acoustic_energy_one_field(msh, hho_di, assembler, t, p_dof_n, v_dof_n, simulation_log);
            }
            
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
