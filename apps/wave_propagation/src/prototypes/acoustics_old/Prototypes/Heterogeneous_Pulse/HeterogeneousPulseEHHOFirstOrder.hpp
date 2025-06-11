
//  Contributions by Omar Durán and Romain Mottier

#ifndef HeterogeneousPulseEHHOFirstOrder_hpp
#define HeterogeneousPulseEHHOFirstOrder_hpp

void HeterogeneousPulseEHHOFirstOrder(char **argv){
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_acoustics_lua_file(argv);
    sim_data.print_simulation_data();
    
    // Building a cartesian mesh
    timecounter tc;
    tc.tic();

    RealType lx = 1.0;
    RealType ly = 1.0;
    size_t nx = 2;
    size_t ny = 2;
    
    
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
        
    }
    
    else {
    
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
    
    // Time controls : Final time value 0.25
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    RealType ti = 0.0;
    RealType tf = 0.25;
    RealType dt = tf/nt;
    
    auto null_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y;
            x = pt.x();
            y = pt.y();
            return 0.0;
    };
    
    auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        return {0,0};
    };
    
    auto vel_fun = [](const mesh_type::point_type& pt) -> RealType {
            RealType x,y,xc,yc,r,wave;
            x = pt.x();
            y = pt.y();
            xc = 0.5;
            yc = 0.25;
            r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
            wave = 0.1*(-4*std::sqrt(10.0/3.0)*(-1 + 1600.0*r*r))/(std::exp(800*r*r)*std::pow(M_PI,0.25));
            return wave;
    };
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

    // Solving a primal HHO mixed problem
    boundary_type bnd(msh);
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = acoustic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);

    auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> mat_data(2);
        RealType rho, vp;
        rho = 1.0;
        if (y < 0.5) {
            vp = 5.0;
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
    assembler.project_over_cells(msh, x_dof, vel_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, vel_fun);
    
    
    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::string silo_file_name = "e_inhomogeneous_scalar_mixed_";
        postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("inhomogeneous_acoustic_two_fields_explicit.txt");
        
    std::ofstream sensor_top_log("top_sensor_e_acoustic_two_fields.csv");
    std::ofstream sensor_bot_log("bot_sensor_e_acoustic_two_fields.csv");
    typename mesh_type::point_type top_pt(0.5, 2.0/3.0);
    typename mesh_type::point_type bot_pt(0.5, 1.0/3.0);
    std::pair<typename mesh_type::point_type,size_t> top_pt_cell = std::make_pair(top_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> bot_pt_cell = std::make_pair(bot_pt, -1);
    
    postprocessor<mesh_type>::record_data_acoustic_two_fields(0, top_pt_cell, msh, hho_di, x_dof, sensor_top_log);
    postprocessor<mesh_type>::record_data_acoustic_two_fields(0, bot_pt_cell, msh, hho_di, x_dof, sensor_bot_log);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }
    
    // Solving a first order equation HDG/HHO propagation problem
    int s = 4;
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    erk_butcher_tableau::erk_tables(s, a, b, c);

    tc.tic();
    assembler.assemble(msh, null_fun);
    tc.toc();
    std::cout << bold << cyan << "Stiffness and rhs assembly completed: " << tc << " seconds" << reset << std::endl;
    size_t n_face_dof = assembler.get_n_face_dof();
    tc.tic();
    erk_hho_scheme<RealType> erk_an(assembler.LHS,assembler.RHS,assembler.MASS,n_face_dof);
    erk_an.Kcc_inverse(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
    if(sim_data.m_hdg_stabilization_Q){
        erk_an.Sff_inverse(std::make_pair(assembler.get_n_faces(), assembler.get_face_basis_data()));
    }
    else {
      if (sim_data.m_iterative_solver_Q) {
        erk_an.setIterativeSolver();
      }
      erk_an.DecomposeFaceTerm();
    }
    tc.toc();
    std::cout << bold << cyan << "ERK analysis created: " << tc << " seconds" << reset << std::endl;
    
    erk_an.refresh_faces_unknowns(x_dof);
    Matrix<RealType, Dynamic, 1> x_dof_n;
    timecounter simulation_tc;
    simulation_tc.tic();
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        
        RealType tn = dt*(it-1)+ti;
        // ERK step
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
//                    RealType t = tn + c(i,0) * dt;
//                    auto exact_vel_fun      = functions.Evaluate_v(t);
//                    auto rhs_fun            = functions.Evaluate_f(t);
//                    assembler.get_bc_conditions().updateDirichletFunction(exact_vel_fun, 0);
//                    assembler.assemble_rhs(msh, rhs_fun);
//                    assembler.apply_bc(msh);
//                    erk_an.SetFg(assembler.RHS);
                    erk_an.erk_weight(yn, ki);
                }

                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        tc.toc();
        std::cout << bold << cyan << "ERK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        RealType t = tn + dt;
        
        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "e_inhomogeneous_scalar_mixed_";
            postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, vel_fun, null_flux_fun, false);
        }
        
        postprocessor<mesh_type>::record_data_acoustic_two_fields(it, top_pt_cell, msh, hho_di, x_dof, sensor_top_log);
        postprocessor<mesh_type>::record_data_acoustic_two_fields(it, bot_pt_cell, msh, hho_di, x_dof, sensor_bot_log);
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }
    }
    simulation_tc.toc();
    simulation_log << "Simulation time : " << simulation_tc << " seconds" << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Number of ERK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
}

#endif 
