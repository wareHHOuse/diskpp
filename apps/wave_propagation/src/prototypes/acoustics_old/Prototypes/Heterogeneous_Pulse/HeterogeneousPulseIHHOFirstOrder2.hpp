
//  Contributions by Omar Durán and Romain Mottier

#ifndef HeterogeneousPulseIHHOFirstOrder_hpp
#define HeterogeneousPulseIHHOFirstOrder_hpp

void HeterogeneousPulseIHHOFirstOrder(char **argv){
    
    ///////////////////////////////////////////////////////////////////////////////
    //////////////////////////////   Préprocessing   //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    using RealType = double;
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
        
        // A commenter pour la lecture des différents maillages polyhédriques réguliers 
        /* polygon_2d_mesh_reader<RealType> mesh_builder;
        std::string mesh_file = "meshes/mexican_hat_polymesh_adapted_nel_8512.txt";
        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_file);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);*/
        
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

    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////////   Time controls   ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    // Time controls : Final time value 0.50
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt *= 2;
    }
    size_t it = 0;
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = tf/nt;

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////   Initial Condition   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
        
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
    bnd.addDirichletEverywhere(null_fun);
    tc.tic();
    auto assembler = acoustic_two_fields_assembler<mesh_type>(msh, hho_di, bnd);

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
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    tc.toc();
    std::cout << bold << cyan << "Assembler generation: " << tc << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "Mass Assembly completed: " << tc << " seconds" << reset << std::endl;
        
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////   Projecting initial data   /////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun);
    assembler.project_over_faces(msh, x_dof, null_fun);
 
    //  // Source Term
    typename mesh_type::point_type pt_source(0.5,0.5);
    std::pair<typename mesh_type::point_type,size_t> source_cell = std::make_pair(pt_source, -1);
    // source_term<mesh_type>::ricker_fluid(it,dt,source_cell,msh,hho_di,x_dof);

    /////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////   Postprocessing  /////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    
    if (sim_data.m_render_silo_files_Q) {
        std::string silo_file_name = "inhomogeneous_scalar_mixed_";
        postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, null_fun, null_flux_fun, false);
    }
    
    std::ofstream simulation_log("inhomogeneous_acoustic_two_fields.txt");
    
    std::ofstream sensor_log1("sensor1.csv");
    std::ofstream sensor_log2("sensor2.csv");
    std::ofstream sensor_log3("sensor3.csv");
    typename mesh_type::point_type pt1(1.0/3, 2.0/3);
    typename mesh_type::point_type pt2(1.0/2, 2.0/3);
    typename mesh_type::point_type pt3(2.0/3, 2.0/3);
    std::pair<typename mesh_type::point_type,size_t> pt1_cell = std::make_pair(pt1, -1);
    std::pair<typename mesh_type::point_type,size_t> pt2_cell = std::make_pair(pt2, -1);
    std::pair<typename mesh_type::point_type,size_t> pt3_cell = std::make_pair(pt3, -1);
    postprocessor<mesh_type>::record_data_acoustic_two_fields(0, pt1_cell, msh, hho_di, x_dof, sensor_log1);
    postprocessor<mesh_type>::record_data_acoustic_two_fields(0, pt2_cell, msh, hho_di, x_dof, sensor_log2);
    postprocessor<mesh_type>::record_data_acoustic_two_fields(0, pt3_cell, msh, hho_di, x_dof, sensor_log3);
    
    if (sim_data.m_report_energy_Q) {
        postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, ti, x_dof, simulation_log);
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///////   Solving a first order equation HDG/HHO propagation problem   ////////
    ///////////////////////////////////////////////////////////////////////////////
    
    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////   DIRK(s) schemes   //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////    
    
    int s = 3;
    bool is_sdirk_Q = true;
    
    if (is_sdirk_Q) {
        dirk_butcher_tableau::sdirk_tables(s, a, b, c);
    }else{
        dirk_butcher_tableau::dirk_tables(s, a, b, c);
    }
    
    tc.tic();
    assembler.assemble(msh, null_fun);
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

    Matrix<double, Dynamic, 1> x_dof_n;
    timecounter simulation_tc;
    simulation_tc.tic();
        
    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////   Temporal loop   ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
   	
    for(size_t it = 1; it <= nt; it++){

        std::cout << bold << yellow << "Time step number : " << it << " being executed." << reset << std::endl;
        RealType tn = dt*(it-1)+ti;


        
        // DIRK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
            xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
            
            Matrix<double, Dynamic, 1> yn, ki;

            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {
                
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }
                
                {
                    assembler.RHS.setZero();
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
        source_term<mesh_type>::punctual_source(it,dt,source_cell,msh,hho_di,x_dof_n);
        x_dof = x_dof_n;

        
        RealType t = tn + dt;
        
        /////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////   Postprocessing  //////////////////////////////////
    	/////////////////////////////////////////////////////////////////////////////////
        
        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "inhomogeneous_scalar_mixed_";
            postprocessor<mesh_type>::write_silo_two_fields(silo_file_name, it, msh, hho_di, x_dof, null_fun, null_flux_fun, false);
        }
        
        postprocessor<mesh_type>::record_data_acoustic_two_fields(it, pt1_cell, msh, hho_di, x_dof, sensor_log1);
        postprocessor<mesh_type>::record_data_acoustic_two_fields(it, pt2_cell, msh, hho_di, x_dof, sensor_log2);
        postprocessor<mesh_type>::record_data_acoustic_two_fields(it, pt3_cell, msh, hho_di, x_dof, sensor_log3);
        
        if (sim_data.m_report_energy_Q) {
            postprocessor<mesh_type>::compute_acoustic_energy_two_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
        }

    }
    simulation_tc.toc();
    simulation_log << "Simulation time : " << simulation_tc << " seconds" << std::endl;
    simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
    simulation_log << "Number of DIRK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log.flush();
}

#endif 
