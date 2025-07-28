

//  Created by Romain Mottier

// source /opt/intel/oneapi/setvars.sh intel64
// ../../../elastoacoustic -k 3 -s 0 -r 0 -c 0 -p 0 -l 0 -n 1400 -f 0 -e 0 

void ConicWavesEHHOFirstOrder(int argc, char **argv);

void ConicWavesEHHOFirstOrder(int argc, char **argv){
    
    // ######################################################################
    // ###################################################################### Simulation paramaters 
    // ######################################################################
  
    std::cout << std::endl << bold << red << "   EXPLICIT PULSE - COUPLING - GEOPHYSIC" << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, cpu;
    cpu.tic();
    DBSetDeprecateWarnings(0);
    
    // ##################################################
    // ################################################## Mesh generation 
    // ##################################################
  
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
        
        mesh_files.push_back("../../meshes/pulse/simplices/simplex_l2_0.4.txt");    // l = 0
        mesh_files.push_back("../../meshes/pulse/simplices/simplex_l3_0.21.txt");   // l = 1 
        mesh_files.push_back("../../meshes/pulse/simplices/simplex_l4_0.096.txt");  // l = 2
        mesh_files.push_back("../../meshes/pulse/simplices/simplex_l5_0.0485.txt"); // l = 3
        mesh_files.push_back("../../meshes/pulse/simplices/simplex_l6_0.024.txt");  // l = 4
        
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l2.txt");   // -l 0
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l3.txt");   // -l 1 
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l4.txt");   // -l 2
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l5.txt");   // -l 3
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l6.txt");   // -l 4
        // mesh_files.push_back("../../meshes/pulse/poly/poly_l7.txt");   // -l 5
        
        // Reading the polygonal mesh
        mesh_builder.set_poly_mesh_file(mesh_files[l]);
        mesh_builder.build_mesh();
        mesh_builder.move_to_mesh_storage(msh);
        mesh_builder.remove_duplicate_points();
    }
    else {
        // RealType lx = 5000.0;  
        // RealType ly = 3500.0;  
        RealType lx = 7500.0;  
        RealType ly = 5250.0;   
        // size_t nx = 140;
        // size_t ny = 140;
        size_t nx = 322/2; // k= 4 -> 224
        size_t ny = 322/2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx, ly, nx, ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        // mesh_builder.set_translation_data(-2500.0, -2500.0);
        mesh_builder.set_translation_data(-3750.0, -3750.0);
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
    std::cout << tc << " seconds" << reset << std::endl << std::endl;

    // ######################################################################
    // ###################################################################### Time controls 
    // ######################################################################
    
    size_t nt = sim_data.m_nt_divs;
    RealType ti = 0.0;
    RealType tf = 0.625;
    RealType dt = (tf-ti)/nt;
    RealType t  = ti;
  
    // ######################################################################
    // ###################################################################### HHO setting 
    // ######################################################################
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q)
        cell_k_degree++;
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);
    
    // ##################################################
    // ################################################## Material data 
    // ##################################################

    auto water_mat_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> acoustic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp;
        rho = 1025.0;            
        vp  = 1500.0;     
        acoustic_material_data<double> material(rho,vp);
        return material;
    };
    
    auto granit_mat_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> elastic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp, vs;
        rho = 2690.0; 
        vp  = 6000.0;    
        vs  = 3000.0; 
        elastic_material_data<double> material(rho,vp,vs);
        return material;
    };

    // ###################################################################### 
    // ###################################################################### Structure setting 
    // ###################################################################### 

    std::map<size_t,elastic_material_data<RealType>>  e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    
    RealType eps = 1.0e-10;
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (std::fabs(bar.y()) < eps) {
            interface_face_indexes.insert(fc_id);
            continue;
        }      
    }
    
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        
        // Assigning the material properties
        if (bar.y() > 0) {
            acoustic_material_data<RealType> material = water_mat_fun(bar);
            a_material.insert(std::make_pair(cell_ind,material));
        }
        else {
            elastic_material_data<RealType> material = granit_mat_fun(bar);
            e_material.insert(std::make_pair(cell_ind,material));
        }
        
        // Detection of faces on the interfaces
        auto cell_faces = faces(msh,cell);
        for (auto face :cell_faces) {
            auto fc_id = msh.lookup(face);
            bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
            if (is_member_Q) {
                if (bar.y() > 0) {
                    interface_cell_pair_indexes[fc_id].second = cell_ind;
                }
                else {
                    interface_cell_pair_indexes[fc_id].first = cell_ind;
                }
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
        if (is_member_Q) {
        }
        else {
            if (bar.y() > 0) {
                acoustic_internal_faces.insert(fc_id);
            }
            else {
                elastic_internal_faces.insert(fc_id);
            }
        }
    }

    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
        auto face = *face_it;
        mesh_type::point_type bar = barycenter(msh, face);
        auto fc_id = msh.lookup(face);
        if (bar.y() > 0) {
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
    
    auto null_s_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> double {
      return 0.0;
    }; 

    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
      disk::static_vector<double, 2> f{0,0};
      return f;
    };

    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_matrix<double,2,2> {
      double x,y;
      x = pt.x();
      y = pt.y();
      disk::static_matrix<double, 2, 2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
      return sigma;
    };

    // Boundary Condition
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
    
    // ##################################################
    // ################################################## Solving a primal HHO mixed problem 
    // ##################################################
  
    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if (sim_data.m_scaled_stabilization_Q)
        assembler.set_scaled_stabilization();
    
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
    std::cout << tc << " seconds" << reset << std::endl;    
  
    // ######################################################################
    // ###################################################################### Projecting initial data 
    // ######################################################################

    // Acoustic pulse - Heterogeneous domain - Fluide = Eau
    auto pulse_geophysic_v = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {    
        double x, y, xc, yc, r, wave, vp, lp, fc, c, vx, vy;
        x    = pt.x();
        y    = pt.y();
        xc   = 0.0;
        yc   = 100.0;
        fc   = 10.0;
        c    = 10.0;
        vp   = 1500.0;
        lp   = vp/fc;
        r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
        vx   = wave*(x-xc);
        vy   = wave*(y-yc);
        disk::static_vector<double, 2> v{vx,vy};
        return v;
    };
  
    Matrix<RealType, Dynamic, 1> x_dof;
    // Acoustic pulse intialized in pressure 
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, pulse_geophysic_v);
    assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
    
    // ##################################################
    // ################################################## Solving a first order equation HDG/HHO propagation problem
    // ##################################################

    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
  
    // ERK(s) schemes
    int s = 4;
    erk_butcher_tableau::erk_tables(s, a, b, c);
    std::cout << std::endl << bold << red << "   ASSEMBLY 2 : " << std::endl;
    std::cout << bold << cyan << "      First stiffness assembly completed: ";
    tc.tic();
    assembler.assemble(msh, null_fun, null_s_fun, true);
    tc.toc();
    std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
    assembler.LHS += assembler.COUPLING;
    
    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
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
    filename << "Explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".txt";
    std::string filename_str = filename.str();
    std::ofstream simulation_log(filename_str);
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Space step = " << h << std::endl;
    simulation_log.flush();

    size_t it = 0;
    std::ostringstream filename_silo;
    filename_silo << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
    std::string silo_file_name = filename_silo.str();
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);

    std::ostringstream filename_e;
    filename_e << "energy_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
    std::string filename_e_str = filename_e.str();
    std::ofstream energy_file(filename_e_str);
    if (sim_data.m_report_energy_Q) 
        postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);
    
    // ##################################################
    // ################################################## Sensors
    // ##################################################
  
    bool e_side_Q = true;
    bool a_side_Q = false;

    std::ostringstream filename_acou;
    filename_acou << "A_implicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_acou_str = filename_acou.str();
    std::ofstream Acoustic_sensor_1_log(filename_acou_str);
    typename mesh_type::point_type Acoustic_s1_pt(-2300,  100);
    std::pair<typename mesh_type::point_type,size_t> Acoustic_s1_pt_cell  = std::make_pair(Acoustic_s1_pt, -1);

    std::ostringstream filename_int;
    filename_int <<  "I_implicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_int_str = filename_int.str();
    std::ofstream Interface_sensor_1_log(filename_int_str);    
    typename mesh_type::point_type Interface_s1_pt(-1300, 0.0);
    std::pair<typename mesh_type::point_type,size_t> Interface_s1_pt_cell = std::make_pair(Interface_s1_pt, -1);

    std::ostringstream filename_ela;
    filename_ela <<  "E_implicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_ela_str = filename_ela.str();
    std::ofstream Elastic_sensor_1_log(filename_ela_str);
    typename mesh_type::point_type Elastic_s1_pt(-2300,  -200);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s1_pt_cell = std::make_pair(Elastic_s1_pt, -1);

    bool sensors = true;
    if (sensors) {
        postprocessor<mesh_type>::record_acoustic_data_elasto_acoustic_four_fields(0, Acoustic_s1_pt_cell, msh, hho_di, assembler, x_dof, a_side_Q, Acoustic_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Interface_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Interface_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
    }

    std::cout << std::endl;

    // ##################################################
    // ################################################## Time marching
    // ##################################################
    
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++) {
    
        tcit.tic();
        std::cout << bold << red << "   Time step number " << it << ": t = " << t << reset << std::endl;
        RealType tn = dt*(it-1)+ti;
        
        // ERK step
        tc.tic();
        {
            size_t n_dof = x_dof.rows();
            Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
            Matrix<RealType, Dynamic, 1> Fg, Fg_c, xd;
            xd = Matrix<RealType, Dynamic, 1>::Zero(n_dof, 1);
            
            Matrix<RealType, Dynamic, 1> yn, ki;
            x_dof_n = x_dof;
            for (int i = 0; i < s; i++) {
                yn = x_dof;
                for (int j = 0; j < s - 1; j++) {
                    yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
                }
                erk_an.erk_weight(yn, ki);
                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;      
            }
        }
        tc.toc();
        std::cout << bold << cyan << "      ERK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;

        // ##################################################
        // ################################################## Last postprocess
        // ##################################################

        int silo_mod = static_cast<int>(std::round(nt / 50.0)); // Number of silo files 
        if ((!sim_data.m_render_silo_files_Q && (it == 1 || it == std::round(nt/2) || it == nt)) || ((sim_data.m_render_silo_files_Q && (it%silo_mod == 0)) || it == nt)) {
            std::ostringstream filename;
            filename << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
            std::string silo_file_name = filename.str();
            postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
        }
        
        if (sensors) {
            postprocessor<mesh_type>::record_acoustic_data_elasto_acoustic_four_fields(it, Acoustic_s1_pt_cell, msh, hho_di, assembler, x_dof, a_side_Q, Acoustic_sensor_1_log);
            // postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Interface_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Interface_sensor_1_log);
            postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
        }

        if (sim_data.m_report_energy_Q) 
            postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);

        t += dt;
        
        tcit.toc();
        std::cout << bold << cyan << "      Iteration completed in " << tcit << " seconds" << reset << std::endl << std::endl;

    }
    
    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
    
}

