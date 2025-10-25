

//  Created by Romain Mottier
// ../wave_propagation -k 2 -s 0 -r 0 -c 0 -p 0 -l 6 -n 2500 -f 1 -e 0

void HeterogeneousERK4_LTS_HHO_FirstOrder(int argc, char **argv);

void HeterogeneousERK4_LTS_HHO_FirstOrder(int argc, char **argv){
    
    // ######################################################################
    // ###################################################################### Simulation paramaters 
    // ######################################################################
    
    std::cout << std::endl << bold << red << "   EXPLICIT PULSE - COUPLING" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, cpu;
    cpu.tic();
    
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
        RealType lx = 1;  
        RealType ly = 1;          
        size_t nx = 2;
        size_t ny = 2;
        cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
        mesh_builder.refine_mesh(sim_data.m_n_divs);
        mesh_builder.set_translation_data(-0.5, -0.5);
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
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        nt = sim_data.m_nt_divs;
    }
    
    RealType ti = 0.0;
    RealType tf = 0.25;
    RealType dt = (tf-ti)/nt;
    RealType t  = ti;
    
    // ######################################################################
    // ###################################################################### HHO setting 
    // ######################################################################
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);
    
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
    
    auto water_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> acoustic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp;
        rho = 1.0;            
        vp  = 1.0;     
        acoustic_material_data<double> material(rho,vp);
        return material;
    };
    
    auto granit_mat_fun_adi = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> elastic_material_data<double> {
        double x,y;
        x = pt.x();
        y = pt.y();
        double rho, vp, vs;
        rho = 2.624390244; 
        vp  = 4.0;    
        vs  = 2.0; 
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
            // acoustic_material_data<RealType> material = water_mat_fun_adi(bar); 
            acoustic_material_data<RealType> material = acoustic_mat_fun(bar); 
            a_material.insert(std::make_pair(cell_ind,material));
        }
        else {
            // elastic_material_data<RealType> material = granit_mat_fun_adi(bar); 
            elastic_material_data<RealType> material = elastic_mat_fun(bar); 
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

    // Boundary condition
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
    if (sim_data.m_scaled_stabilization_Q) {
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
  
    // ######################################################################
    // ###################################################################### Projecting initial data 
    // ######################################################################
    
    auto v_fun_adi_acoustic = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> disk::static_vector<double, 2> {
        double x, y, xc, yc, r, wave, vx, vy, c, lp, fc, vp;
        x    = pt.x();
        y    = pt.y();
        xc   = 0.0;
        yc   = 0.125; // 0.1;
        fc   = 10.0;
        c    = 10;
        vp   = std::sqrt(1.0);
        lp   = vp/fc;
        r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
        vx   = wave*(x-xc);
        vy   = wave*(y-yc);
        disk::static_vector<double, 2> v{vx,vy}; 
        return v;
    };

    auto p_fun_ricker = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> double {
        double x, y, xc, yc, r, fc, c, vp, lp, p0, pi2;
        x  = pt.x();
        y  = pt.y();
        xc = 0.0;
        yc = 0.15;
        fc = 10.0;
        c  = 10.0;
        vp = std::sqrt(1.0);
        lp  = vp / fc;           // longueur caractéristique
        p0  = 1.0;
        pi2 = M_PI * M_PI;
        
        r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        
        double arg = pi2 * r*r / (lp*lp);
        
        double p = p0 * (1.0 - 2.0 * arg) * std::exp(-arg);
        return p;
    };

    
    Matrix<RealType, Dynamic, 1> x_dof;
    // Acoustic pulse intialized in pressure 
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, p_fun_ricker, null_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, p_fun_ricker);
    // // Acoustic pulse intialized in velocity 
    // assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_fun_adi_acoustic);
    // assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
    // Elastic pulse intialized in velocity 
    // assembler.project_over_cells(msh, x_dof, v_fun, null_flux_fun, null_s_fun, null_fun);
    // assembler.project_over_faces(msh, x_dof, v_fun, null_s_fun);
  
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
    std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
    assembler.LHS += assembler.COUPLING; 
    assembler.assemble_P(msh, 1);

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
    simulation_log.flush();

    if (sim_data.m_render_silo_files_Q) {
        size_t it = 0;
        std::ostringstream filename;
        filename << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
        std::string silo_file_name = filename.str();
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
    }

    // ##################################################
    // ################################################## Sensors
    // ##################################################

    bool e_side_Q = true;
    bool a_side_Q = false;

    std::ostringstream filename_acou;
    filename_acou << "A_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_acou_str = filename_acou.str();
    std::ofstream Acoustic_sensor_1_log(filename_acou_str);
    typename mesh_type::point_type Acoustic_s1_pt(-0.15,  0.1);
    std::pair<typename mesh_type::point_type,size_t> Acoustic_s1_pt_cell  = std::make_pair(Acoustic_s1_pt, -1);

    std::ostringstream filename_int;
    filename_int <<  "I_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_int_str = filename_int.str();
    std::ofstream Interface_sensor_1_log(filename_int_str);    
    typename mesh_type::point_type Interface_s1_pt(-0.15, 0.0);
    std::pair<typename mesh_type::point_type,size_t> Interface_s1_pt_cell = std::make_pair(Interface_s1_pt, -1);

    std::ostringstream filename_ela;
    filename_ela <<  "E_explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << ".csv";
    std::string filename_ela_str = filename_ela.str();
    std::ofstream Elastic_sensor_1_log(filename_ela_str);
    typename mesh_type::point_type Elastic_s1_pt(-0.15,  -0.1);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s1_pt_cell = std::make_pair(Elastic_s1_pt, -1);

    bool sensors = false;
    if (sensors) {
        postprocessor<mesh_type>::record_acoustic_data_elasto_acoustic_four_fields(0, Acoustic_s1_pt_cell, msh, hho_di, assembler, x_dof, a_side_Q, Acoustic_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Interface_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Interface_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
    }

    std::cout << std::endl;

    // ##################################################
    // ################################################## Time marching
    // ##################################################
    
    size_t p = 1;
    size_t dtau = dt / p;
    auto l0 = x_dof;
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++) {
        
        tcit.tic();
        std::cout << bold << red << "   Time step number " << it << ": t = " << t << reset << std::endl;
        x_dof_n = x_dof;
        size_t n_dof = x_dof.rows();

        ////////////////////////////////////////////////////////////// COARSE PREDICTOR
        std::vector<Matrix<RealType, Dynamic, 1>>  w, k;           
        erk_an.compute_wi(x_dof_n, s, assembler.IminusP, w, k);

        for (int m = 0; m < p; ++m) {   

            auto w_stage_1 =   w[0] + (m*dtau)*w[1] + (0.5*m*m*dtau*dtau)*w[2] + (m*m*m*dtau*dtau*dtau)*w[3]/6.0;
            auto w_stage_2_3 = w[0] + ((m+0.5)*dtau)*w[1] + (0.5*(m+0.5)*(m+0.5)*dtau*dtau)*w[2] + ((m+0.5)*(m+0.5)*(m+0.5)*dtau*dtau*dtau)*w[3]/6.0;
            auto w_stage_4 =   w[0] + ((m+1.0)*dtau)*w[1] + (0.5*(m+1.0)*(m+1.0)*dtau*dtau)*w[2] + ((m+1.0)*(m+1.0)*(m+1.0)*dtau*dtau*dtau)*w[3]/6.0;

            auto y_stage_1 = assembler.P*x_dof;
            erk_an.erk_lts_weight(y_stage_1, k[0], w_stage_1);
            
            auto y_stage_2 = assembler.P*(x_dof + dtau/2.0 * k[0]);
            erk_an.erk_lts_weight(y_stage_2, k[1], w_stage_2_3);
            
            auto y_stage_3 = assembler.P*(x_dof + dtau/2.0 * k[1]);
            erk_an.erk_lts_weight(y_stage_3, k[2], w_stage_2_3);
            
            auto y_stage_4 = assembler.P*(x_dof + dtau * k[2]);
            erk_an.erk_lts_weight(y_stage_4, k[3], w_stage_4);
            
            x_dof_n = x_dof + dtau/6 * (k[0] + 2*k[1] + 2*k[2] + k[3]);

        }
        x_dof = x_dof_n;

        // ##################################################
        // ################################################## Last postprocess
        // ##################################################
        
        if (sim_data.m_render_silo_files_Q) {
            std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
            postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
        }
        
        postprocessor<mesh_type>::record_acoustic_data_elasto_acoustic_four_fields(it, Acoustic_s1_pt_cell, msh, hho_di, assembler, x_dof, a_side_Q, Acoustic_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Interface_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Interface_sensor_1_log);
        postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
        
        t += dt;
        
        tcit.toc();
        std::cout << bold << cyan << "      Iteration completed in " << tcit << " seconds" << reset << std::endl << std::endl;
        
        if (it == nt) {
            std::ostringstream silo;
            silo << "Silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
            std::string silo_str = silo.str();
            postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_str, it, msh, hho_di, x_dof, e_material, a_material, false);      
        }
        
        
    }
    
    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
    
}





























