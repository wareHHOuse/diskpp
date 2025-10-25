

//  Created by Romain Mottier

// ../wave_propagation -s0 -k3 -r0 -c0 -p0 -l4 -n300 -i0 -f0 -e0

void ERK4_LTS_conv_test(int argc, char **argv);

void ERK4_LTS_conv_test(int argc, char **argv){
  
    // #############################################################################################
    // ############################## Simulation paramaters ######################################## 
    // #############################################################################################
    
    std::cout << std::endl << bold << red << "   EXPLICIT ACOUSTIC CONV TEST" << reset << std::endl;
    
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, cpu;
    
    // #############################################################################################
    // ############################## Mesh generation ##############################################
    // #############################################################################################
    
    cpu.tic();

    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;
    
    if (sim_data.m_polygonal_mesh_Q) {
        size_t l = sim_data.m_n_divs;
        polygon_2d_mesh_reader<RealType> mesh_builder;
        std::vector<std::string> mesh_files;
        bool use_poly_mesh = false; 
        bool use_simp_mesh = false; 
        if (use_poly_mesh) {
            for (int i = 0; i <= 9; ++i) {
                mesh_files.push_back("../../meshes/conv_test/poly/poly_" + std::to_string(32 * (1 << i)) + ".txt");
            }
        } 
        else if (use_simp_mesh) {
            std::vector<double> conv_vals = {1.0, 0.35, 0.15, 0.07, 0.035, 0.026, 0.017, 0.0125, 0.0085, 0.005};
            for (int i = 0; i < conv_vals.size(); ++i) {
                mesh_files.push_back(
                    "../../meshes/conv_test/simplices/unstructured/l" + std::to_string(i) + "_conv_test_" + std::to_string(conv_vals[i]) + ".txt");
            }
        }
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

    
    // #############################################################################################
    // ################################ Time controls ##############################################
    // #############################################################################################
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) 
        nt  = sim_data.m_nt_divs;
    
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = (tf-ti)/nt;
    RealType t = ti;
    
    // #############################################################################################
    // ############################## Manufactured solution ########################################
    // #############################################################################################

    scal_vec_analytic_functions functions;
    functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionCubicInTimeAcoustic);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuarticInTimeAcoustic);
    // functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpaceAcoustic);
    
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

    // Acoustic analytical functions
    auto s_u_fun    = functions.Evaluate_s_u(t);
    auto s_v_fun    = functions.Evaluate_s_v(t);
    auto s_a_fun    = functions.Evaluate_s_a(t);
    auto s_f_fun    = functions.Evaluate_s_f(t);
    auto s_flux_fun = functions.Evaluate_s_q(t);
    
    // #############################################################################################
    // ################################## HHO setting ##############################################
    // #############################################################################################

    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if (sim_data.m_hdg_stabilization_Q) {
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree, sim_data.m_k_degree);
    
    // #############################################################################################
    // ################################ Material data ##############################################
    // #############################################################################################
    
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
    
    // #############################################################################################
    // ############################## Boundary conditions ##########################################
    // #############################################################################################

    std::map<size_t,elastic_material_data<RealType>> e_material;
    std::map<size_t,acoustic_material_data<RealType>> a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;

    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        mesh_type::point_type bar = barycenter(msh, cell);
        // Assigning the material properties
        acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
        a_material.insert(std::make_pair(cell_ind,material));
    }
    
    // Internal faces structure 
    std::set<size_t> elastic_internal_faces;
    std::set<size_t> acoustic_internal_faces;
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

    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);
    
    // #############################################################################################
    // ###################################### Assembly #############################################
    // #############################################################################################

    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_coupling_stabilization();
    if (sim_data.m_scaled_stabilization_Q) {
        assembler.set_scaled_stabilization();
    }    
    assembler.assemble_mass(msh);
    assembler.assemble_coupling_terms(msh);
    
    // #############################################################################################
    // ###################### Projecting initial data ##############################################
    // #############################################################################################

    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, s_v_fun, s_flux_fun);
    assembler.project_over_faces(msh, x_dof, null_fun, s_v_fun);
    
    // #############################################################################################
    // ###################################### Solving ##############################################
    // #############################################################################################

    Matrix<RealType, Dynamic, Dynamic> a;
    Matrix<RealType, Dynamic, 1> b;
    Matrix<RealType, Dynamic, 1> c;
    
    // ERK schemes
    int s = 4;
    erk_butcher_tableau::erk_tables(s, a, b, c);
    assembler.assemble(msh, null_fun, s_f_fun, true);
    assembler.LHS += assembler.COUPLING; 
    assembler.assemble_P(msh, 1);

    size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
    size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
    size_t e_face_dofs = assembler.get_e_face_dof();
    size_t a_face_dofs = assembler.get_a_face_dof();
    
    erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, assembler.COUPLING, elastic_cell_dofs, acoustic_cell_dofs, e_face_dofs, a_face_dofs);
    erk_an.Mcc_inverse(assembler.get_elastic_cells(), assembler.get_acoustic_cells(), assembler.get_e_cell_basis_data(), assembler.get_a_cell_basis_data());
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);//assembler.get_interfaces());
    erk_an.refresh_faces_unknowns(x_dof);
    
    // ##################################################
    // ################################################## Preprocessor
    // ##################################################  
    
    std::ostringstream filename;
    filename << "explicit_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_discret_" << sim_data.m_hdg_stabilization_Q << ".txt";
    std::string filename_str = filename.str();
    std::ofstream simulation_log(filename_str);
    sim_data.write_simulation_data(simulation_log);
    simulation_log << "Number of ERK steps =  " << s << std::endl;
    simulation_log << "Number of time steps =  " << nt << std::endl;
    simulation_log << "Step size =  " << dt << std::endl;
    simulation_log << "Number of equations : " << assembler.RHS.rows() << std::endl;
    simulation_log << "Space step = " << h << std::endl;
    simulation_log.flush();
    std::cout << std::endl << std::endl;

    size_t it = 0;
    std::ostringstream filename_silo;
    filename_silo << "silo_l_" << sim_data.m_n_divs << "_n_" << sim_data.m_nt_divs << "_k_" << sim_data.m_k_degree << "_s_" << s << "_";
    std::string silo_file_name = filename_silo.str();
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);

    // ##################################################
    // ################################################## Time marching
    // ##################################################
    
    size_t p = 1;
    size_t dtau = dt / p;
    auto l0 = x_dof;

    // DISCRTIZATION INFOS
    std::cout << bold << red << "   TIME LOOP: " << std::endl;
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++) {
        
        std::cout << bold << cyan << "      Time step number " << it << ": t = " << t << reset << std::endl;
        RealType tn = dt*(it-1)+ti;
        
        tc.tic();

        size_t n_dof = x_dof.rows();
        Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
        Matrix<double, Dynamic, Dynamic> k_T = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
        Matrix<double, Dynamic, Dynamic> k_F = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
        Matrix<double, Dynamic, 1> Fg, Fg_c, xd;
        xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
            
        Matrix<double, Dynamic, 1> yn, ki;
        x_dof_n = x_dof;
        
        //////////////////////////////////////////////////////////// SOURCE TERMS
        auto s_v_fun    = functions.Evaluate_s_v(t);
        auto s_f_fun_tn = functions.Evaluate_s_f(t);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
        assembler.assemble_rhs(msh, null_fun, s_f_fun, true);
        erk_an.SetFg(assembler.RHS);
        auto F0 = erk_an.invMc() * erk_an.Fc();

        auto t_dt = t+dt;
        s_v_fun = functions.Evaluate_s_v(t_dt);
        s_f_fun = functions.Evaluate_s_f(t_dt);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
        assembler.assemble_rhs(msh, null_fun, s_f_fun, true);
        erk_an.SetFg(assembler.RHS);
        auto F1 = erk_an.invMc() * erk_an.Fc();

        auto t_dt_2 = t+dt/2.0;
        s_v_fun = functions.Evaluate_s_v(t_dt_2);
        s_f_fun = functions.Evaluate_s_f(t_dt_2);
        assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
        assembler.assemble_rhs(msh, null_fun, s_f_fun, true);
        erk_an.SetFg(assembler.RHS);
        auto F2 = erk_an.invMc() * erk_an.Fc();

        ////////////////////////////////////////////////////////////// COARSE PREDICTOR
        Matrix<double, Dynamic, 1> w0, w1, w2, w3;
        Eigen::VectorXd tmp;
        tmp = assembler.IminusP * x_dof_n;
        // w0 
        w0 = erk_an.apply_B(tmp) + assembler.IminusP * F0;
        // w1
        tmp = assembler.IminusP * (erk_an.apply_B(x_dof_n) + F0);
        w1 = erk_an.apply_B(tmp) + assembler.IminusP * ((-3.0*F0 + 4.0*F2 - F1)/dt);
        // w2
        tmp = assembler.IminusP * (erk_an.apply_B_power(x_dof_n,2) + erk_an.apply_B(F0) + (-3.0*F0 + 4.0*F2 - F1)/dt);
        w2 = erk_an.apply_B(tmp) + assembler.IminusP * ((4.0*F0 - 8.0*F2 + 4.0*F1)/(dt*dt));

        // w3
        auto tmp1 = (-3*F0+4*F2-F1)/dt;
        tmp = assembler.IminusP * (erk_an.apply_B_power(x_dof_n,3) + erk_an.apply_B_power(F0,2) + erk_an.apply_B(tmp1) + ((4.0*F0 - 8.0*F2 + 4.0*F1)/(dt*dt)) );
        w3 = erk_an.apply_B(tmp);
        
        // LOOP ON LOCAL REFINEMENT RATIOS
        for (int m = 0; m < p; ++m) {            
            
            auto tm = tn + m * dtau;
            auto tm1 = tn + (m+1)*dtau;
            auto tm2 = tn + (m+0.5)*dtau;

            auto tmp1 = assembler.P*l0;     
            auto k1 = w0 + m*dtau*w1 + (m*m/2)*dtau*dtau*w2 + (m*m*m/6)*dtau*dtau*dtau*w3 + erk_an.apply_B(tmp1);

            auto tmp2 = assembler.P*(l0 + (dtau/2)*k1);     
            auto k2 = w0 + (m+0.5)*dtau*w1 + 0.5*(m+0.5)*(m+0.5)*dtau*dtau*w2 + (1/6)*(m+0.5)*(m+0.5)*(m+0.5)*dtau*dtau*dtau*w3 + erk_an.apply_B(tmp2);

            auto tmp3 = assembler.P*(l0 + (dtau/2)*k2);
            auto k3 = w0 + (m+0.5)*dtau*w1 + 0.5*(m+0.5)*(m+0.5)*dtau*dtau*w2 + (1/6)*(m+0.5)*(m+0.5)*(m+0.5)*dtau*dtau*dtau*w3 + erk_an.apply_B(tmp3);

            auto tmp4 = assembler.P*(l0 + dtau*k3);
            auto k4 = w0 + (m+1)*dtau*w1 + 0.5*(m+1)*(m+1)*dtau*dtau*w2 + (1/6)*(m+1)*(m+1)*(m+1)*dtau*dtau*dtau*w3 + erk_an.apply_B(tmp4);
            
            // ACCUMULATED SOLUTION
            auto l1 = l0 + (dtau/6) * (k1 + 2*k2 + 2*k3 + k4);
            l0 = l1;
            
        }
        
        tc.toc();
        std::cout << bold << yellow << "         LTS-ERK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof_n = l0;
        x_dof = x_dof_n;
        t = tn + dt;
        
        if(it == nt) {
            std::cout << std::endl;
            postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler, x_dof, null_fun, null_flux_fun, s_v_fun, s_flux_fun, simulation_log);
            postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic_energy_norm(msh, hho_di, assembler, x_dof, null_fun, null_flux_fun, s_v_fun, s_flux_fun, simulation_log);
        }

    }
    
    bool mesh_quality = false;
    if (mesh_quality) {
        std::ofstream mesh_file("mesh_file.txt");
        postprocessor<mesh_type>::mesh_quality(msh, assembler, mesh_file);
    }

    cpu.toc();
    simulation_log << "TOTAL CPU TIME: " << cpu << std::endl;
    std::cout << bold << red << "   TOTAL CPU TIME: " << cpu << std::endl << std::endl;
    
}

