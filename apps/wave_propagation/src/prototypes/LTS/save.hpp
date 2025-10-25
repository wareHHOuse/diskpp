

//  Created by Romain Mottier

// ../../elastoacoustic -k 1 -s 0 -r 0 -c 1 -p 1 -l 8 -n 1000 -f 0 -e 0
void ELTSAcousticFirstOrder(int argc, char **argv);

void ELTSAcousticFirstOrder(int argc, char **argv){
  
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
                mesh_files.push_back(
                    "../../meshes/conv_test/poly/poly_" + std::to_string(32 * (1 << i)) + ".txt"
                );
            }
        } 
        else if (use_simp_mesh) {
            std::vector<double> conv_vals = {1.0, 0.35, 0.15, 0.07, 0.035, 0.026, 0.017, 0.0125, 0.0085, 0.005};
            for (int i = 0; i < conv_vals.size(); ++i) {
                mesh_files.push_back(
                    "../../meshes/conv_test/simplices/unstructured/l" + std::to_string(i) +
                    "_conv_test_" + std::to_string(conv_vals[i]) + ".txt"
                );
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
        if (h_l < h) 
            h = h_l;
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
    
    // DISCRTIZATION INFOS
    std::cout << bold << red << "   TIME LOOP: " << std::endl;
    Matrix<RealType, Dynamic, 1> x_dof_n;
    for(size_t it = 1; it <= nt; it++) {
        
        std::cout << bold << cyan << "      Time step number " << it << ": t = " << t << reset << std::endl;
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
                    // Manufactured solution
                    auto s_v_fun    = functions.Evaluate_s_v(t);
                    auto s_f_fun    = functions.Evaluate_s_f(t);
                    
                    assembler.get_e_bc_conditions().updateDirichletFunction(null_fun, 0);
                    assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
                    assembler.assemble_rhs(msh, null_fun, s_f_fun, true);
                    erk_an.SetFg(assembler.RHS);
                    erk_an.erk_weight(yn, ki);
                }
                
                // Accumulated solution
                x_dof_n += dt*b(i,0)*ki;
                k.block(0, i, n_dof, 1) = ki;
            }
        }
        
        tc.toc();
        std::cout << bold << yellow << "         ERK step completed: " << tc << " seconds" << reset << std::endl;
        x_dof = x_dof_n;
        
        t = tn + dt;
        auto s_v_fun    = functions.Evaluate_s_v(t);
        auto s_flux_fun = functions.Evaluate_s_q(t);
        
        if(it == nt){
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

