

//  Created by Romain Mottier

void HeterogeneousEHHOFirstOrderCFL(int argc, char **argv);

void HeterogeneousEHHOFirstOrderCFL(int argc, char **argv){
    
    // ###################################################################### Simulation paramaters 
    // ######################################################################
    
    std::cout << std::endl << bold << red << "   EXPLICIT PULSE - CFL" << std::endl << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, cpu;
    cpu.tic();
    
    // ###################################################################### Mesh generation 
    // ######################################################################
  
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

    // ###################################################################### Time controls 
    // ######################################################################
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        // nt *= 2;
        nt = sim_data.m_nt_divs;
    }
    RealType ti = 0.0;
    RealType tf = 0.25;
    RealType dt = (tf-ti)/nt;
    RealType t  = ti;
    
    // Loop over level of time step option -n  
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {

        RealType dt = (tf-ti)/nt;

        // ###################################################################### HHO setting 
        // ######################################################################
        
        // Creating HHO approximation spaces and corresponding linear operator
        size_t cell_k_degree = sim_data.m_k_degree;
        if(sim_data.m_hdg_stabilization_Q){
            cell_k_degree++;
        }
        disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);
        
        // ###################################################################### Assigning Material Data
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
                acoustic_material_data<RealType> material = acoustic_mat_fun(bar); 
                a_material.insert(std::make_pair(cell_ind,material));
            }
            else {
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
                        // interface_face_pair_indexes[fc_id].second = cp_fc;
                    }
                    else {
                        interface_cell_pair_indexes[fc_id].first = cell_ind;
                        // interface_face_pair_indexes[fc_id].first = cp_fc;
                    }
                }
                // cp_fc = cp_fc + 1;
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
                    // std::cout << fc_id << std::endl;
                }
                else {
                    elastic_internal_faces.insert(fc_id);
                    // std::cout << fc_id << std::endl;
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
        
        // Detect interface elastic - acoustic
        e_boundary_type e_bnd(msh);
        a_boundary_type a_bnd(msh);
        e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
        a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
        
        // ###################################################################### Solving a primal HHO mixed problem 
        // ######################################################################
        
        tc.tic();
        auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, e_material, a_material);
        assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
        
        // Stabilization type 
        if(sim_data.m_hdg_stabilization_Q){
            assembler.set_hdg_stabilization();
        }
        if(sim_data.m_scaled_stabilization_Q){
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
        std::cout << tc << " seconds" << reset << std::endl;    
        
        // ###################################################################### Projecting initial data 
        // ######################################################################
        
        Matrix<RealType, Dynamic, 1> x_dof;
        
        // Acoustic pulse intialized in velocity 
        assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_fun_adi_acoustic);
        assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
        
        // Elastic pulse intialized in velocity 
        // assembler.project_over_cells(msh, x_dof, v_fun_adi, null_flux_fun, null_s_fun, null_fun);
        // assembler.project_over_faces(msh, x_dof, v_fun_adi, null_s_fun);
        
        ////////// Post process of the initial data 
        if (sim_data.m_render_silo_files_Q) {
            size_t it = 0;
            std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
            postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
        }
        
        
        // Solving a first order equation HDG/HHO propagation problem
        Matrix<RealType, Dynamic, Dynamic> a;
        Matrix<RealType, Dynamic, 1> b;
        Matrix<RealType, Dynamic, 1> c;
        
  // ERK(s) schemes
  int s = 2;
  erk_butcher_tableau::erk_tables(s, a, b, c);
  
  std::cout << std::endl << std::endl;
  std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
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

  
  if(sim_data.m_hdg_stabilization_Q) {
    erk_an.Sff_inverse(assembler.get_elastic_faces(), assembler.get_acoustic_faces(), assembler.get_e_face_basis_data(), assembler.get_a_face_basis_data(), assembler.get_e_compress(), assembler.get_a_compress(), elastic_internal_faces, acoustic_internal_faces, interface_face_indexes);
  }
  else {
    erk_an.DecomposeFaceTerm();  
  }
  tc.toc();
  std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;
  tc.tic();
  erk_an.refresh_faces_unknowns(x_dof);
  tc.toc();
  std::cout << bold << cyan << "      Inverse of Sff + Coupling in: " << tc << " seconds" << reset << std::endl << std::endl;

  std::cout << bold << red << "   CFL (dt/h) =  " << dt/(lx/mesh_builder.get_nx()) << std::endl;

  // ##################################################
  // ################################################## Time marching
  // ##################################################
  
  std::cout << std::endl << std::endl;
  
  Matrix<RealType, Dynamic, 1> x_dof_n;

    dt = (tf-ti)/nt;        
    t = ti;
    
    bool approx_fail_check_Q = false;
    bool        fail_check_Q = false;

    std::ofstream energy_file("energy_file.txt");
    auto energy_0 = postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);  

    RealType energy = energy_0;

    for(size_t it = 1; it <= nt; it++) {
                
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
            
            // for (int i = 0; i < elastic_cell_dofs; ++i) {
            //     for (int j = 0; j < x_dof_n.cols(); ++j) {
            //         std::cout << std::setw(4) << x_dof_n(i, j) << "    ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl << std::endl << std::endl << std::endl;
        
        }
      }
      tc.toc();
      // std::cout << bold << cyan << "      ERK step completed: " << tc << " seconds" << reset << std::endl;
      x_dof = x_dof_n;

      // ##################################################
      // ################################################## Last postprocess
      // ##################################################


      RealType energy_n = postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);
      RealType relative_energy   = (energy_n - energy)   / energy;
      RealType relative_energy_0 = (energy_n - energy_0) / energy_0;

      if (sim_data.m_render_silo_files_Q) {
        std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
        postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
      }

      t += dt;

      bool unstable_check_Q = (relative_energy > 1e-2) || (relative_energy_0 >= 0.5e-2);
      // bool dissipation_check_Q = (relative_energy < -15e-2) || (relative_energy_0 <= -15e-2);

      // if (dissipation_check_Q) { // energy is increasing
      //   fail_check_Q = true;
      //   break;
      // }      
          
      if (unstable_check_Q) { // energy is increasing
        approx_fail_check_Q = true;
        break;
      }
          
      energy = energy_n;

    }
    
    // if (fail_check_Q) {
    //   std::cout << bold << red << "   Simulation dissipates too much energy" << reset << std::endl;
    //   break;
    // }

    if (approx_fail_check_Q) {
      std::cout << bold << red << "   Simulation is unstable" << reset << std::endl << std::endl;
      simulation_log << std::endl;
      simulation_log << "      Simulation is unstable for :"<< std::endl;
      simulation_log << "      Number of equations : " << assembler.RHS.rows() << std::endl;
      simulation_log << "      Number of ERK steps =  " << s << std::endl;
      simulation_log << "      Number of time steps =  " << nt << std::endl;
      simulation_log << "      dt size =  " << dt << std::endl;
      simulation_log << "      h size =  " << lx/mesh_builder.get_nx() << std::endl;
      simulation_log << "      CFL (dt/h) =  " << dt/(lx/mesh_builder.get_nx()) << std::endl;
      simulation_log << std::endl;
      simulation_log.flush();
      break;
    }
    else {
      simulation_log << "      Simulation is stable for :"<< std::endl;
      simulation_log << "      Number of equations : " << assembler.RHS.rows() << std::endl;
      simulation_log << "      Number of ERK steps =  " << s << std::endl;
      simulation_log << "      Number of time steps =  " << nt << std::endl;
      simulation_log << "      dt size =  " << dt << std::endl;
      simulation_log << "      h size =  " << lx/mesh_builder.get_nx() << std::endl;
      simulation_log << "      CFL (dt/h) =  " << dt/(lx/mesh_builder.get_nx()) << std::endl;
      simulation_log << std::endl;
      simulation_log.flush();
      std::cout << bold << red << "Simulation is stable for: " << nt << reset << std::endl << std::endl;
      nt -= 10;
      std::cout << bold << red << "Number of time steps: " << nt << reset << std::endl << std::endl;
      continue;
    }

  }



}





























