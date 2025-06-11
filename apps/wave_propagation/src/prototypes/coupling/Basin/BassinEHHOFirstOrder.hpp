
//  Created by Romain Mottier

void BassinEHHOFirstOrder(int argc, char **argv);

void BassinEHHOFirstOrder(int argc, char **argv){
    
    // ######################################################################
    // ###################################################################### Simulation paramaters 
    // ######################################################################
  
    std::cout << std::endl << bold << red << "   EXPLICIT SEDIMENTARY BASIN" << std::endl;
    using RealType = double;
    simulation_data sim_data = preprocessor::process_args(argc, argv);
    sim_data.print_simulation_data();
    timecounter tc, tcit, cpu;
    cpu.tic();

    // ######################################################################
    // ###################################################################### Mesh generation 
    // ######################################################################
  
    tc.tic();
    typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
    typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
    typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
    mesh_type msh;
    
    auto validate_l = [](size_t l) -> size_t {
        if ((0 <= l) && (l < 1) ) {
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
    
    mesh_files.push_back("../../meshes/basin.txt");    // l = 0
    
    // Reading the polygonal mesh
    mesh_builder.set_poly_mesh_file(mesh_files[l]);
    mesh_builder.build_bassin();
    mesh_builder.move_to_mesh_storage(msh);
    mesh_builder.remove_duplicate_points();
    
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
    std::cout << tc.to_double() << " seconds" << reset << std::endl << std::endl;

    // ######################################################################
    // ###################################################################### Time controls 
    // ######################################################################
    
    size_t nt = 10;
    for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
        // nt *= 2;
        nt = sim_data.m_nt_divs;
    }
    
    RealType ti = 0.0;
    RealType tf = 1.0;
    RealType dt = (tf-ti)/nt;
    RealType t  = ti;
  
    // ######################################################################
    // ###################################################################### HHO setting 
    // ######################################################################
    
    // Creating HHO approximation spaces and corresponding linear operator
    size_t cell_k_degree = sim_data.m_k_degree;
    if(sim_data.m_hdg_stabilization_Q){
        cell_k_degree++;
    }
    disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);
    
    // ##################################################
    // ################################################## Material data 
    // ##################################################
  
    auto air_mat_fun = []() -> acoustic_material_data<double> {
        double rho, vp;
        rho = 1.225;            
        vp  = 343.0;     
        acoustic_material_data<double> material(rho,vp);
        return material;
    };

    auto basin_elastic_mat_fun = []() -> elastic_material_data<double> {
        double rho, vp, vs;
        rho = 2570.0; 
        vp  = 5350.0;    
        vs  = 3090.0; 
        elastic_material_data<double> material(rho,vp,vs);
        return material;
    };
    
    auto basin_sediment_mat_fun = []() -> elastic_material_data<double> {
        double rho, vp, vs;
        rho = 1200.0; 
        vp  = 3400.0;    
        vs  = 1400.0; 
        elastic_material_data<double> material(rho,vp,vs);
        return material;
    };

    // ###################################################################### 
    // ###################################################################### Structure setting -- A debug 
    // ###################################################################### 

    std::map<size_t, elastic_material_data<RealType>>  basin_e_material;
    std::map<size_t, acoustic_material_data<RealType>> basin_a_material;
    std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes;
    std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
    
    for (auto & cell : msh ) {
        auto cell_ind = msh.lookup(cell);
        if (mesh_builder.polygons[cell_ind].m_material == 1) {
            elastic_material_data<RealType> material = basin_elastic_mat_fun();
            basin_e_material.insert(std::make_pair(cell_ind,material));  
        }     
        if (mesh_builder.polygons[cell_ind].m_material == 2) {
            elastic_material_data<RealType> material = basin_sediment_mat_fun();
            basin_e_material.insert(std::make_pair(cell_ind,material));
        } 
        if (mesh_builder.polygons[cell_ind].m_material == 3) {
            acoustic_material_data<RealType> material = air_mat_fun();
            basin_a_material.insert(std::make_pair(cell_ind,material));
        } 
        
    }

    // Detect interfaces 
    disk::connectivity<disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>> conn(msh);
    for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
        const auto face = *face_it;
        auto fc_id = msh.lookup(face);
        auto connected_cells = conn.connected_cells(face);
        std::vector<int> cell_inds;
        std::vector<int> materials;
        for (const auto& cell : connected_cells) {
            auto cell_ind = msh.lookup(cell);
            cell_inds.push_back(cell_ind);
            materials.push_back(mesh_builder.polygons[cell_ind].m_material);
        }
        if (materials[0] == 3) {
            if (materials[1] == 1 || materials[1] == 2) {
                interface_cell_pair_indexes[fc_id].first = cell_inds[1];
                interface_cell_pair_indexes[fc_id].second = cell_inds[0];
                continue;
            }
        }
        if (materials[1] == 3) {
            if (materials[0] == 1 || materials[0] == 2) {
                interface_cell_pair_indexes[fc_id].first = cell_inds[0];
                interface_cell_pair_indexes[fc_id].second = cell_inds[1];
                continue;
            }
        }
    }
    
    auto null_s_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> double {
      return 0.0;
    }; 

    auto null_fun = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> static_vector<double, 2> {
      static_vector<double, 2> f{0,0};
      return f;
    };

    auto null_flux_fun = [](const typename disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> static_matrix<double,2,2> {
      double x,y;
      x = pt.x();
      y = pt.y();
      static_matrix<double, 2, 2> sigma = static_matrix<double,2,2>::Zero(2,2);
      return sigma;
    };
    
    // Boundary faces
    size_t bc_elastic_id  = 0;
    size_t bc_acoustic_id = 1;
    // Detect interface elastic - acoustic
    e_boundary_type e_bnd(msh);
    a_boundary_type a_bnd(msh);
    e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
    a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
  
    // ##################################################
    // ################################################## Solving a primal HHO mixed problem 
    // ##################################################
  
    tc.tic();
    auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd, basin_e_material, basin_a_material);
    assembler.set_interface_cell_indexes(interface_cell_pair_indexes);
    assembler.set_hdg_stabilization();
    if(sim_data.m_scaled_stabilization_Q){
        assembler.set_scaled_stabilization();
    }
    
    tc.toc();
    std::cout << bold << red << "   ASSEMBLY 1 : " << std::endl;
    std::cout << bold << cyan << "      Assembler generation : ";
    std::cout << tc.to_double() << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_mass(msh);
    tc.toc();
    std::cout << bold << cyan << "      Mass Assembly : ";
    std::cout << tc.to_double() << " seconds" << reset << std::endl;
    
    tc.tic();
    assembler.assemble_coupling_terms(msh);
    tc.toc();
    std::cout << bold << cyan << "      Coupling assembly : ";
    std::cout << tc.to_double() << " seconds" << reset << std::endl << std::endl;    
  
    // ######################################################################
    // ###################################################################### Projecting initial data 
    // ######################################################################
    
    auto pulse_basin = [](const disk::mesh<double, 2, disk::generic_mesh_storage<double, 2>>::point_type& pt) -> static_vector<double, 2> {
        double x,y,xc,yc,r,wave,vx,vy,c,lp, fc, vp;
        x    = pt.x();
        y    = pt.y();
        xc   = 0.0;
        yc   = 750.0;
        fc   = 25.0;
        c    = 10;
        vp   = 5350.0;
        lp   = vp/fc;
        r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
        wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
        vx   = wave*(x-xc);
        vy   = wave*(y-yc);
        static_vector<double, 2> v{vx,vy};
        return v;
    };
    
    Matrix<RealType, Dynamic, 1> x_dof;
    assembler.project_over_cells(msh, x_dof, pulse_basin, null_flux_fun, null_s_fun, null_fun);
    assembler.project_over_faces(msh, x_dof, pulse_basin, null_s_fun);
  
    // ##################################################
    // ################################################## Solving a first order equation HDG/HHO propagation problem
    // ##################################################
  
  
  // Solving a first order equation HDG/HHO propagation problem
  int s = 3;
  Matrix<RealType, Dynamic, Dynamic> a;
  Matrix<RealType, Dynamic, 1> b;
  Matrix<RealType, Dynamic, 1> c;
  erk_butcher_tableau::erk_tables(s, a, b, c);
  
  std::cout << std::endl << std::endl;
  std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
  std::cout << bold << cyan << "      First stiffness assembly completed: ";
  tc.tic();
  assembler.assemble(msh, null_fun, null_s_fun, true);
  tc.toc();
  std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
  assembler.LHS += assembler.COUPLING; 
  
  tc.tic();

  size_t elastic_cell_dofs  = assembler.get_e_n_cells_dof();
  size_t acoustic_cell_dofs = assembler.get_a_n_cells_dof();
  size_t face_dofs = assembler.get_n_face_dof();

  size_t e_cells = assembler.get_elastic_cells();
  size_t a_cells = assembler.get_acoustic_cells();
  size_t n_cells = e_cells + a_cells;

  ////////// Post process of the initial data 
  if (sim_data.m_render_silo_files_Q) {
    size_t it = 0;
    std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, basin_e_material, basin_a_material, false);
  }
  
  std::ofstream simulation_log("elasto_acoustic_inhomogeneous_four_fields.txt");
  std::ofstream energy_file("energy_file.txt");
  std::ofstream dissipation_file("dissipation_file.txt");

  if (sim_data.m_report_energy_Q) {
    postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field_bassin(mesh_builder, msh, hho_di, assembler, t, x_dof, energy_file);
  }  

  // erk_coupling_hho_scheme<RealType> erk_an(assembler.LHS, assembler.RHS, assembler.MASS, elastic_cell_dofs, acoustic_cell_dofs, face_dofs);
  // erk_an.Kcc_inverse(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
  // if (sim_data.m_hdg_stabilization_Q) {
  //   erk_an.Sff_inverse(std::make_pair(assembler.get_n_faces(), assembler.get_face_basis_data()));
  // }
  // else {
  //   if (sim_data.m_iterative_solver_Q) {
  //     erk_an.setIterativeSolver();
  //   }
  //   erk_an.DecomposeFaceTerm();
  // }
  // tc.toc();
  // std::cout << bold << cyan << "      ERK analysis created: " << tc << " seconds" << reset << std::endl;

  // erk_an.refresh_faces_unknowns(x_dof);
  
  // // ##################################################
  // // ################################################## Time marching
  // // ##################################################
  
  // std::cout << std::endl << std::endl;
  
  // Matrix<RealType, Dynamic, 1> x_dof_n;
  
  // for(size_t it = 1; it <= nt; it++) {
    

  //   tcit.tic();
  //   std::cout << bold << red << "   Time step number " << it << ": t = " << t 
  //             << reset;
              
  //   RealType tn = dt*(it-1)+ti;
    
  //   // ERK step
  //   tc.tic();
  //   {
  //     size_t n_dof = x_dof.rows();
  //     Matrix<RealType, Dynamic, Dynamic> k = Matrix<RealType, Dynamic, Dynamic>::Zero(n_dof, s);
  //     Matrix<RealType, Dynamic, 1> Fg, Fg_c,xd;
  //     xd = Matrix<RealType, Dynamic, 1>::Zero(n_dof, 1);
            
  //     Matrix<RealType, Dynamic, 1> yn, ki;

  //     x_dof_n = x_dof;
  //     for (int i = 0; i < s; i++) {                
  //       yn = x_dof;
  //       for (int j = 0; j < s - 1; j++) {
  //         yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
  //       }
  //       erk_an.erk_weight(yn, ki);
  //       // Accumulated solution
  //       x_dof_n += dt*b(i,0)*ki;
  //       k.block(0, i, n_dof, 1) = ki;
  //     }
  //   }
  //   tc.toc();
  //   // std::cout << std::endl << bold << cyan << "      DIRK step completed: " << tc << " seconds" 
  //   //                                                            << reset << std::endl;

  //   x_dof = x_dof_n;

  //   // ##################################################
  //   // ################################################## Last postprocess
  //   // ##################################################

  //   if (sim_data.m_render_silo_files_Q) {
  //     std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
  //     postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, bassin_e_material, bassin_a_material, false);
  //   }

  //   t += dt;
    
  //   if (sim_data.m_report_energy_Q) {
  //     postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field_bassin(mesh_builder, msh, hho_di, assembler, t, x_dof, energy_file);
  //     postprocessor<mesh_type>::compute_elasto_acoustic_dissipation_four_field(assembler, t, x_dof, dissipation_file);
  //   } 

  //   std::cout << std::endl;
  //   tcit.toc();
  //   std::cout << bold << cyan << "      Iteration completed in " << tcit << " seconds" << reset << std::endl << std::endl;


  // }
  // sim_data.write_simulation_data(simulation_log);
  // simulation_log << "Number of ERK steps =  " << s << std::endl;
  // simulation_log << "Number of time steps =  " << nt << std::endl;
  // simulation_log << "Step size =  " << dt << std::endl;
  // simulation_log.flush();

}





























