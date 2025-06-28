

//  Created by Romain Mottier

void ElasticIHHOFirstOrder(int argc, char **argv);

void ElasticIHHOFirstOrder(int argc, char **argv){
    
  // ###################################################################### Simulation paramaters 
  // ######################################################################
  
  using RealType = double;
  simulation_data sim_data = preprocessor::process_args(argc, argv);
  sim_data.print_simulation_data();
  timecounter tc, tcit;

  // ###################################################################### Mesh generation 
  // ######################################################################
  
  typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
  typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
  typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
  mesh_type msh;
  Mesh_generation(sim_data, msh);

  // ###################################################################### Time controls 
  // ######################################################################
  
  size_t nt = 10;
  for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
    nt *= 2;
  }
  RealType ti = 0.0;
  RealType tf = 0.5;
  RealType dt = (tf-ti)/nt;
  RealType t  = ti;
  
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

  for (auto & cell : msh ) {
    auto cell_ind = msh.lookup(cell);
    mesh_type::point_type bar = barycenter(msh, cell);
    elastic_material_data<RealType> material = elastic_mat_fun(bar);
    e_material.insert(std::make_pair(cell_ind,material));
  }
  
  size_t bc_elastic_id  = 0;
  size_t bc_acoustic_id = 1;
  for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
    auto face = *face_it;
    mesh_type::point_type bar = barycenter(msh, face);
    auto fc_id = msh.lookup(face);
    disk::boundary_descriptor bi{bc_elastic_id, true};
    msh.backend_storage()->boundary_info.at(fc_id) = bi;
    elastic_bc_face_indexes.insert(fc_id);
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
  
  // ###################################################################### Projecting initial data 
  // ######################################################################
  
  Matrix<RealType, Dynamic, 1> x_dof;
  assembler.project_over_cells(msh, x_dof, sin_elastic, null_flux_fun, null_s_fun, null_fun);
  assembler.project_over_faces(msh, x_dof, sin_elastic, null_s_fun);
  
  ////////// Post process of the initial data 
  if (sim_data.m_render_silo_files_Q) {
    size_t it = 0;
    std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
  }
  
  std::ofstream simulation_log("elasto_acoustic_inhomogeneous_four_fields.txt");
  std::ofstream energy_file("energy_file.txt");
  std::ofstream dissipation_file("dissipation_file.txt");

    std::ofstream Elastic_sensor_1_log("Elastic_s1_four_fields_h.csv");
    std::ofstream Elastic_sensor_2_log("Elastic_s2_four_fields_h.csv");
    std::ofstream Elastic_sensor_3_log("Elastic_s3_four_fields_h.csv");
    std::ofstream Elastic_sensor_4_log("Elastic_s4_four_fields_h.csv");
    bool e_side_Q = true;
    bool a_side_Q = false;
    typename mesh_type::point_type Elastic_s1_pt(-7.5, -2.5);
    typename mesh_type::point_type Elastic_s2_pt(-5.0, -2.5);
    typename mesh_type::point_type Elastic_s3_pt(-2.5, -2.5);
    typename mesh_type::point_type Elastic_s4_pt( 0.0, -2.5);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s1_pt_cell   = std::make_pair(Elastic_s1_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s2_pt_cell   = std::make_pair(Elastic_s2_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s3_pt_cell   = std::make_pair(Elastic_s3_pt, -1);
    std::pair<typename mesh_type::point_type,size_t> Elastic_s4_pt_cell   = std::make_pair(Elastic_s4_pt, -1);
    std::cout << bold << cyan << "      " << "Elastic sensor at (-7.5,-2.5)" << reset << std::endl; 
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
    std::cout << bold << cyan << "      " << "Elastic sensor at (-5.0,-2.5)" << reset << std::endl; 
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s2_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_2_log);
    std::cout << bold << cyan << "      " << "Elastic sensor at (-2.5,-2.5)" << reset << std::endl; 
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s3_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_3_log);  
    std::cout << bold << cyan << "      " << "Elastic sensor at (0.0,-2.5)" << reset << std::endl; 
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(0, Elastic_s4_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_4_log);  


  if (sim_data.m_report_energy_Q) {
    postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, t, x_dof, energy_file);
    postprocessor<mesh_type>::compute_elasto_acoustic_dissipation_four_field(assembler, t, x_dof, dissipation_file);
  }  
  
  // Solving a first order equation HDG/HHO propagation problem
  Matrix<RealType, Dynamic, Dynamic> a;
  Matrix<RealType, Dynamic, 1> b;
  Matrix<RealType, Dynamic, 1> c;
  
  // DIRK(s) schemes
  int s = 3;
  bool is_sdirk_Q = true;
  
  if (is_sdirk_Q) {
    dirk_butcher_tableau::sdirk_tables(s, a, b, c);
  } 
  else {
    dirk_butcher_tableau::dirk_tables(s, a, b, c);
  }
  
  std::cout << std::endl << std::endl;
  std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
  std::cout << bold << cyan << "      First stiffness assembly completed: ";
  tc.tic();
  assembler.assemble(msh, null_fun, null_s_fun);
  tc.toc();
  std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
  assembler.LHS += assembler.COUPLING; 
  dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
  
  if (sim_data.m_sc_Q) {
    std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(2);
    vec_cell_basis_data[0] = std::make_pair(assembler.get_e_material_data().size(), assembler.get_e_cell_basis_data());
    vec_cell_basis_data[1] = std::make_pair(assembler.get_a_material_data().size(), assembler.get_a_cell_basis_data());
    dirk_an.set_static_condensation_data(vec_cell_basis_data, assembler.get_n_face_dof());
  }
  
  if (is_sdirk_Q) {
    double scale = a(0,0) * dt;
    dirk_an.SetScale(scale);
    std::cout << bold << cyan << "      Matrix decomposed: "; 
    tc.tic();
    dirk_an.ComposeMatrix();
    //        dirk_an.setIterativeSolver();
    dirk_an.DecomposeMatrix();
    tc.toc();
    std::cout << tc << " seconds" << reset << std::endl;
  }
  
  // ##################################################
  // ################################################## Time marching
  // ##################################################
  
  std::cout << std::endl << std::endl;
  
  Matrix<RealType, Dynamic, 1> x_dof_n;
  
  for(size_t it = 1; it <= nt; it++) {
    

    tcit.tic();
    std::cout << bold << red << "   Time step number " << it << ": t = " << t 
              << reset;
              
    RealType tn = dt*(it-1)+ti;
    
    // DIRK step
    tc.tic();
    {
      size_t n_dof = x_dof.rows();
      Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
      Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
      xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
      
      RealType t;
      Matrix<double, Dynamic, 1> yn, ki;
      
      x_dof_n = x_dof;
      for (int i = 0; i < s; i++) {	
        yn = x_dof;
        for (int j = 0; j < s - 1; j++) {
          yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
        }
        t = tn + c(i,0) * dt;
        assembler.RHS.setZero();
        dirk_an.SetFg(assembler.RHS);
        dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);
        // Accumulated solution
        x_dof_n += dt*b(i,0)*ki;
        k.block(0, i, n_dof, 1) = ki;
      }
    }
    tc.toc();
    // std::cout << std::endl << bold << cyan << "      DIRK step completed: " << tc << " seconds" 
    //                                                            << reset << std::endl;

    x_dof = x_dof_n;

    // ##################################################
    // ################################################## Last postprocess
    // ##################################################

    if (sim_data.m_render_silo_files_Q) {
      std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
      postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
    }

      postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_1_log);
      postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s2_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_2_log);
      postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s3_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_3_log);
      postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, Elastic_s4_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, Elastic_sensor_4_log);

    t += dt;
    
    if (sim_data.m_report_energy_Q) {
      postprocessor<mesh_type>::compute_elasto_acoustic_energy_four_field(msh, hho_di, assembler, 
                                                                          t, x_dof, energy_file);
      postprocessor<mesh_type>::compute_elasto_acoustic_dissipation_four_field(assembler, t, x_dof, dissipation_file);
    } 

    std::cout << std::endl;
    tcit.toc();
    std::cout << bold << cyan << "      Iteration completed in " << tcit << " seconds" << reset << std::endl << std::endl;


  }
  sim_data.write_simulation_data(simulation_log);
  simulation_log << "Number of equations: " << dirk_an.DirkAnalysis().n_equations() << std::endl;
  simulation_log << "Number of DIRK steps =  " << s << std::endl;
  simulation_log << "Number of time steps =  " << nt << std::endl;
  simulation_log << "Step size =  " << dt << std::endl;
  simulation_log.flush();

}


























