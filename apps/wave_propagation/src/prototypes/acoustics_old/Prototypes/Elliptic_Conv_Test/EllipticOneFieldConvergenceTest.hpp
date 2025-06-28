
//  Contributions by Omar Durán and Romain Mottier

#ifndef EllipticOneFieldConvergenceTest_hpp
#define EllipticOneFieldConvergenceTest_hpp

void EllipticOneFieldConvergenceTest(char **argv){

    simulation_data sim_data = preprocessor::process_convergence_test_lua_file(argv);
    sim_data.print_simulation_data();

    // Manufactured exact solution
    bool quadratic_function_Q = sim_data.m_exact_functions == 1;
    auto exact_scal_fun = [quadratic_function_Q](const mesh_type::point_type& pt) -> RealType {
        if(quadratic_function_Q){
            return (1.0-pt.x())*pt.x() * (1.0-pt.y())*pt.y();
        }else{
            return std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        }
        
    };

    auto exact_flux_fun = [quadratic_function_Q](const typename mesh_type::point_type& pt) -> std::vector<RealType> {
        double x,y;
        x = pt.x();
        y = pt.y();
        std::vector<RealType> flux(2);
        if(quadratic_function_Q){
            flux[0] = (1 - x)*(1 - y)*y - x*(1 - y)*y;
            flux[1] = (1 - x)*x*(1 - y) - (1 - x)*x*y;
            return flux;
        }else{
            flux[0] =  M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y());
            flux[1] =  M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y());
            return flux;
        }

    };

    auto rhs_fun = [quadratic_function_Q](const typename mesh_type::point_type& pt) -> RealType {
        double x,y;
        x = pt.x();
        y = pt.y();
        if(quadratic_function_Q){
            return -2.0*((x - 1)*x + (y - 1)*y);
        }else{
            return 2.0*M_PI*M_PI*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
        }
    };
    
    if (sim_data.m_polygonal_mesh_Q) {
        if (sim_data.m_n_divs > 8) {
            sim_data.m_n_divs = 8;
        }
    }

    // simple material
    RealType rho = 1.0;
    RealType vp = 1.0;
    acoustic_material_data<RealType> material(rho,vp);
    
    std::ofstream error_file("steady_scalar_polygon_error.txt");
    
    for(size_t k = 0; k <= sim_data.m_k_degree; k++){
        std::cout << bold << cyan << "Running an approximation with k : " << k << reset << std::endl;
        error_file << "Approximation with k : " << k << std::endl;
        for(size_t l = 0; l < sim_data.m_n_divs; l++){
            
            // Reading the polygonal mesh
            timecounter tc;
            
            tc.tic();
            mesh_type msh;
            if (sim_data.m_polygonal_mesh_Q) {
                polygon_2d_mesh_reader<RealType> mesh_builder;
                std::vector<std::string> mesh_files;
                mesh_files.push_back("meshes/unit_square_polymesh_nel_20.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_40.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_80.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_160.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_320.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_640.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_1280.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_2560.txt");
                mesh_files.push_back("meshes/unit_square_polymesh_nel_5120.txt");
                
                // Reading the polygonal mesh
                mesh_builder.set_poly_mesh_file(mesh_files[l]);
                mesh_builder.build_mesh();
                mesh_builder.move_to_mesh_storage(msh);
            }else{
                RealType lx = 1.0;
                RealType ly = 1.0;
                size_t nx = 2;
                size_t ny = 2;
                
                cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
                mesh_builder.refine_mesh(l);
                mesh_builder.build_mesh();
                mesh_builder.move_to_mesh_storage(msh);
            }
            std::cout << bold << cyan << "Mesh generation: " << tc << " seconds" << reset << std::endl;
            
            // Creating HHO approximation spaces and corresponding linear operator
            size_t cell_k_degree = k;
            if(sim_data.m_hdg_stabilization_Q){
                cell_k_degree++;
            }
            disk::hho_degree_info hho_di(cell_k_degree,k);

            // Solving a scalar primal HHO problem
            boundary_type bnd(msh);
            bnd.addDirichletEverywhere(exact_scal_fun);
            tc.tic();
            auto assembler = acoustic_one_field_assembler<mesh_type>(msh, hho_di, bnd);
            if(sim_data.m_hdg_stabilization_Q){
                assembler.set_hdg_stabilization();
            }
            assembler.load_material_data(msh,material);
            assembler.assemble(msh, rhs_fun);
            assembler.apply_bc(msh);
            tc.toc();
            std::cout << bold << cyan << "Assemble in : " << tc << " seconds" << reset << std::endl;
            
            // Solving LS
            Matrix<RealType, Dynamic, 1> x_dof;
            if (sim_data.m_sc_Q) {
                tc.tic();
                linear_solver<RealType> analysis(assembler.LHS,assembler.get_n_face_dof());
                analysis.condense_equations(std::make_pair(msh.cells_size(), assembler.get_cell_basis_data()));
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                if(sim_data.m_iterative_solver_Q){
                  analysis.set_iterative_solver(true);
                }
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                error_file << "Number of equations (SC) : " << analysis.n_equations() << std::endl;
            }else{
                tc.tic();
                linear_solver<RealType> analysis(assembler.LHS);
                tc.toc();
                std::cout << bold << cyan << "Create analysis in : " << tc << " seconds" << reset << std::endl;
                
                if(sim_data.m_iterative_solver_Q){
                  analysis.set_iterative_solver(true);
                }
                
                tc.tic();
                analysis.factorize();
                tc.toc();
                std::cout << bold << cyan << "Factorized in : " << tc << " seconds" << reset << std::endl;
                
                tc.tic();
                x_dof = analysis.solve(assembler.RHS);
                tc.toc();
                std::cout << bold << cyan << "Linear Solve in : " << tc << " seconds" << reset << std::endl;
                error_file << "Number of equations : " << analysis.n_equations() << std::endl;
            }
            
            // Computing errors
            postprocessor<mesh_type>::compute_errors_one_field(msh, hho_di, assembler, x_dof, exact_scal_fun, exact_flux_fun,error_file);
            
            if (sim_data.m_render_silo_files_Q) {
                std::string silo_file_name = "steady_scalar_k" + std::to_string(k) + "_";
                postprocessor<mesh_type>::write_silo_one_field(silo_file_name, l, msh, hho_di, x_dof, exact_scal_fun, false);
            }
        }
        error_file << std::endl << std::endl;
    }
    error_file.close();
}

#endif 
