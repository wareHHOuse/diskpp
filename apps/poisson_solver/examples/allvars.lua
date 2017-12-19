-- Method to use: can be one of the following
--  * fem:          solve a 2D poisson problem with homogeneous Dirichlet bc using FEM
--  * hho:          solve a 2D poisson problem with homogeneous Dirichlet bc using HHO
--  * fem_eigs      solve a 2D poisson eigenvalue problem with hom. Dirichlet bc using FEM
--  * hho_eigs      solve a 2D poisson eigenvalue problem with hom. Dirichlet bc using HHO
--  * reference     analytical solution of the eigenvalue problem. Requires unit square mesh
config.method = "fem_eigs"                      -- "fem" or "hho"

-- Method degree: for now valid only for HHO
config.degree = "1"

-- Mesh input file. Only triangular for now.
config.input_mesh = "/Users/matteo/Desktop/strangeshape/circle.mesh2d"

-- Output solution file (readable with VisIt or ParaView)
config.visit_output = "reference.silo"

-- The following two are not used yet
config.gnuplot_output = "reference.dat"
config.hdf5_output = "reference.h5"

-- Text file with numerical values of computed eigenvalues
config.eigval_output = "eigenvalues.txt"

-- Solver config for "fem" and "hho". Solver can be "cg" or "pardiso" for now.
solver.solver_type = "cg"
solver.cg.max_iter = 2000
solver.cg.rr_tol = 1e-9
solver.cg.rr_max = 100

-- Solver config for eigenvalue problems
solver.feast.verbose = true         -- Print solver logs
solver.feast.tolerance = 8          -- Tolerance: 8 means 10^-8
solver.feast.min_eigval = 1         -- Minimum value for eigenvalue search range
solver.feast.max_eigval = 200       -- Maximum value for eigenvalue search range
solver.feast.subspace_size = 64     -- Subspace size



-- Other stuff

tensor_ofs = 1
tensor_amplitude = 100
tensor_eps = 0.02


function material_parameter(x, y)
--    c = math.cos(math.pi*x/tensor_eps);
--    s = math.sin(math.pi*y/tensor_eps);
--    e = math.exp( (x*x + y*y)/2.0 )
--    return tensor_ofs + tensor_amplitude*c*c*s*s + e;
    return 3.0;
end    

function right_hand_side(x, y)
--    return math.sin(x)*math.sin(y);
    --return 2*math.pi*math.pi*math.sin(math.pi*x)*math.sin(math.pi*y);
    return math.sin(y)
end
