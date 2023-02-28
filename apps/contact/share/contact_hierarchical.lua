config.degree_cell = 2
config.degree_face = 1
config.ref_degree_cell = 1
config.ref_degree_face = 1

config.input_mesh = "../../../diskpp/meshes/2D_triangles/netgen/tri01.mesh2d"
config.visit_output = "eigs_hho.silo"
config.eigval_output = "eigs_hho.txt"
config.compute_reference = true -- Run computations for reference solution and save it; otherwise it is read it

hs.levels = 7
hs.sol_level_min = 0
hs.sol_level_max = 6

hs.fem_theta = "p" -- Theta: p means 1; n means -1: z means 0
hs.fem_gamma = 100 -- Gamma: 5 means 10^5

hs.hho_theta = "n" -- Theta: p means 1; n means -1: z means 0
hs.hho_gamma = 100 -- Gamma: 5 means 10^5
hs.hho_trace = "l" -- Trace-type: f means face-based trace ; l means cell-based trace

-- Solver config for eigenvalue problems
solver.feast.verbose = false -- Print solver logs
solver.feast.tolerance = 10 -- Tolerance: 8 means 10^-8
--solver.feast.min_eigval = 1 -- Minimum value for eigenvalue search range
--solver.feast.max_eigval = 100 -- Maximum value for eigenvalue search range
--solver.feast.subspace_size = 32 -- Subspace sizeit was just in c
