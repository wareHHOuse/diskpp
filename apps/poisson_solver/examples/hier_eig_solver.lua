config.degree = 0
config.input_mesh = "../../../eigvals_paper_meshes/l_shaped/netgen/l_shape3.mesh2d"
config.visit_output = "eigs_hho.silo"
config.eigval_output = "eigs_hho.txt"

hs.levels = 4
hs.eigvec = 2

hs.sol_level_min = 0
hs.sol_level_max = 3

hs.stab_weight = 2*config.degree + 3
-- Solver config for eigenvalue problems
solver.feast.verbose = true         -- Print solver logs
solver.feast.tolerance = 8          -- Tolerance: 8 means 10^-8
solver.feast.min_eigval = 1         -- Minimum value for eigenvalue search range
solver.feast.max_eigval = 100       -- Maximum value for eigenvalue search range
solver.feast.subspace_size = 32     -- Subspace size
