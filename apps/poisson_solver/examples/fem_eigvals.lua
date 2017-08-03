config.method = "fem_eigs"
config.input_mesh = "/Users/matteo/Desktop/strangeshape/circle.mesh2d"
config.visit_output = "eigs_fem.silo"
config.eigval_output = "eigs_fem.txt"

-- Solver config for eigenvalue problems
solver.feast.verbose = true         -- Print solver logs
solver.feast.tolerance = 8          -- Tolerance: 8 means 10^-8
solver.feast.min_eigval = 1         -- Minimum value for eigenvalue search range
solver.feast.max_eigval = 200       -- Maximum value for eigenvalue search range
solver.feast.subspace_size = 64     -- Subspace size
