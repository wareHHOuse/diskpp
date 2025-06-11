
-- The file performs a simulation for a heterogeneous acoustics problem with analytical solution:

config.fac_k_deg = 3 -- <int>:  Face polynomial degree: default 0
config.num_l_ref = 1 -- <int>:  Number of uniform spatial refinements: default 0
config.num_t_ref = 7 -- <int>:  Number of uniform time refinements: default 0
config.stab_type = 1 -- <0-1>:  Stabilization type 0 (HHO), 1 (HDG-like): default 0
config.stab_scal = 0 -- <0-1>:  Stabilization scaling 0 O(1), 1 O(1/h_{f}): default 0
config.stat_cond = 1 -- <0-1>:  Static condensation: default 0
config.iter_solv = 0 -- <0-1>:  Iterative solver : default 0
config.writ_ener = 1 -- <0-1>:  Report energy at each time value : default 0
config.poly_mesh = 1 -- <0-1>:  Use of polynoal meshes : default 0
config.silo_output = 1 -- <0-1>:  Write silo files : default 0

