
-- The file performs a convergence test for an elliptic problem:
-- 	The approximation is constructed with implemented 
-- 	assemblers for the acoustic/elastic wave equation

config.max_k_deg = 4 -- <int>:  Maximum face polynomial degree: default 0
config.max_l_ref = 4 -- <int>:  Maximum number of uniform spatial refinements: default 0
config.stab_type = 0 -- <0-1>:  Stabilization type 0 (HHO), 1 (HDG-like): default 0
config.stab_scal = 1 -- <0-1>:  Stabilization scaling 0 O(1), 1 O(1/h_{f}): default 0
config.stat_cond = 1 -- <0-1>:  Static condensation: default 0
config.iter_solv = 0 -- <0-1>:  Iterative solver : default 0
config.exac_func = 0 -- <0-1>:  Manufactured function type 0 (non-polynomial), 1 (quadratic): default 0
config.poly_mesh = 0 -- <0-1>:  Use of polynoal meshes : default 0
config.silo_output = 0 -- <0-1>:  Write silo files : default 0

