config.degree_cell = 0
config.degree_face = 0

config.input_mesh = "../../../diskpp/meshes/2D_triangles/netgen/disk_int_1_1.mesh2d"


vpst.alpha = 0.0001; --ALG augmentation parameter
vpst.Lref = 1;    --Reference dimension
vpst.Vref = 1;    --Reference velocity
vpst.mu = 1;      --Viscosity
vpst.Bn = 0.9;    --Bingham number
vpst.f = 1;       --force
vpst.problem = "circular" --Choose only "circular","annulus", or "square" 

