config.degree_cell = 0
config.degree_face = 0

config.input_mesh = "../../../diskpp/meshes/2D_triangles/netgen/square_tri1.mesh2d"


vpst.alpha = 0.1; --ALG augmentation parameter
vpst.Lref = 1;    --Reference dimension
vpst.Vref = 1;    --Reference velocity
vpst.mu = 1;      --Viscosity
vpst.Bn = 0.3;    --Bingham number
vpst.f = 1;       --force
vpst.problem = "square" --Choose only "circular","annulus", or "square" 

