config.degree_cell = 0
config.degree_face = 0

--config.input_mesh ="../../../diskpp/meshes/2D_quads/medit/square_h00125.medit2d"
config.input_mesh ="../../../diskpp/meshes/2D_quads/diskpp/testmesh-256-256.quad"

bi.hname = "256"; 
bi.alpha = 100;   --ALG augmentation parameter
bi.Lref  = 1;     --Reference dimension
bi.Vref  = 1;     --Reference velocity
bi.mu = 1;        --Viscosity
bi.Bn = 2.;       --Bingham number
bi.f  = 0;        --force
bi.problem = "DRIVEN" --Choose only "circular","annulus", or "square" 
