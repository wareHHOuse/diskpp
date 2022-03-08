config.degree_cell = 0
config.degree_face = 0

--config.input_mesh ="../../../diskpp/meshes/2D_quads/medit/square_h00125.medit2d"
config.input_mesh ="../../../diskpp/meshes/2D_triangles/medit/vane_sym_h0125.medit2d"


bi.alpha = 1;   --ALG augmentation parameter
bi.Lref  = 4.02;     --Reference dimension R 
bi.Vref  = 4.02;     --Reference velocity V = omega * R
bi.mu = 1;        -- Viscosity
bi.Bn = 1;        -- Bingham number
bi.f  = 0;        -- force
bi.problem = "VANE" --Choose only "circular","annulus", or "square" 
