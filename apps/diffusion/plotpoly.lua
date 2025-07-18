local sfe = require("sfexp");

cfg = sfe.config:new();
cfg.vertex = 1;
cfg.degree = 0;
cfg.variant = 1 
cfg.silo_fn = "explore.silo"
cfg.use_stabfree = true
cfg.eps = 0.01
cfg.nsamples = 101 

print(cfg)

local poly = make_point_vector();
poly:add( sfe.point:new(0.1, 0.0) )
poly:add( sfe.point:new(0.9, 0.0) )
poly:add( sfe.point:new(1.0, 1.0) )
poly:add( sfe.point:new(0.5, 1.0) )
poly:add( sfe.point:new(0.0, 1.0) )

run_poly(cfg, poly)

--cfg.num_vertices = 5
--run(cfg)
