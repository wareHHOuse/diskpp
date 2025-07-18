local sfe = require("sfexp");

cfg = sfe.config:new();
cfg.vertex = 3;
cfg.degree = 0;
cfg.variant = 0
cfg.silo_fn = "explore.silo"

print(cfg)

local poly = make_point_vector();
poly:add( sfe.point:new(0.0, 0.0) )
poly:add( sfe.point:new(1.0, 0.0) )
poly:add( sfe.point:new(1.0, 1.0) )
poly:add( sfe.point:new(0.0, 1.0) )

run_poly(cfg, poly)


