sfe = require("sfexp");

cfg = sfe.config:new();
cfg.num_vertices = 5;
cfg.vertex = 3;
cfg.degree = 0;
cfg.variant = 0
cfg.silo_fn = "exp.silo"

print(cfg)

run(cfg)


