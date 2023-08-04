mesh.source = "internal";
mesh.type = "triangles";
mesh.refinement_level = 3;

hho.order = 2;
hho.variant = "mixed_order_high";
hho.use_stabfree = true;

boundary[1] = "dirichlet";
boundary[2] = "dirichlet";
boundary[3] = "dirichlet";
boundary[4] = "dirichlet";

function dirichlet_data(bnd_num, x, y)
    return 0;
end

function neumann_data(bnd_num, x, y)
    return 0;
end

function right_hand_side(domain_num, x, y)
    return 0;
end

function diffusion_coefficient(domain_num, x, y)
    return 1;
end

function solution_process()
    assemble();
    return 0;
end

local dp = diffusion_parameter_2D:new(7)
print(dp)
dp:entry(1,1,4.2);
print(dp:entry(1,1));
print(dp)

for ord = 0,3 do
    hho.order = ord;
    run()
end
