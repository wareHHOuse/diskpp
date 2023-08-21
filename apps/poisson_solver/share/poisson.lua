mesh.source = "internal";
mesh.type = "hexagons";

hho.variant = "mixed_order_high";
hho.use_stabfree = true;

boundary[0] = "dirichlet";
boundary[1] = "dirichlet";
boundary[2] = "dirichlet";
boundary[3] = "dirichlet";


function dirichlet_data(bnd_num, x, y)
    return 0;
end

function solution(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    return sx*sy;
end

function right_hand_side(domain_num, x, y, z)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local sz = 1.0;
    if sim.dimension == 3 then
        sz = math.sin(math.pi*z);
    end
    return sim.dimension*math.pi*math.pi*sx*sy*sz;
end

function diffusion_coefficient(domain_num, x, y)
    return 1;
end

sol_infos = {}

function solution_process()
    local sol_info = {};

    local tc = timecounter:new();

    --print(sim.dimension);
    
    tc:tic();
    assemble();
    sol_info.asmtime = tc:toc();
    print("Assembly time: " .. sol_info.asmtime .. " seconds");

    solve();

    export_to_visit();

    sol_info.h = mesh_h();
    sol_info.L2err, sol_info.Aerr = check();
    
    print("h: " .. sol_info.h .. ", L2 error: " .. 
        sol_info.L2err .. ", A error: " ..
        sol_info.Aerr);

    return 0;
end

--local dp = diffusion_parameter_2D:new(7)
--print(dp)
--dp:entry(1,1,4.2);
--print(dp:entry(1,1));
--print(dp)


mesh.refinement_level = 5;
hho.order = 0;

for i = 2,4 do
    mesh.refinement_level = i;
    run()
end

