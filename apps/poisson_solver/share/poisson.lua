mesh.source = "internal";
--mesh.source = "file";
mesh.type = "triangles";
--mesh.type = "quadrangles";
--mesh.type = "hexagons";

hho.variant = "mixed_order_low";
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

hho.order = 2;
local k00 = 1.0;
local k11 = 0.0000001;

function right_hand_side(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    return (k00+k11)*math.pi*math.pi*sx*sy;
end

local dp = tensor_2D:new()
dp:entry(0, 0, k00)
dp:entry(1, 1, k11)

function diffusion_coefficient(domain_num, x, y)
    return dp;
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

    table.insert(sol_infos, sol_info)

    return 0;
end

for i = 2,4 do
    mesh.refinement_level = 2;
    run()
end

print(" *** L2 error report *** ")

for i,si in ipairs(sol_infos) do
    io.write(i .. " " .. si.h .. " " .. si.L2err)
    if i > 1 then
        local num = math.log(prev_err/si.L2err);
        local den = math.log(prev_h/si.h)
        io.write(" " .. (num/den))
    end
    io.write("\n");
    prev_h = si.h;
    prev_err = si.L2err;
end

print(" *** Energy error report *** ")

for i,si in ipairs(sol_infos) do
    io.write(i .. " " .. si.h .. " " .. si.Aerr)
    if i > 1 then
        local num = math.log(prev_err/si.Aerr);
        local den = math.log(prev_h/si.h)
        io.write(" " .. (num/den))
    end
    io.write("\n");
    prev_h = si.h;
    prev_err = si.Aerr;
end