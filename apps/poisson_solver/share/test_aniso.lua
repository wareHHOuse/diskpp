mesh.source = "internal";
--mesh.source = "file";
--mesh.type = "triangles";
--mesh.type = "quadrangles";
mesh.type = "hexagons";

hho.variant = "equal_order";
hho.use_stabfree = false;
hho.dt_in_stab = false;
hho.order = 0;

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

local inv_aniso_ratio = 1;
local K = tensor_2D:new()

function right_hand_side(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local cx = math.cos(math.pi*x);
    local cy = math.cos(math.pi*y);
    local k00 = K:entry(0, 0);
    local k01 = K:entry(0, 1);
    local k10 = K:entry(1, 0);
    local k11 = K:entry(1, 1);
    local pi2 = math.pi*math.pi;
    return (k00+k11)*pi2*sx*sy - (k01+k10)*pi2*cx*cy;
end

function diffusion_coefficient(domain_num, x, y)
    return K;
end

--[[
function right_hand_side(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local cx = math.cos(math.pi*x);
    local cy = math.cos(math.pi*y);
    local s3x = math.sin(3*math.pi*x);
    local sx3 = sx*sx*sx;
    local cx2 = cx*cx;
    local pi = math.pi;
    local pi2 = pi*pi;

    local k00 = K:entry(0, 0);
    local k11 = K:entry(1, 1);
    
    local a = -3*k00*pi2*cx2*sx*sy;
    local b = pi2*(k11*sy/2 + (k11-k00)*cy/2)*(3*s3x - sx);
    local c = pi2*(k11-k00)*cx2*sx*cy;
    local d = -pi2*k00*sx3*sy;
    local e = -pi2*k11*sx*cx2*sy;
    return -(a+b+c+d+e);
end

function diffusion_coefficient(domain_num, x, y)
    local k00 = K:entry(0, 0);
    local k11 = K:entry(1, 1);

    local sx = math.sin(math.pi*x);
    local cx = math.cos(math.pi*x);
    local sx2 = sx*sx;
    local cx2 = cx*cx;

    local rK00 = k00*cx2 + k11*sx2;
    local rK01 = (k11-k00)*sx*cx;
    local rK10 = rK01;
    local rK11 = k00*sx2 + k11*cx2;

    local rK = tensor_2D:new();
    rK:entry(0, 0, rK00);
    rK:entry(0, 1, rK01);
    rK:entry(1, 0, rK10);
    rK:entry(1, 1, rK11);
    return rK;
end
]]--

sol_infos = {}

function solution_process()
    local sol_info = {};

    local tc = timecounter:new();
    tc:tic();
    assemble();
    sol_info.asmtime = tc:toc();
    print("Assembly time: " .. sol_info.asmtime .. " seconds");

    solve();

    export_to_visit();

    sol_info.k11 = inv_aniso_ratio;
    sol_info.h = mesh_h();
    sol_info.L2err, sol_info.Aerr = check();
    
    print("h: " .. sol_info.h .. ", L2 error: " .. 
        sol_info.L2err .. ", A error: " ..
        sol_info.Aerr);

    table.insert(sol_infos, sol_info)

    return 0;
end

function shorten_variant_name(long_vname)
    local vname = "";
    if long_vname == "mixed_order_low" then
        vname = "mol"
    end
    if long_vname == "equal_order" then
        vname = "eo"
    end
    if long_vname == "mixed_order_high" then
        vname = "moh"
    end
    return vname;
end

local rot = math.pi/8;
local r_c = math.cos(rot);
local r_s = math.sin(rot);

function test_mesh_refinement()
    for i = 1,5 do
        mesh.refinement_level = i
        run()
    end

    local vname = shorten_variant_name(hho.variant);
    
    local filename = "aniso/data/aniso_".. inv_aniso_ratio ..
        "_conv_" .. mesh.type .. "_" .. vname .. "_" ..
        hho.order .. ".txt";

    file = io.open(filename, "w");
    for i,si in ipairs(sol_infos) do
        file:write(si.h .. " " .. si.L2err .. " " .. si.Aerr .. "\n");
    end
    file:close();
end


local mesh_types = {"triangles", "quadrangles", "hexagons"}
local hho_vars = {"mixed_order_low", "equal_order", "mixed_order_high"}
local k11_vals = {1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001}

--local mesh_types = {"quadrangles"}
--local hho_vars = {"equal_order"}
--local k11_vals = {0.1}


for key0,mt in pairs(mesh_types) do
    mesh.type = mt;
    for key1,var in pairs(hho_vars) do
        hho.variant = var;

        local order_min = 0;
        local order_max = 3;
        if var == "mixed_order_low" then
            order_min = 1;
        end

        for key2,k11 in pairs(k11_vals) do
            inv_aniso_ratio = k11;
            local dx = 1;
            local dy = inv_aniso_ratio;
            K:entry(0,0, r_c * dx * r_c + r_s * dy * r_s);
            K:entry(0,1, r_s * dy * r_c - r_c * dx * r_s);
            K:entry(1,0, r_s * dy * r_c - r_c * dx * r_s);
            K:entry(1,1, r_s * dx * r_s + r_c * dy * r_c);
            for order = order_min,order_max do
                sol_infos = {}
                hho.order = order;
                test_mesh_refinement();
            end
        end
    end
end

