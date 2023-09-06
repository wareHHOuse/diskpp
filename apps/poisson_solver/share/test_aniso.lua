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

local dp = tensor_2D:new()

function right_hand_side(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local k00 = dp:entry(0, 0);
    local k11 = dp:entry(1, 1);
    return (k00+k11)*math.pi*math.pi*sx*sy;
end

function diffusion_coefficient(domain_num, x, y)
    return dp;
end

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

    sol_info.k11 = dp:entry(1,1);
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

function test_mesh_refinement()
    for i = 1,6 do
        mesh.refinement_level = i
        run()
    end

    local vname = shorten_variant_name(hho.variant);
    
    local k11 = dp:entry(1, 1);
    local filename = "aniso/data/aniso_".. k11 .. "_conv_" ..
        mesh.type .. "_" .. vname .. "_" ..
        hho.order .. ".txt";

    file = io.open(filename, "w");
    for i,si in ipairs(sol_infos) do
        file:write(si.h .. " " .. si.L2err .. " " .. si.Aerr .. "\n");
    end
    file:close();
end

function test_anisotropy()
    mesh.refinement_level = 4;
    local k11 = 1.0;
    dp:entry(1, 1, k11);
    
    for i = 1,12 do
        run()
        k11 = k11/10.0;
        dp:entry(1, 1, k11);
    end


    print(" *** L2 error report *** ")

    for i,si in ipairs(sol_infos) do
        io.write(i .. " " .. si.k11 .. " " .. si.L2err)
        if i > 1 then
            local num = math.log(prev_err/si.L2err);
            local den = math.log(prev_k11/si.k11)
            io.write(" " .. (num/den))
        end
        io.write("\n");
        prev_k11 = si.k11;
        prev_err = si.L2err;
    end

    print(" *** Energy error report *** ")

    for i,si in ipairs(sol_infos) do
        io.write(i .. " " .. si.k11 .. " " .. si.Aerr)
        if i > 1 then
            local num = math.log(prev_err/si.Aerr);
            local den = math.log(prev_k11/si.k11)
            io.write(" " .. (num/den))
        end
        io.write("\n");
        prev_k11 = si.k11;
        prev_err = si.Aerr;
    end

--[[
    local vname = short_variant_name();
    
    local filename = "aniso/data/aniso_kyy_" .. mesh.type .. "_" .. vname .. "_" .. hho.order .. ".txt";

    file = io.open(filename, "w");
    for i,si in ipairs(sol_infos) do
        file:write(si.k11 .. " " .. si.L2err .. " " .. si.Aerr .. "\n");
    end
    file:close();
--]]
end

--test_anisotropy();

local mesh_types = {"triangles", "quadrangles", "hexagons"}
local hho_vars = {"mixed_order_low", "equal_order", "mixed_order_high"}
local k11_vals = {1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001}

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
            dp:entry(1,1,k11);
            for order = order_min,order_max do
                sol_infos = {}
                hho.order = order;
                test_mesh_refinement();
            end
        end
    end
end

