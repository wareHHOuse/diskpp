
sim.mesh_filename = "msh2.geo3s"
sim.frequency = 300e6
sim.order = 1

silo.filename = "output.silo"

materials.epsilon = function (tag, x, y, z)
    return 1, 0;
end

materials.mu = function (tag, x, y, z)
    return 1, 0
end

materials.sigma = function (tag, x, y, z)
    return 0
end

function plane_wave(tag, x, y, z)
    local Fx = complex.new(0, 0);
    local Fy = complex.new(2, 0);
    local Fz = complex.new(0, 0);
    return complex3.new(Fx, Fy, Fz);
end

function reference_solution(tag, x, y, z)
    re = math.cos(2*math.pi*z);
    im = -math.sin(2*math.pi*z);

    local Fx = complex.new(0, 0);
    local Fy = complex.new(re, im);
    local Fz = complex.new(0, 0);
    return complex3.new(Fx, Fy, Fz);
end

local neumann_bcs = {1,2,7,8}
for i,v in ipairs(neumann_bcs) do
    boundary[v] = {}
    boundary[v].kind = "neumann";
end

local impedance_bcs = {5, 11}
for i,v in ipairs(impedance_bcs) do
    boundary[v] = {}
    boundary[v].kind = "impedance";
end

boundary[6] = {}
boundary[6].kind = "tfsf";
boundary[6].source = plane_wave;

--boundary[11].value = 377/2;

domain = {}
domain[1] = {}
domain[1].scattered_field = true;

--impedances = {377, 377/2, 377/4}
--frequencies_MHz = {75, 150, 300, 350, 400, 450, 500, 550, 600};

impedances = { 377/2 }
frequencies_MHz = { 75 }

function xxx()
    for i,v in ipairs(frequencies_MHz) do
        for j, z in ipairs(impedances) do
            sim.frequency = v*1e6;
            boundary[11].value = z;
            assemble();
            solve();
            silo.filename = "output_" .. v .. ".silo";
            save_to_silo();
            compute_return_loss(6);

            gamma = (377 - z)/(377 + z);
            local rl = 20*math.log10(gamma);
            print("Expected RL: " .. rl)
            --local err = compute_error();
            --print("Error: ", err);
        end
    end
end

