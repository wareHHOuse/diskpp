
sim.mesh_filename = "msh.geo3s"
sim.frequency = 300e6
sim.order = 2

silo.filename = "output.silo"

materials.epsilon = function (tag, x, y, z)
    return 1, 0
end

materials.mu = function (tag, x, y, z)
    return 1, 0
end

materials.sigma = function (tag, x, y, z)
    return 0
end

function plane_wave(tag, x, y, z)
    local Fx = complex.new(0, 0);
    local Fy = complex.new(1, 0);
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

boundary[1] = {}
boundary[1].kind = "neumann";

boundary[2] = {}
boundary[2].kind = "neumann";

boundary[5] = {}
boundary[5].kind = "impedance";
boundary[5].source = plane_wave;

boundary[6] = {}
boundary[6].kind = "impedance";

frequencies_MHz = {300};

function xxx()
    for i,v in ipairs(frequencies_MHz) do
        sim.frequency = v*1e6;
        assemble();
        solve();
        silo.filename = "output_" .. v .. ".silo";
        save_to_silo();
        local err = compute_error();
        print("Error: ", err);
    end
end

