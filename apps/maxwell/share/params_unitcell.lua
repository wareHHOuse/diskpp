sim.mesh_filename = "unit_cell.geo3s"
sim.frequency = 200e6;
sim.order = 1;
silo.filename = "unit_cell.silo"

materials.epsilon = function(tag, x, y, z)
    if (tag == 2) then
        return 1.173, -0.1;
    end

    if (tag == 4) then
        return 11, 0
    end

    return 1, 0;
end

materials.mu = function(tag, x, y, z)
    --if (tag == 4) then
    --   return 1, -12.5
    --end

    return 1, 0;
end

materials.sigma = function(tag, x, y, z)
    return 0;
end

function plane_wave(tag, x, y, z)
    local Fx = complex.new(1, 0);
    local Fy = complex.new(0, 0);
    local Fz = complex.new(0, 0);
    return complex3.new(Fx, Fy, Fz);
end

boundary[5] = {};
boundary[5].kind = "impedance"
boundary[5].source = plane_wave;

boundary[2] = {};
boundary[2].kind = "neumann";

boundary[4] = {};
boundary[4].kind = "neumann";

boundary[26] = {};
boundary[26].kind = "neumann";

boundary[28] = {};
boundary[28].kind = "neumann";

boundary[42] = {};
boundary[42].kind = "neumann";

boundary[44] = {};
boundary[44].kind = "neumann";

boundary[64] = {};
boundary[64].kind = "neumann";

boundary[66] = {};
boundary[66].kind = "neumann";

function xxx()
    assemble();
    solve();
    save_to_silo();
end

