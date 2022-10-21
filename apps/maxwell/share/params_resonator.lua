
sim.mesh_filename = "msh.geo3s"
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

function reference_solution(tag, x, y, z)
    local re = math.sin(math.pi*x)*math.sin(math.pi*y);
    local im = 0;

    local Fx = complex.new(0, 0);
    local Fy = complex.new(0, 0);
    local Fz = complex.new(re, im);
    return complex3.new(Fx, Fy, Fz);
end

function reference_reconstruction(tag, x, y, z)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local cx = math.cos(math.pi*x);
    local cy = math.cos(math.pi*y);

    local Fx = complex.new( math.pi*sx*cy, 0);
    local Fy = complex.new(-math.pi*cx*sy, 0);
    local Fz = complex.new( 0, 0);

    return complex3.new(Fx, Fy, Fz);
end

function kappa2(tag, x, y, z)
    local omega = 2*math.pi*sim.frequency;
    local mu = materials.mu(tag, x, y, z)*const.mu0;
    local epsilon = materials.epsilon(tag, x, y, z)*const.eps0;
    return (omega*epsilon)*(omega*mu);
end

function source(tag, x, y, z)
    local k2 = kappa2(tag, x, y, z);
    local p = 2*math.pi*math.pi - k2;

    local re = p*math.sin(math.pi*x)*math.sin(math.pi*y);
    local im = 0;

    local Fx = complex.new(0, 0);
    local Fy = complex.new(0, 0);
    local Fz = complex.new(re, im);
    return complex3.new(Fx, Fy, Fz);

end

domain[1] = {};
domain[1].source = source;

function xxx()
    assemble();
    solve();
    save_to_silo();
    local err_e = compute_error();
    print("E error: " .. err_e);
    local err_h = compute_reconstruction_error();
    print("H error: " .. err_h);
end

