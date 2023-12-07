
sim.mesh_filename = "full_validation.geo3s"
sim.frequency = 270e6
sim.order = 1

silo.filename = "output.silo"

epsr = 64;

materials.epsilon = function (tag, x, y, z)
    if tag == 2 then
        return epsr, 0
    end
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
    local Fy = complex.new(1, 0);
    local Fz = complex.new(0, 0);
    return complex3.new(Fx, Fy, Fz);
end

--function err_enabled(tag)
--    if tag == 2 then
--        return true;
--    end
--
--    return false;
--end

function impedance(tag)
    local epsr = materials.epsilon(tag, 0, 0, 0);
    local mur = materials.mu(tag, 0, 0, 0);

    local eps = epsr * const.eps0;
    local mu = mur * const.mu0;

    return math.sqrt(mu/eps);
end

function propagation_constant(tag)
    local epsr = materials.epsilon(tag, 0, 0, 0);
    local mur = materials.mu(tag, 0, 0, 0);

    local eps = epsr * const.eps0;
    local mu = mur * const.mu0;

    local omega = 2*math.pi*sim.frequency;

    return omega*math.sqrt(mu*eps);
end

local Z1 = impedance(1);
local Z2 = impedance(2);

local gamma_12 = (Z2 - Z1)/(Z2 + Z1);
local tau_12 = 2*Z2/(Z2 + Z1);

print("Reflection coeff at 1-2 interface: " .. gamma_12);
print("Transmission coeff at 1-2 interface: " .. tau_12);

function cexp(e)
    return math.cos(e), math.sin(e);
end

function reference_solution(tag, x, y, z)
    local Fx = complex.new(0, 0);
    local Fz = complex.new(0, 0);

    if tag == 1 then
        local kappa_1 = propagation_constant(1);
        local ref, imf = cexp( -kappa_1*z );
        local reb, imb = cexp( kappa_1*(z-2) )
        local Fy = complex.new( ref + gamma_12*reb, imf + gamma_12*imb);
        return complex3.new(Fx, Fy, Fz);
    end

    if tag == 2 then
        local kappa_1 = propagation_constant(1);
        local kappa_2 = propagation_constant(2);
        local ref, imf = cexp( -kappa_1*1 - kappa_2*(z-1) );
        local Fy = complex.new(tau_12*ref, tau_12*imf);
        return complex3.new(Fx, Fy, Fz);
    end

    local kappa_1 = propagation_constant(1);
    local reb, imb = cexp( kappa_1*(z-2) )
    local Fy = complex.new( gamma_12*reb, gamma_12*imb);
    return complex3.new(Fx, Fy, Fz);
end

local neumann_bcs = {1,2,7,8,12,13}
for i,v in ipairs(neumann_bcs) do
    boundary[v] = {}
    boundary[v].kind = "neumann";
end

local impedance_bcs = {11, 16}
for i,v in ipairs(impedance_bcs) do
    boundary[v] = {}
    boundary[v].kind = "impedance";
end

boundary[5] = {}
boundary[5].kind = "tfsf";
boundary[5].source = plane_wave;

--boundary[11].value = 377/2;

domain = {}
domain[3] = {}
domain[3].scattered_field = true;

frequencies_MHz = {100, 150, 200, 250, 300};
--frequencies_MHz = {300};

function fif(condition, if_true, if_false)
  if condition then return if_true else return if_false end
end

function xxx()
    local fn = "validation_" .. epsr .. ".txt"
    local fh = io.open(fn, "w");
    
    for i,v in ipairs(frequencies_MHz) do
       --silo.filename = "validation_" .. v .. ".silo";
       --save_to_silo();

       sim.frequency = v*1e6;

       hho.classical_stabparam = true
       assemble();
       solve();
       local rl_csp = compute_return_loss(5);
       local err_csp = compute_error();
       print("Error (csp): ", err_csp);
       silo.filename = "csp_" .. v .. ".silo";
       save_to_silo();

       hho.classical_stabparam = false;
       assemble();
       solve();
       local rl_wsp = compute_return_loss(5);
       local err_wsp = compute_error();
       print("Error (wsp): ", err_wsp);
       silo.filename = "wsp_" .. v .. ".silo";
       save_to_silo();
        fh:write(v .. " " .. rl_csp .. " " .. err_csp .. " " .. rl_wsp .. " " .. err_wsp .. "\n")
    end
    fh:close();
end

