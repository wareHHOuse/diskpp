function epsilon(tag)
    if (tag == 14 or tag == 2) then
        return 1, 0;
    end
    return 24, 0;
end

function mu(tag)
    return 1, 0;
end

function sigma(tag)
    return 0;
end
--[[
function pw_source_1183(tag, x, y, z)
    local field = E_complex_field.new();
    field.Ez_re = 1.0;
    return field;
end

function e_field_1196(tag, x, y, z)
    local field = E_complex_field.new();
    field.Ez_re = math.sin(y*math.pi/0.0229);
    return field;
end


boundary[1196] = {};
boundary[1196].kind = "plane_wave"
boundary[1196].source = pw_source_1183;

--boundary[1196] = {};
--boundary[1196].kind = "e_field";
--boundary[1196].source = e_field_1196;

boundary[1205] = {};
boundary[1205].kind = "impedance";
boundary[1205].value = 700;

boundary[1208] = {};
boundary[1208].kind = "impedance";
boundary[1208].value = 700;
--]]

function waveguide_impedance(freq, a, b, n, m)
    local eps = 8.85e-12;
    local mu = 4*math.pi*1e-7;
    local Z0 = math.sqrt(mu/eps);
    local c0 = 1/math.sqrt(mu*eps);
    local k0 = 2*math.pi*freq*math.sqrt(mu*eps);
    local lambda0 = c0/freq;
    local aaa = 2*math.pi/lambda0;
    local bbb = n*math.pi/a;
    local ccc = m*math.pi/b;
    local beta = math.sqrt(aaa*aaa - bbb*bbb - ccc*ccc);
    print(freq, " ", n, " ", m , " ", k0*Z0/beta);
    return k0*Z0/beta;
end

order = 2


boundary[1189] = {};
boundary[1189].kind = "impedance";

boundary[1195] = {};
boundary[1195].kind = "impedance";

freq_start = 7.0e9;
freq_stop = 11.0e9;
freq_step = 100e6;

frequency = freq_start;
count = 1;

fh = io.open("s11.txt", "w")

while (frequency <= freq_stop) do
    print("Running simulation at frequency ", frequency);
    boundary[1189].value = waveguide_impedance(frequency, 0.0229, 0.011, 1, 0); -- TE20
    boundary[1195].value = waveguide_impedance(frequency, 0.0229, 0.011, 1, 0); -- TE10
    wgz = waveguide_impedance(frequency, 0.0229, 0.011, 1, 0);
    silo_fn = "maxwell_" .. count .. ".silo";

    run_solver()

    fh:write(frequency, " ", s11, "\n")
    fh:flush()

    frequency = frequency + freq_step;
    count = count + 1;
end

fh:close()

