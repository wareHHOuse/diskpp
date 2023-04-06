
sim.order = 2;
sim.mesh_filename = "modeconv.geo3g"


materials.epsilon = function (tag, x, y, z)
    if (tag == 14 or tag == 2) then
        return 1, 0;
    end
    return 24, 0;
end

materials.mu = function (tag, x, y, z)
    return 1, 0
end

materials.sigma = function (tag, x, y, z)
    return 0
end

function waveguide_impedance(freq, a, b, n, m)
    local eps = const.eps0;
    local mu = const.mu0;
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

function source_mode_te10(tag, x, y, z)
    local Ez_re = math.sin(y*math.pi/0.0229);

    local Fx = complex.new(0, 0);
    local Fy = complex.new(0, 0);
    local Fz = complex.new(Ez_re, 0);
    return complex3.new(Fx, Fy, Fz);
end

local impedance_bcs = {1189, 1195}
for i,v in ipairs(impedance_bcs) do
    boundary[v] = {}
    boundary[v].kind = "impedance";
end

boundary[1188] = {}
boundary[1188].kind = "tfsf";
boundary[1188].source = source_mode_te10;


domain = {}
domain[2] = {}
domain[2].scattered_field = true;

freq_start = 7.0e9;
freq_stop = 11.0e9;
freq_step = 100e6;

local range_low = {};
range_low.start = 7.0e9;
range_low.stop = 11.0e9;
range_low.step = 100e6;
range_low.m = 1;
range_low.n = 0;
range_low.name = "low";

local range_high = {};
range_high.start = 13.0e9;
range_high.stop = 17.0e9;
range_high.step = 100e6;
range_high.m = 2;
range_high.n = 0;
range_high.name = "high";

local ranges = {range_low, range_high};
--local ranges = {range_high};

function xxx()
   print("Running order " .. sim.order)

   for iran, range in ipairs(ranges) do
       local data_fn = "s11_order_" .. sim.order .. "_range_" .. range.name .. ".txt";

       fh = io.open(data_fn, "w")

       local freq_start = range.start;
       local freq_stop = range.stop;

       local frequency = freq_start;
       local count = 0;

       while (frequency <= freq_stop) do
           sim.frequency = frequency;

           print("Running simulation at frequency ", frequency);
           boundary[1189].value = waveguide_impedance(frequency, 0.0229, 0.011, range.m, range.n); -- TE20
           boundary[1195].value = waveguide_impedance(frequency, 0.0229, 0.011, 1, 0); -- TE10
           boundary[1188].value = waveguide_impedance(frequency, 0.0229, 0.011, 1, 0);

           assemble();
           solve();
           silo.filename = "modeconv_" .. "r" .. iran .. "o" .. sim.order .. "_" .. count .. ".silo";
           save_to_silo();
           local s11 = compute_return_loss(1188);

           fh:write(frequency, " ", s11, "\n")
           fh:flush()

           frequency = frequency + range.step;
           count = count + 1;
       end
       fh:close();

    end
end

