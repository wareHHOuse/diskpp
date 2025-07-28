local sfe = require("sfexp");

-- cat Domain_6.csv | tail -n +2 | sed 's/\(.*\);\(.*\);/points:add( sfe.point:new(\1,\2) )/g'

function concave4()
    local points = make_point_vector()
    points:add( sfe.point:new(0.5, 0.0) )
    points:add( sfe.point:new(1.3, 0.75) )
    points:add( sfe.point:new(1.0, 0.0) )
    points:add( sfe.point:new(1.3, -0.75) )
    return points
end

function concave5()
    local points = make_point_vector()
    points:add( sfe.point:new(1.0606601717798214e+00,1.0606601717798212e+00) )
    points:add( sfe.point:new(3.0616169978683830e-17,5.0000000000000000e-01) )
    points:add( sfe.point:new(-1.0606601717798212e+00,1.0606601717798214e+00) )
    points:add( sfe.point:new(-2.5000000000000022e-01,-4.3301270189221919e-01) )
    points:add( sfe.point:new(2.5000000000000006e-01,-4.3301270189221930e-01) )
    return points
end

function concave6()
    local points = make_point_vector()
    points:add( sfe.point:new(1.5000000000000000e+00,0.0000000000000000e+00) )
    points:add( sfe.point:new(2.5000000000000006e-01,4.3301270189221930e-01) )
    points:add( sfe.point:new(-7.4999999999999967e-01,1.2990381056766580e+00) )
    points:add( sfe.point:new(-5.0000000000000000e-01,6.1232339957367660e-17) )
    points:add( sfe.point:new(-7.5000000000000067e-01,-1.2990381056766576e+00) )
    points:add( sfe.point:new(2.4999999999999967e-01,-4.3301270189221952e-01) )
    return points
end

function concave7()
    local points = make_point_vector()
    points:add( sfe.point:new(1.0000000000000000e+00,0.0000000000000000e+00) )
    points:add( sfe.point:new(3.4202014332566882e-01,9.3969262078590832e-01) )
    points:add( sfe.point:new(-1.9101299543362335e-01,1.0832885283134288e+00) )
    points:add( sfe.point:new(-2.4999999999999989e-01,4.3301270189221935e-01) )
    points:add( sfe.point:new(-1.8000000000000000e+00,2.2043642384652358e-16) )
    points:add( sfe.point:new(-1.9101299543362338e-01,-1.0832885283134288e+00) )
    points:add( sfe.point:new(6.9282032302755070e-01,-4.0000000000000036e-01) )
    return points
end

function concave8()
    local points = make_point_vector()
    points:add( sfe.point:new(1.5000000000000000e+00,0.0000000000000000e+00) )
    points:add( sfe.point:new(3.5355339059327379e-01,3.5355339059327373e-01) )
    points:add( sfe.point:new(9.1848509936051484e-17,1.5000000000000000e+00) )
    points:add( sfe.point:new(-3.5355339059327373e-01,3.5355339059327379e-01) )
    points:add( sfe.point:new(-1.5000000000000000e+00,1.8369701987210297e-16) )
    points:add( sfe.point:new(-3.5355339059327384e-01,-3.5355339059327373e-01) )
    points:add( sfe.point:new(-2.7554552980815448e-16,-1.5000000000000000e+00) )
    points:add( sfe.point:new(3.5355339059327368e-01,-3.5355339059327384e-01) )
    return points
end

function concave9()
    local points = make_point_vector()
    points:add( sfe.point:new(9.2893050912796538e-01,9.2322275130987996e-01) )
    points:add( sfe.point:new(1.6132706602562455e-01,6.0017494279869843e-01) )
    points:add( sfe.point:new(-3.3511200065087798e-01,8.5890271820246022e-01) )
    points:add( sfe.point:new(-5.1785603009758219e-01,2.2938253790856930e-01) )
    points:add( sfe.point:new(-1.0857312392458165e+00,-3.5715494120122565e-01) )
    points:add( sfe.point:new(-2.7236996293005655e-01,-7.4765111568082843e-01) )
    points:add( sfe.point:new(1.4258574776721900e-01,-5.5152266323297305e-01) )
    points:add( sfe.point:new(8.0137781025594423e-01,-3.6881802466479097e-01) )
    points:add( sfe.point:new(8.2271009605504108e-01,5.2920749583350120e-16) )
    return points
end

function make_poly(n)
    local points = make_point_vector()

    for i = 0,n-1 do
        local x = math.cos(i*2*math.pi/n)
        local y = math.sin(i*2*math.pi/n)
        points:add( sfe.point:new(x,y) )
    end

    return points;
end

cfg = sfe.config:new();
cfg.vertex = 1;
cfg.degree = 1;
cfg.variant = 0 
cfg.silo_fn = "explore.silo"
cfg.use_stabfree = true
cfg.eps = 0.04
cfg.nsamples = 201 
cfg.multi_thread = true

print(cfg)

--local poly = make_point_vector();
--poly:add( sfe.point:new(0.0, 0.0) )
--poly:add( sfe.point:new(0.7, 0.0) )
--poly:add( sfe.point:new(1.0, 1.0) )
--poly:add( sfe.point:new(0.2, 0.8) )
--poly:add( sfe.point:new(-0.3, 1.0) )

local poly = concave9()
run_poly(cfg, poly)

--local poly2 = make_poly(6)
--run_poly(cfg, poly2)

--cfg.num_vertices = 8
--run(cfg)
