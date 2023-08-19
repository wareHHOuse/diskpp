mesh.source = "internal";
mesh.type = "triangles";

hho.variant = "equal_order";
hho.use_stabfree = false;

boundary[0] = "dirichlet";
boundary[1] = "dirichlet";
boundary[2] = "dirichlet";
boundary[3] = "dirichlet";


function dirichlet_data(bnd_num, x, y)
   
    return 0;
end

function neumann_data(bnd_num, x, y)
    return 0;
end

function solution(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    return sx*sy;
end

function right_hand_side(domain_num, x, y, z)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    local sz = 1.0;
    if sim.dimension == 3 then
        sz = math.sin(math.pi*z);
    end
    return sim.dimension*math.pi*math.pi*sx*sy*sz;
end

function diffusion_coefficient(domain_num, x, y)
    return 1;
end

function solution_process()
    local tc = timecounter:new();

    --print(sim.dimension);
    
    tc:tic();
    assemble();
    tc:toc();
    print(tc);

    solve();

    export_to_visit();

    check();
    return 0;
end

--local dp = diffusion_parameter_2D:new(7)
--print(dp)
--dp:entry(1,1,4.2);
--print(dp:entry(1,1));
--print(dp)


mesh.refinement_level = 5;
hho.order = 0;

for i = 2,5 do
    mesh.refinement_level = i;
    run()
end

