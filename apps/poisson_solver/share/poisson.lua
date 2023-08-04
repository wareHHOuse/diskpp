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

function right_hand_side(domain_num, x, y)
    local sx = math.sin(math.pi*x);
    local sy = math.sin(math.pi*y);
    return 2.0*math.pi*math.pi*sx*sy;
end

function diffusion_coefficient(domain_num, x, y)
    return 1;
end

function solution_process()
    assemble();
    solve();
    return 0;
end

--local dp = diffusion_parameter_2D:new(7)
--print(dp)
--dp:entry(1,1,4.2);
--print(dp:entry(1,1));
--print(dp)


        mesh.refinement_level = 4;
        hho.order = 1;
        run()

