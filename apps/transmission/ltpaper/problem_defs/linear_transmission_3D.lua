solver = "pardiso"

--[[ TEST CASE 1: regular function ]]--
testcase_regular_1 = {};

testcase_regular_1.u1 = function (x, y, z)
    return math.cos(math.pi*x) * math.cos(math.pi*y) * math.cos(math.pi*z);
end

testcase_regular_1.u2 = function (x, y, z)
    return math.sin(math.pi*x) * math.sin(math.pi*y) * (z-1.5);
end

testcase_regular_1.grad_u1 = function (x, y, z)
    local gx = - math.pi * math.sin(math.pi*x) * math.cos(math.pi*y) * math.cos(math.pi*z);
    local gy = - math.pi * math.cos(math.pi*x) * math.sin(math.pi*y) * math.cos(math.pi*z);
    local gz = - math.pi * math.cos(math.pi*x) * math.cos(math.pi*y) * math.sin(math.pi*z);
    return gx, gy, gz;
end

testcase_regular_1.grad_u2 = function (x, y, z)
    local gx = math.pi * math.cos(math.pi*x) * math.sin(math.pi*y) * (z-1.5);
    local gy = math.pi * math.sin(math.pi*x) * math.cos(math.pi*y) * (z-1.5);
    local gz = math.sin(math.pi*x) * math.sin(math.pi*y);
    return gx, gy, gz;
end

testcase_regular_1.f1 = function (x, y, z)
    local threepisq = 3.0 * math.pi * math.pi;
    return threepisq * math.cos(math.pi*x) * math.cos(math.pi*y) * math.cos(math.pi*z);
end

testcase_regular_1.f2 = function (x, y, z)
    local twopisq = 2.0 * math.pi * math.pi;
    return twopisq * math.sin(math.pi*x) * math.sin(math.pi*y) * (z-1.5);
end

testcase_regular_1.testcase_name = "Regular 3D";


--[[ TEST CASE 2: singularity ]]--
testcase_singularity = {};

testcase_singularity.omega1 = 3;
testcase_singularity.omega2 = 5;

testcase_singularity.u1 = function (x, y, z)
    local num = 1 - x*x - y*y - z*z;
    local den = (x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2;
    local c = 0.680553933007535;
    return num/den - c;
end

testcase_singularity.u2 = function (x, y, z)
    local twopi = 2*math.pi
    return math.sin(twopi*x) * math.sin(twopi*y) * math.sin(twopi*z)
end

testcase_singularity.grad_u1 = function (x, y, z)
    local gx = -2*x/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2) + (1.1 - 2*x)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2
    local gy = -2*y/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2) + (1.1 - 2*y)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2
    local gz = -2*z/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2) + (1.1 - 2*z)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2
    return gx, gy, gz;
end

testcase_singularity.grad_u2 = function (x, y, z)
    local twopi = 2*math.pi
    local gx = twopi * math.cos(twopi*x) * math.sin(twopi*y) * math.sin(twopi * z)
    local gy = twopi * math.sin(twopi*x) * math.cos(twopi*y) * math.sin(twopi * z)
    local gz = twopi * math.sin(twopi*x) * math.sin(twopi*y) * math.cos(twopi * z)
    return gx, gy, gz
end

testcase_singularity.f1 = function (x, y, z)
    return -(-4*x*(1.1 - 2*x)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2 - 4*y*(1.1 - 2*y)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2 - 4*z*(1.1 - 2*z)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2 + (1.1 - 2*x)*(2.2 - 4*x)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^3 + (1.1 - 2*y)*(2.2 - 4*y)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^3 + (1.1 - 2*z)*(2.2 - 4*z)*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^3 - 6/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2) - 6*(-x^2 - y^2 - z^2 + 1)/((x - 0.55)^2 + (y - 0.55)^2 + (z - 0.55)^2)^2)
end

testcase_singularity.f2 = function (x, y, z)
    local twopi = 2*math.pi
    return 3 * twopi^2 * math.sin(twopi*x) * math.sin(twopi*y) * math.sin(twopi*z)
end

testcase_singularity.testcase_name = "Singularity 3D";

--[[ TEST CASE 3: nonsmooth ]]--
testcase_nonsmooth = {};
testcase_nonsmooth.omega1 = 3;
testcase_nonsmooth.omega2 = 5;

testcase_nonsmooth.u1 = function(x, y, z)
    local rtarg = x*x + y*y + z*z;
    return math.pow(rtarg, 1.0/4.0)-0.968110352555118;
end

testcase_nonsmooth.u2 = function(x, y, z)
    local twopi = 2*math.pi
    return math.sin(twopi*x) * math.sin(twopi*y) * math.sin(twopi*z)
end

testcase_nonsmooth.grad_u1 = function(x, y, z)
    local rtarg = x*x + y*y + z*z;
    local sq = 0.5*math.pow(rtarg, -3.0/4.0);
    local gx = x*sq;
    local gy = y*sq;
    local gz = z*sq;
    return gx, gy, gz
end

testcase_nonsmooth.grad_u2 = function (x, y, z)
    local twopi = 2*math.pi
    local gx = twopi * math.cos(twopi*x) * math.sin(twopi*y) * math.sin(twopi * z)
    local gy = twopi * math.sin(twopi*x) * math.cos(twopi*y) * math.sin(twopi * z)
    local gz = twopi * math.sin(twopi*x) * math.sin(twopi*y) * math.cos(twopi * z)
    return gx, gy, gz
end

testcase_nonsmooth.f1 = function (x, y, z)
    local rtarg = x*x + y*y + z*z;
    return 0.75*math.pow(rtarg, -3.0/4.0);
end

testcase_nonsmooth.f2 = function (x, y, z)
    local twopi = 2*math.pi
    return 3 * twopi^2 * math.sin(twopi*x) * math.sin(twopi*y) * math.sin(twopi*z)
end

testcase_nonsmooth.testcase_name = "Nonsmooth 3D";



--[[ SELECT TEST CASE ]]--
--testcase = testcase_regular_1;
--testcase = testcase_singularity;
testcase = testcase_nonsmooth;

--[[ DO NOT TOUCH ]]--
omega1 = testcase.omega1
omega2 = testcase.omega2
u1 = testcase.u1;
u2 = testcase.u2;
grad_u1 = testcase.grad_u1;
grad_u2 = testcase.grad_u2;
f1 = testcase.f1;
f2 = testcase.f2;
testcase_name = testcase.testcase_name;
