solver = "pardiso"

--[[ TEST CASE 1: regular function ]]--
testcase_regular_1 = {};
testcase_regular_1.testcase_name = "Regular 2D";

testcase_regular_1.u1 = function (x, y)
    return math.cos(math.pi*x) * math.cos(math.pi*y);
end

testcase_regular_1.u2 = function (x, y)
    return math.sin(math.pi*x) * math.sin(math.pi*y);
end

testcase_regular_1.grad_u1 = function (x, y)
    local gx = - math.pi * math.sin(math.pi*x) * math.cos(math.pi*y);
    local gy = - math.pi * math.cos(math.pi*x) * math.sin(math.pi*y);
    return gx, gy;
end

testcase_regular_1.grad_u2 = function (x, y)
    local gx = math.pi * math.cos(math.pi*x) * math.sin(math.pi*y);
    local gy = math.pi * math.sin(math.pi*x) * math.cos(math.pi*y);
    return gx, gy;
end

testcase_regular_1.f1 = function (x, y)
    local twopisq = 2.0 * math.pi * math.pi;
    return twopisq * math.cos(math.pi*x) * math.cos(math.pi*y);
end

testcase_regular_1.f2 = function (x, y)
    local twopisq = 2.0 * math.pi * math.pi;
    return twopisq * math.sin(math.pi*x) * math.sin(math.pi*y);
end


--[[ TEST CASE 2: singularity ]]--
testcase_singularity = {};
testcase_singularity.testcase_name = "Singularity 2D";

testcase_singularity.u1 = function (x, y)
    local num = x*y;
    local den = (x-0.55)^2 + y^2;
    return num/den; 
end

testcase_singularity.u2 = function (x, y)
    return (x-y)/(x^2 + y^2);
end

testcase_singularity.grad_u1 = function (x, y)
    local gx = (x*y*(1.1 - 2*x))/((y^2 + (x - 0.55)^2)^2) + y/(y^2 + (x - 0.55)^2);
    local gy = (-2*x*y^2)/((y^2 + (x - 0.55)^2)^2) + x/(y^2 + (x - 0.55)^2);
    return gx, gy;
end

testcase_singularity.grad_u2 = function (x, y)
    local gx = -2*x*(x - y)/(x^2 + y^2)^2 + 1/(x^2 + y^2);
    local gy = -2*y*(x - y)/(x^2 + y^2)^2 - 1/(x^2 + y^2);
    return gx, gy;
end

testcase_singularity.f1 = function (x, y)
    local a = 8*x*y^3/(y^2 + (x - 0.55)^2)^3;
    local b = x*y*(1.1 - 2*x)*(2.2 - 4*x)/(y^2 + (x - 0.55)^2)^3;
    local c = - 8*x*y/(y^2 + (x - 0.55)^2)^2;
    local d = 2*y*(1.1 - 2*x)/(y^2 + (x - 0.55)^2)^2;
    return a + b + c + d;
end

testcase_singularity.f2 = function (x, y)
    local a = 8*x^2*(x - y)/(x^2 + y^2)^3;
    local b = - 4*x/(x^2 + y^2)^2;
    local c = 8*y^2*(x - y)/(x^2 + y^2)^3;
    local d = 4*y/(x^2 + y^2)^2;
    local e = - 4*(x - y)/(x^2 + y^2)^2;
    return a + b + c + d + e;
end


--[[ TEST CASE 3: nonsmooth ]]--
testcase_nonsmooth = {};
testcase_nonsmooth.testcase_name = "Nonsmooth 2D";

testcase_nonsmooth.u1 = function(x, y)
    return x*y/(x^2 + y^2) - 5*math.log(5)/14 + 13*math.log(2)/14;
end

testcase_nonsmooth.grad_u1 = function(x, y)
    local gx = -2*x^2*y/(x^2 + y^2)^2 + y/(x^2 + y^2);
    local gy = -2*x*y^2/(x^2 + y^2)^2 + x/(x^2 + y^2);
    return gx, gy;
end

testcase_nonsmooth.f1 = function(x, y)
    return 8*x^3*y/(x^2 + y^2)^3 + 8*x*y^3/(x^2 + y^2)^3 - 12*x*y/(x^2 + y^2)^2;
end

testcase_nonsmooth.omega1 = 10;
testcase_nonsmooth.omega2 = 11;


--[[ SELECT TEST CASE ]]--
testcase = testcase_regular_1;

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

