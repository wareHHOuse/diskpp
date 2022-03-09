// Mode converter GMSH model

SetFactory("OpenCASCADE");

Mesh.Algorithm = 1;
Mesh.Algorithm3D = 1;

cl = 0.005;

x_base = 0;
y_base = 0;
z_base = 0;

wg_xlen = 0.11;
wg_ylen = 0.0229;
wg_zlen = 0.01;

rods_start_x = 0.10;

rod_spacing = 0.00954;

rod1_d = 0.01145;  rod1_r = 0.000380;  rod1_x = x_base + rods_start_x - rod_spacing*0;  rod1_y = y_base + rod1_d;   rod1_z = z_base;
rod2_d = 0.01145;  rod2_r = 0.000550;  rod2_x = x_base + rods_start_x - rod_spacing*1;  rod2_y = y_base + rod2_d;   rod2_z = z_base;
rod3_d = 0.01145;  rod3_r = 0.000550;  rod3_x = x_base + rods_start_x - rod_spacing*2;  rod3_y = y_base + rod3_d;   rod3_z = z_base;
rod4_d = 0.01145;  rod4_r = 0.000550;  rod4_x = x_base + rods_start_x - rod_spacing*3;  rod4_y = y_base + rod4_d;   rod4_z = z_base;
rod5_d = 0.01045;  rod5_r = 0.000515;  rod5_x = x_base + rods_start_x - rod_spacing*4;  rod5_y = y_base + rod5_d;   rod5_z = z_base;
rod6_d = 0.00900;  rod6_r = 0.000570;  rod6_x = x_base + rods_start_x - rod_spacing*5;  rod6_y = y_base + rod6_d;   rod6_z = z_base;
rod7_d = 0.00700;  rod7_r = 0.000640;  rod7_x = x_base + rods_start_x - rod_spacing*6;  rod7_y = y_base + rod7_d;   rod7_z = z_base;
rod8_d = 0.00500;  rod8_r = 0.000690;  rod8_x = x_base + rods_start_x - rod_spacing*7;  rod8_y = y_base + rod8_d;   rod8_z = z_base;
rod9_d = 0.00300;  rod9_r = 0.000760;  rod9_x = x_base + rods_start_x - rod_spacing*8;  rod9_y = y_base + rod9_d;   rod9_z = z_base;
rod10_d = 0.00125; rod10_r = 0.000920; rod10_x = x_base + rods_start_x - rod_spacing*9; rod10_y = y_base + rod10_d; rod9_z = z_base;
rod11_d = 0.02165; rod11_r = 0.000960; rod11_x = x_base + rods_start_x - rod_spacing*9; rod11_y = y_base + rod11_d; rod9_z = z_base;

cl_rod = 0.0005;
xb = 100; rx = rod1_x; ry = rod1_y; rz = z_base; rr = rod1_r; Include "modeconv.i1";
xb = 200; rx = rod2_x; ry = rod2_y; rz = z_base; rr = rod2_r; Include "modeconv.i1";
xb = 300; rx = rod3_x; ry = rod3_y; rz = z_base; rr = rod3_r; Include "modeconv.i1";
xb = 400; rx = rod4_x; ry = rod4_y; rz = z_base; rr = rod4_r; Include "modeconv.i1";
xb = 500; rx = rod5_x; ry = rod5_y; rz = z_base; rr = rod5_r; Include "modeconv.i1";
xb = 600; rx = rod6_x; ry = rod6_y; rz = z_base; rr = rod6_r; Include "modeconv.i1";
xb = 700; rx = rod7_x; ry = rod7_y; rz = z_base; rr = rod7_r; Include "modeconv.i1";
xb = 800; rx = rod8_x; ry = rod8_y; rz = z_base; rr = rod8_r; Include "modeconv.i1";
xb = 900; rx = rod9_x; ry = rod9_y; rz = z_base; rr = rod9_r; Include "modeconv.i1";
xb = 1000; rx = rod10_x; ry = rod10_y; rz = z_base; rr = rod10_r; Include "modeconv.i1";
xb = 1100; rx = rod11_x; ry = rod11_y; rz = z_base; rr = rod11_r; Include "modeconv.i1";

Point(1) = {           x_base,           y_base, z_base, cl/3 };
Point(2) = { x_base + wg_xlen,           y_base, z_base, cl };
Point(3) = { x_base + wg_xlen, y_base + wg_ylen, z_base, cl };
Point(4) = {           x_base, y_base + wg_ylen, z_base, cl/3 };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(5) = {1,2,3,4};

Plane Surface(6) = {5};

Extrude {0, 0, wg_zlen} {
    Surface { 6, 121, 221, 321, 421, 521, 621, 721, 821, 921, 1021, 1121 };
    //Surface {6};
}

BooleanDifference(100) = { Volume{1}; Delete; } { Volume{2}; };
BooleanDifference(101) = { Volume{100}; Delete; } { Volume{3}; };
BooleanDifference(102) = { Volume{101}; Delete; } { Volume{4}; };
BooleanDifference(103) = { Volume{102}; Delete; } { Volume{5}; };
BooleanDifference(104) = { Volume{103}; Delete; } { Volume{6}; };
BooleanDifference(105) = { Volume{104}; Delete; } { Volume{7}; };
BooleanDifference(106) = { Volume{105}; Delete; } { Volume{8}; };
BooleanDifference(107) = { Volume{106}; Delete; } { Volume{8}; };
BooleanDifference(108) = { Volume{107}; Delete; } { Volume{9}; };
BooleanDifference(109) = { Volume{108}; Delete; } { Volume{10}; };
BooleanDifference(110) = { Volume{109}; Delete; } { Volume{11}; };
BooleanDifference(111) = { Volume{110}; Delete; } { Volume{12}; };

