clear;
clf;

rows = 2;
cols = 4;

ymin = 1e-7;
ymax = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o2/bicgstab_eo.txt");
bicgstab_k1 = load("l1_s4_k1_o2/bicgstab_eo.txt");
bicgstab_k2 = load("l1_s4_k2_o2/bicgstab_eo.txt");
bicgstab_k3 = load("l1_s4_k3_o2/bicgstab_eo.txt");

subplot(rows,cols,1);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Equal order, overlap = 2")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o4/bicgstab_eo.txt");
bicgstab_k1 = load("l1_s4_k1_o4/bicgstab_eo.txt");
bicgstab_k2 = load("l1_s4_k2_o4/bicgstab_eo.txt");
bicgstab_k3 = load("l1_s4_k3_o4/bicgstab_eo.txt");

subplot(rows,cols,2);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Equal order, overlap = 4")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o8/bicgstab_eo.txt");
bicgstab_k1 = load("l1_s4_k1_o8/bicgstab_eo.txt");
bicgstab_k2 = load("l1_s4_k2_o8/bicgstab_eo.txt");
bicgstab_k3 = load("l1_s4_k3_o8/bicgstab_eo.txt");

subplot(rows,cols,3);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Equal order, overlap = 8")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o16/bicgstab_eo.txt");
bicgstab_k1 = load("l1_s4_k1_o16/bicgstab_eo.txt");
bicgstab_k2 = load("l1_s4_k2_o16/bicgstab_eo.txt");
bicgstab_k3 = load("l1_s4_k3_o16/bicgstab_eo.txt");

subplot(rows,cols,4);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Equal order, overlap = 16")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o2/bicgstab_mo.txt");
bicgstab_k1 = load("l1_s4_k1_o2/bicgstab_mo.txt");
bicgstab_k2 = load("l1_s4_k2_o2/bicgstab_mo.txt");
bicgstab_k3 = load("l1_s4_k3_o2/bicgstab_mo.txt");

subplot(rows,cols,5);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Mixed order, overlap = 2")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o4/bicgstab_mo.txt");
bicgstab_k1 = load("l1_s4_k1_o4/bicgstab_mo.txt");
bicgstab_k2 = load("l1_s4_k2_o4/bicgstab_mo.txt");
bicgstab_k3 = load("l1_s4_k3_o4/bicgstab_mo.txt");

subplot(rows,cols,6);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Mixed order, overlap = 4")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o8/bicgstab_mo.txt");
bicgstab_k1 = load("l1_s4_k1_o8/bicgstab_mo.txt");
bicgstab_k2 = load("l1_s4_k2_o8/bicgstab_mo.txt");
bicgstab_k3 = load("l1_s4_k3_o8/bicgstab_mo.txt");

subplot(rows,cols,7);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Mixed order, overlap = 8")
ylim([ymin,ymax]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bicgstab_k0 = load("l1_s4_k0_o16/bicgstab_mo.txt");
bicgstab_k1 = load("l1_s4_k1_o16/bicgstab_mo.txt");
bicgstab_k2 = load("l1_s4_k2_o16/bicgstab_mo.txt");
bicgstab_k3 = load("l1_s4_k3_o16/bicgstab_mo.txt");

subplot(rows,cols,8);
semilogy(bicgstab_k0(:,1));
hold on;
semilogy(bicgstab_k1(:,1));
semilogy(bicgstab_k2(:,1));
semilogy(bicgstab_k3(:,1));
title("Mixed order, overlap = 16")
ylim([ymin,ymax]);


