clear;
clf;

rows = 2;
cols = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o2/error_eo.txt");
error_k1 = load("l1_s16_k1_o2/error_eo.txt");
error_k2 = load("l1_s16_k2_o2/error_eo.txt");
error_k3 = load("l1_s16_k3_o2/error_eo.txt");

subplot(rows,cols,1);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Equal order, overlap = 2")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o4/error_eo.txt");
error_k1 = load("l1_s16_k1_o4/error_eo.txt");
error_k2 = load("l1_s16_k2_o4/error_eo.txt");
error_k3 = load("l1_s16_k3_o4/error_eo.txt");

subplot(rows,cols,2);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Equal order, overlap = 4")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o8/error_eo.txt");
error_k1 = load("l1_s16_k1_o8/error_eo.txt");
error_k2 = load("l1_s16_k2_o8/error_eo.txt");
error_k3 = load("l1_s16_k3_o8/error_eo.txt");

subplot(rows,cols,3);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Equal order, overlap = 8")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o16/error_eo.txt");
error_k1 = load("l1_s16_k1_o16/error_eo.txt");
error_k2 = load("l1_s16_k2_o16/error_eo.txt");
error_k3 = load("l1_s16_k3_o16/error_eo.txt");

subplot(rows,cols,4);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Equal order, overlap = 16")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o2/error_mo.txt");
error_k1 = load("l1_s16_k1_o2/error_mo.txt");
error_k2 = load("l1_s16_k2_o2/error_mo.txt");
error_k3 = load("l1_s16_k3_o2/error_mo.txt");

subplot(rows,cols,5);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Mixed order, overlap = 2")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o4/error_mo.txt");
error_k1 = load("l1_s16_k1_o4/error_mo.txt");
error_k2 = load("l1_s16_k2_o4/error_mo.txt");
error_k3 = load("l1_s16_k3_o4/error_mo.txt");

subplot(rows,cols,6);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Mixed order, overlap = 4")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o8/error_mo.txt");
error_k1 = load("l1_s16_k1_o8/error_mo.txt");
error_k2 = load("l1_s16_k2_o8/error_mo.txt");
error_k3 = load("l1_s16_k3_o8/error_mo.txt");

subplot(rows,cols,7);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Mixed order, overlap = 8")
ylim([1e-8,10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_k0 = load("l1_s16_k0_o16/error_mo.txt");
error_k1 = load("l1_s16_k1_o16/error_mo.txt");
error_k2 = load("l1_s16_k2_o16/error_mo.txt");
error_k3 = load("l1_s16_k3_o16/error_mo.txt");

subplot(rows,cols,8);
semilogy(error_k0(:,1));
hold on;
semilogy(error_k1(:,1));
semilogy(error_k2(:,1));
semilogy(error_k3(:,1));
title("Mixed order, overlap = 16")
ylim([1e-8,10]);


