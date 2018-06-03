sol = [1       1           1       1
0.975	0.81581     0.9875	0.72838
0.95	0.64279     0.975	0.49123
0.9     0.34121     0.95	0.13617
0.85	0.1116      0.925	-0.01946
0.8     -0.03834	0.9     -0.03123
0.75	-0.10703	0.85	-0.04142
0.7     -0.14267	0.8     -0.05148
0.65	-0.1781     0.75	-0.06148
0.6     -0.20606	0.7     -0.07142
0.55	-0.20679	0.65	-0.08125
0.5     -0.18898	0.6     -0.0908
0.45	-0.16033	0.575	-0.09275
0.4 	-0.12674	0.55	-0.08027
0.35	-0.09261	0.525	-0.05286
0.3     -0.0612     0.5     -0.02406
0.25	-0.03485	0.45	-0.00137
0.2     -0.01523	0.4     -0.00069
0.15	-0.00356	0.3     -0.00038
0.1     -0.00058	0.2     -0.00021
0.05	-0.00027	0.1     -0.0001
0           0       0       0 ];


figure(11)
hold on;
e_y = load('plot_over_y_driven_k0_a1_Bi2mu1_J_50_50.data');
plot(e_y(:,2), e_y(:,3))
e_y = load('plot_over_y_driven_k0_a1_Bi2mu1_J_50_50.data');
plot(e_y(:,2), e_y(:,3))
e_y = load('plot_over_y_driven_k0_a1_Bi50mu1_J_50_50.data');
plot(e_y(:,2), e_y(:,3))


e_y = load('plot_over_y_driven_k0_a1_Bi2mu1_J_tri_h00125.data');
plot(e_y(:,2), e_y(:,3))
e_y = load('plot_over_y_driven_k0_a1_Bi2mu2_J_tri_h00125.data');
plot(e_y(:,2), e_y(:,3))
e

plot(sol(:,1), sol(:,2), '*');
set(gca,'box','on'); set(gcf,'color','w');
set(gca,'fontsize',20);
xlabel('x','FontSize',20);
ylabel('u','FontSize',20);

figure(12)

e_y1 = load("plot_over_y_driven_k0_a1_Bi2mu1_20_20.data"); 
e_y2 = load("plot_over_y_driven_k0_a1_Bi2mu1_40_40.data"); 
e_y3 = load("plot_over_y_driven_k0_a1_Bi2mu1_80_80.data");
e_y4 = load("plot_over_y_driven_k0_a1_Bi2mu1_100_100.data");

plot(e_y1(:,2), e_y1(:,3), 'b',e_y2(:,2), e_y2(:,3), 'k', e_y3(:,2), e_y3(:,3),'r', e_y4(:,2), e_y4(:,3),'g');

legend("20x20","40x40", "80x80","100x100");

figure(14)
e_y1 = load("plot_over_y_driven_k0_a1_Bi2mu1_80_80.data");
e_y2 = load("plot_over_y_driven_k0_a1_Bi2mu2_80_80.data");
e_y3 = load("plot_over_y_driven_k0_a1_Bi2mu1_tol-10_100_100.data");
e_y4 = load("plot_over_y_driven_k0_a1_Bi2mu2_tol-10_100_100.data");

plot(e_y1(:,2), e_y1(:,3), 'b',e_y2(:,2), e_y2(:,3), 'k', e_y3(:,2), e_y3(:,3),'r', e_y4(:,2), e_y4(:,3),'g');

legend("\mu = 1, \sigma_0 = 2","\mu = 2, \sigma_0 = 4", "\mu = 1, \sigma_0 = 2","\mu = 2, \sigma_0 = 4");






CoordsFileName = sprintf('Coords_N%d_A%d_R%d.data', mesh, alpha, k);
    Coords = load("Coords_driven_k0_a1_Bi2mu2_tol-10_100_100.data");
    sc = size(Coords);
    NROWS = sc(2); NCELLS = sc(1);
    NV = NROWS/2;

    X = [Coords(:, 1 : 2 ) Coords(:, 4) Coords(:, 3)];
    Y = Coords(:, NV + 1: 2*NV);

 
    fig = figure;
    F = load("criterion_theta_driven_k0_a1_Bi2mu2_tol-10_100_100.data");   
    F(1, :) = [];
    hp = patch(X', Y', F','EdgeAlpha','none');
        
    figure;
    F = load("criterion_sigma_driven_k0_a1_Bi2mu2_tol-10_100_100.data");
    F(1, :) = [];
    hp = patch(X', Y', F','EdgeAlpha','none');