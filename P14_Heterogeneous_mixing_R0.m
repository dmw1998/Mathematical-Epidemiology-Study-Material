clear all; close all; clc;
%% Part I
% Given data
D = 11;					    % Infectious period 11 days

WAIFW_B = 4.16e-6 * ones(3,3);
WAIFW_B(1,1) = 1.66e-5;

N = [15000; 15000; 30000];      % N = [N_y; N_m, N_o]

% P1
% R_0 = beta*N*D
R_y = WAIFW_B(1,:)*N(1)*D;                  % R_y = [R_yy, R_ym, R_yo]
R_m = WAIFW_B(2,:)*N(2)*D;                  % R_m = [R_my, R_mm, R_mo]
R_o = WAIFW_B(3,:)*N(3)*D;                  % R_o = [R_oy, R_om, R_oo]

NGM = [R_y; R_m; R_o];                      % Next Generated Matrix

% P2
G_y = sum(NGM*[1;0;0]);                     % One infectious young person
G_m = sum(NGM*[0;1;0]);                     % One infectious middle-aged person
G_o = sum(NGM*[0;0;1]);                     % One infectious old person

I_2nd = G_y + G_m + G_o;                    % Total

%% Part II
% P1
prop_I_1st = NGM*[1;0;0]/sum(NGM*[1;0;0]);

% P2
G_y;

% P3
k = 10;                         % Number of iterations
I = [1; 0; 0];                  % Initial infectious is one infectious young person
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
for i = 1 : k                   % For the first 10 interations
    I = NGM*I;                  % I = [I_y; I_m; I_o]
    G_i = sum(I);
    data_coll_G(i) = G_i;
    data_coll_I(:,i) = I;
    data_coll_distr(:,i) = I/G_i;
end

figure
x = 1:10;
plot(x, data_coll_distr);
legend('Young','Middle-aged','Old','Location', 'best')
title('Distribution of the infectious persons in 10 generations (one young infectious)')
% The distribution converges to y:m:o = 0.4336:0.1888:0.3776

% P4
G_y1 = G_i;                     % Total secondary infectious with one young initial infectious after 10 generations

I = [0; 1; 0];                  % Initial infectious is one infectious middle-aged person
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
for i = 1 : k                   % For the first 10 interations
    I = NGM*I;                  % I = [I_y; I_m; I_o]
    G_i = sum(I);
    data_coll_G(i) = G_i;
    data_coll_I(:,i) = I;
    data_coll_distr(:,i) = I/G_i;
end
G_m1 = G_i;                     % Total secondary infectious with one middle-aged initial infectious after 10 generations

I = [0; 0; 1];                  % Initial infectious is one infectious old person
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
for i = 1 : k                   % For the first 10 interations
    I = NGM*I;                  % I = [I_y; I_m; I_o]
    G_i = sum(I);
    data_coll_G(i) = G_i;
    data_coll_I(:,i) = I;
    data_coll_distr(:,i) = I/G_i;
end
G_o1 = G_i;                     % Total secondary infectious with one old initial infectious after 10 generations

avd_G = (G_y1+G_m1+G_o1)/3;     % Average of secondary infectious persons resulting from each infectious person after 10 generations

% P5
I = [20; 50; 30];               % Initial infectious
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
for i = 1 : k                   % For the first 10 interations
    I = NGM*I;                  % I = [I_y; I_m; I_o]
    G_i = sum(I);
    data_coll_G(i) = G_i;
    data_coll_I(:,i) = I;
    data_coll_distr(:,i) = I/G_i;
end

figure
x = 1:10;
subplot(2,1,1);
plot(x, data_coll_distr);
ax = gca;
ax.YLim = [0.17 0.47];
legend('Young','Middle-aged','Old','Location', 'best')
title('Distribution of the infectious persons in 10 generations (20,50,30)')

I = [0.5; 0.2; 0.3];            % Initial infectious
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
for i = 1 : k                   % For the first 10 interations
    I = NGM*I;                  % I = [I_y; I_m; I_o]
    G_i = sum(I);
    data_coll_G(i) = G_i;
    data_coll_I(:,i) = I;
    data_coll_distr(:,i) = I/G_i;
end

subplot(2,1,2);
plot(x, data_coll_distr);
ax = gca;
ax.YLim = [0.17 0.47];
legend('Young','Middle-aged','Old','Location', 'best')
title('Distribution of the infectious persons in 10 generations (0.5,0.2,0.3)')

% As two graphs shown, even they start from different initial values, they converge to the same ratio
% The answer to P3 and P4 does not change.

% P6
[V,D] = eig(NGM,'vector');          % returns vector D of eigenvalues and matrix V whose columns are the corresponding eigenvectors
[val, loc] = max(abs(D));                    % returns maximal value val and index loc
max_eigv = V(:,loc);                    % corresponding eigenvector of the maximal eigenvalue
% val \approx data_coll_G(k)/data_coll_G(k-1)

% P7
I = [1; 0; 0];                  % Initial infectious is one infectious young person
data_coll_distr = [];           % Each column is the ratio of I_y:I_m:I_o
data_coll_I = [];               % Each column is [I_y, I_m, I_o] for each time t
data_coll_G = [];               % Each component is the total infectious for each time t
c = [0.25, 0.5, 0.725, 0.75];
figure
for j = 1 : 4
    for i = 1 : k                   % For the first 10 interations
        I = (1-c(j))*NGM*I;                  % I = [I_y; I_m; I_o]
        G_i = sum(I);
        data_coll_G(i) = G_i;
        data_coll_I(:,i) = I;
        data_coll_distr(:,i) = I/G_i;
    end

    subplot(2,2,j)
    plot(x, data_coll_distr);
    legend('Young','Middle-aged','Old','Location', 'best')
    title('Distribution of the infectious persons with ')
end