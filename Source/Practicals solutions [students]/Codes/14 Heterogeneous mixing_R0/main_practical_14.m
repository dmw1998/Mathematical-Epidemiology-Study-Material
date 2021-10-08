clear all; clc;
%% Part 1

% R0 = beta*N*D
WAIFW = 4.16e-6*ones(3,3);
WAIFW(1,1) = 1.66e-5;

D = 11;

N_y = 15000;
N_m = 15000;
N_o = 30000;

N = [N_y; N_m; N_o];

%P1
R_yy = WAIFW(1,1)*N_y*D;
R_ym = WAIFW(1,2)*N_y*D;
R_yo = WAIFW(1,3)*N_y*D;

% P2

NGM = zeros(3,3);

for i = 1:3
    for j = 1:3
        NGM(i,j) = WAIFW(i,j)*N(i)*D;
    end
end

G2y = sum(NGM*[1;0;0]);
G2m = sum(NGM*[0;1;0]);
G2o = sum(NGM*[0;0;1]);

%% Part 2

% P 1 ~ 4 one infectious young individual

Ini_inf = [1;0;0];

[G, SG, R] = Gen(Ini_inf, NGM);

% P1
py = G(1,1)/SG(1);
pm = G(2,1)/SG(1);
po = G(3,1)/SG(1);

%P2
T1 = SG(1);

%P3
figure(1)
hold on;
plot(1:10, G(1,:)./SG)
plot(1:10, G(2,:)./SG)
plot(1:10, G(3,:)./SG)
legend('young', 'middle', 'old')
hold off;

%P4
%R

%P5
Ini_inf2 = [20; 50; 30];
Ini_inf3 = [0.5; 0.2; 0.3];
[G2, SG2, R2] = Gen(Ini_inf2, NGM);

figure(2)
hold on;
plot(1:10, G2(1,:)./SG2)
plot(1:10, G2(2,:)./SG2)
plot(1:10, G2(3,:)./SG2)
legend('young', 'middle', 'old')
hold off;
title('[20; 50; 30]')

[G3, SG3, R3] = Gen(Ini_inf3, NGM);

figure(3)
hold on;
plot(1:10, G3(1,:)./SG3)
plot(1:10, G3(2,:)./SG3)
plot(1:10, G3(3,:)./SG3)
legend('young', 'middle', 'old')
hold off;
title('[0.5; 0.2; 0.3]')

R2
R3

% P6

[V,D] = eig(NGM);
max_eig = D(1,1);
eig_vec = V(:,1)./sum(V(:,1));


% P7
Ini_inf2 = [20; 50; 30];
figure(4)
title('initial = [20; 50; 30]')
hold on;
for c = [0.25 0.5 0.725 0.75]
     [~,SG,~]=Gen1(Ini_inf2,NGM,c);
    
     plot(1:10, SG);
end
legend('25%', '50%', "72.5%", '75%')
hold off

figure(5)
[~,SG_25,~]=Gen1(Ini_inf2,NGM,0.25);
plot(1:10, SG_25);
title('25% immune')
figure(6)
[~,SG_50,~]=Gen1(Ini_inf2,NGM,0.5);
plot(1:10, SG_50);
title('50% immune')
figure(7)
[~,SG_725,~]=Gen1(Ini_inf2,NGM,0.725);
plot(1:10, SG_725);
title('72.5% immune')
figure(8)
[~,SG_75,~]=Gen1(Ini_inf2,NGM,0.75);
plot(1:10, SG_75);
title('75% immune')

[~,~,R_25]=Gen1(Ini_inf2,NGM,0.25)
[~,~,R_50]=Gen1(Ini_inf2,NGM,0.50)
[~,~,R_725]=Gen1(Ini_inf2,NGM,0.725)
[~,~,R_75]=Gen1(Ini_inf2,NGM,0.75)