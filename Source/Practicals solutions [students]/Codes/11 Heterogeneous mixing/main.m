clear all; clc;

%% Part 1: Calculating WAIFW matrices

%% 1. Fill in the parameter values

N0 = 60000;
N_y = 15000;
N_m = 15000;
N_o = 30000;
b = 1/60/365;
d = 1/30/365;
a_y = 1/15/365;
a_m = 1/15/365;
f = 1/10;
r = 1/11;

%% 2. Calculate Number of infectious
% I = lambda * S / r

lambda_y = 0.000364;
lambda_m = 0.000114;
lambda_o = 0.000114;
S_y0 = 5012;
S_m0 = 3083;
S_o0 = 2739;

I_y0 = lambda_y * S_y0 / r;
I_m0 = lambda_m * S_m0 / r;
I_o0 = lambda_o * S_o0 / r;

%% 3. Write down the expressions for the force of infection (FOI)


%% 4. Calculate appropriate contact parameters

WAIFW_A = zeros(3,3);

beta1 = lambda_y / I_y0;
beta2 = lambda_m / I_m0;
beta3 = lambda_o / I_o0;

WAIFW_A(1,1) = beta1;
WAIFW_A(2,2) = beta2;
WAIFW_A(3,3) = beta3;

WAIFW_B = zeros(3,3);

beta3 = lambda_o / (I_y0+I_m0+I_o0);
beta2 = (lambda_m - I_o0 * beta3) / (I_y0+I_m0);
beta1 = (lambda_y - I_m0 * beta2 - I_o0 * beta3) / I_y0;

WAIFW_B(1:3,1:3) = beta3;
WAIFW_B(1:2,1:2) = beta2;
WAIFW_B(1,1) = beta1;

beta = 1.0530e-05;
WAIFW_H = beta*ones(3,3);

%% Part 2 : The implications of heterogeneous mixing

S_y0 = 5000;
E_y0 = 20;
I_y0 = 20;
R_y0 = N_y-S_y0-I_y0;
S_m0 = 3000;
E_m0 = 3;
I_m0 = 3;
R_m0 = N_m-S_m0-I_m0;
S_o0 = 2700;
E_o0 = 3;
I_o0 = 3;
R_o0 = N_o-S_o0-I_o0;

y0 = [S_y0;E_y0;I_y0;R_y0;S_m0;E_m0;I_m0;R_m0;S_o0;E_o0;I_o0;R_o0];

tspan = 0:200*365;

[~,y_H] = ode45(@(t,y) ODE(WAIFW_H,0,y,t),tspan,y0);
[~,y_A] = ode45(@(t,y) ODE(WAIFW_A,0.86,y,t),tspan,y0);
[t,y_B] = ode45(@(t,y) ODE(WAIFW_B,0.86,y,t),tspan,y0);

[~,y_A_vac] = ode45(@(t,y) ODE(WAIFW_A,0.86,y,t),tspan,y0);
[t,y_B_vac] = ode45(@(t,y) ODE(WAIFW_B,0.86,y,t),tspan,y0);


%% visulization

% Number of new infection(Homogeneous)

NI_H_y = f*(y_H(1:end-1,2) + y_H(2:end,2))/2/N_y*100000;
NI_H_m = f*(y_H(1:end-1,6) + y_H(2:end,6))/2/N_m*100000;
NI_H_o = f*(y_H(1:end-1,10) + y_H(2:end,10))/2/N_o*100000;

% Average number of new infection(Homogeneous)

Avg_NI_H_y = sum(NI_H_y(100*365:end))/length(NI_H_y(100*365:end))
Avg_NI_H_m = sum(NI_H_m(100*365:end))/length(NI_H_m(100*365:end))
Avg_NI_H_o = sum(NI_H_o(100*365:end))/length(NI_H_o(100*365:end))

% Number of new infection per 100000(WAIFW A)

NI_A_y = f*(y_A(1:end-1,2) + y_A(2:end,2))/2/N_y*100000;
NI_A_m = f*(y_A(1:end-1,6) + y_A(2:end,6))/2/N_m*100000;
NI_A_o = f*(y_A(1:end-1,10) + y_A(2:end,10))/2/N_o*100000;

% Average number of new infection(WAIFW A)

Avg_NI_A_y = sum(NI_A_y(100*365:end))/length(NI_A_y(100*365:end));
Avg_NI_A_m = sum(NI_A_m(100*365:end))/length(NI_A_m(100*365:end));
Avg_NI_A_o = sum(NI_A_o(100*365:end))/length(NI_A_o(100*365:end));

% Number of new infection per 100000(WAIFW B)

NI_B_y = f*(y_B(1:end-1,2) + y_B(2:end,2))/2/N_y*100000;
NI_B_m = f*(y_B(1:end-1,6) + y_B(2:end,6))/2/N_m*100000;
NI_B_o = f*(y_B(1:end-1,10) + y_B(2:end,10))/2/N_o*100000;

% Average number of new infection(WAIFW B)

Avg_NI_B_y = sum(NI_B_y(100*365:end))/length(NI_B_y(100*365:end));
Avg_NI_B_m = sum(NI_B_m(100*365:end))/length(NI_B_m(100*365:end));
Avg_NI_B_o = sum(NI_B_o(100*365:end))/length(NI_B_o(100*365:end));

%% Visualization

figure(1)
hold on;
title('Proportion of susceptible(WAIFW A)');
plot(t/365,y_A(:,1)/N_y);
plot(t/365,y_A(:,5)/N_m);
plot(t/365,y_A(:,9)/N_o);
legend('young','middle','old')
ylim([0 1])
hold off;

figure(2)
hold on;
title('Daily number of new infections per 100000(WAIFW A)');
% plot(t(1:end-1)/365,NI_A_y);
plot(t(1:end-1)/365,NI_A_m);
plot(t(1:end-1)/365,NI_A_o);
legend('middle','old')
ylim([0 5])
hold off;

figure(3)
hold on;
title('Proportion of susceptible(WAIFW B)');
plot(t/365,y_B(:,1)/N_y);
plot(t/365,y_B(:,5)/N_m);
plot(t/365,y_B(:,9)/N_o);
legend('young','middle','old')
ylim([0 1])
hold off;

figure(4)
hold on;
title('Daily number of new infections per 100000(WAIFW B)');
% plot(t(1:end-1)/365,NI_B_y);
plot(t(1:end-1)/365,NI_B_m);
plot(t(1:end-1)/365,NI_B_o);
legend('middle','old')
ylim([0 5])
hold off;

