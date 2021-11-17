clear all; close all; clc
%% Setting up varibles (P1)
N0 = 60000;						% Population 60,000 people
[N_y, N_m, N_o] = deal(15000, 15000, 30000); % Population size of young, middle-aged, old group
data_coll = {[15000-1, 0, 1, 0], [15000-1, 0, 1, 0], [30000-1, 0, 1, 0]};
% Initial values
% [S_y, E_y, I_y, R_y], [S_m, E_m, I_m, R_m], [S_o, E_o, I_o, R_o]
% Collect data for each age group. Apply SEIR model for each row of module
% then calculate the trans_rate from younger age to older age. Then apply
% SEIR module for another day.
% Collect the result in each iteration in cooresponding matrix.
f = 1/10; 					    % Pre-infectious period 10 days
r = 1/11;					    % Infectious period 11 days

LE = 60 * 365;                  % Life expectancy 60 years
b = 1/LE;                       % Birth rate
d = 1/(30*365);                 % Death only happend in old group, 30 years

a_y = 1/(15*365);               % rate at which young individuals age
a_m = 1/(15*365);               % rate at which middle-aged individuals age

%% Part I
% Calculate the infectious rate beta for each WAIFW
% By given FOI(force of infection), lambda and Number of infectious, I
lambda = [0.000364; 0.000114; 0.000114];        % P2
S = [5012; 3083; 2739];
I = lambda.*S / r;

% P4
BetaA = lambda./I;              % Beta's for A = diag{beta1, beta2, beta3}; size(3,1)
                                % Beta's for B = [1, 2, 3; 2, 2, 3; 3, 3, 3]
BetaB(3) = lambda(3) / sum(I);
BetaB(2) = (lambda(2) - BetaB(3)*I(3)) / (I(1)+I(2));
BetaB(1) = (lambda(1) - BetaB(2)*I(2) - BetaB(3)*I(3)) / I(1);

% P5 Theoretically WAIFW_B gives more realistic simulation

% P6
trans_period = 100*365;			% model of 100 years
dt = 1;							% for each day

R_0A = 10.9;                    % R_0 for WAIFW A
R_0B = 3.5;                     % R_0 for WAIFW B

WAIFW_A = eye(3).*BetaA;
WAIFW_B = [BetaB(1), BetaB(2), BetaB(3); BetaB(2), BetaB(2), BetaB(3); BetaB(3), BetaB(3), BetaB(3)];

%% Part II
% Given beta for homogeneous

% Vaccination coverage (86%)
% Use the daily number of new infections per 100,000
c = 0.86;                       % level of vaccination coverage 86%
end_time = 200*365;
% y0 = {[S(1), I(1), I(1), N_y-S(1)-I(1)-I(1)], [S(2), I(2), I(2), N_m-S(2)-I(2)-I(2)], [S(3), I(3), I(3), N_o-S(3)-I(3)-I(3)]};
y0 = [5000, 20, 20, N_y-5000-20-20, 3000, 3, 3, N_m-3000-3-3, 2700, 3, 3, N_o-2700-3-3];

vac = @(t) logical(floor(t/(100*365)))*c;

odefunA = @(t,y) [(1-vac(t))*N0*b-WAIFW_A(1,:)*[y(3);y(7);y(11)]*y(1)-a_y*y(1); ...
    WAIFW_A(1,:)*[y(3);y(7);y(11)]*y(1)-f*y(2)-a_y*y(2); ...
    f*y(2)-r*y(3)-a_y*y(3); ...
    vac(t)*b*N0+r*y(3)-a_y*y(4); ...
    a_y*y(1)-WAIFW_A(2,:)*[y(3);y(7);y(11)]*y(5)-a_m*y(5); ...
    WAIFW_A(2,:)*[y(3);y(7);y(11)]*y(5)+a_y*y(2)-f*y(6)-a_m*y(6); ...
    f*y(6)+a_y*y(3)-r*y(7)-a_m*y(7); ...
    r*y(7)+a_y*y(4)-a_m*y(8); ...
    a_m*y(5)-WAIFW_A(3,:)*[y(3);y(7);y(11)]*y(9)-d*y(9); ...
    WAIFW_A(3,:)*[y(3);y(7);y(11)]*y(9)+a_m*y(6)-f*y(10)-d*y(10); ...
    f*y(10)+a_m*y(7)-r*y(11)-d*y(11); ...
    r*y(11)+a_m*y(8)-d*y(12)];

tspan = 0:end_time;

[t,y] = ode45(odefunA,tspan,y0);


figure('WindowState','maximized')
subplot(2,2,1)
plot(t/365, y(:,1)/N_y, t/365, y(:,5)/N_m, t/365, y(:,9)/N_o)
ylim([0 1])
title('Proportion of susceptible(WAIFW A)')
legend('young','middle','old','Location', 'best')

subplot(2,2,2)
plot(t/365, y(:,7)/100000, t/365, y(:,11)/100000)
ylim([0 10E-5])
title('Daily number of new infections / 100,000')
legend('middle','old','Location', 'best')

odefunB = @(t,y) [(1-vac(t))*N0*b-WAIFW_B(1,:)*[y(3);y(7);y(11)]*y(1)-a_y*y(1); ...
    WAIFW_B(1,:)*[y(3);y(7);y(11)]*y(1)-f*y(2)-a_y*y(2); ...
    f*y(2)-r*y(3)-a_y*y(3); ...
    vac(t)*b*N0+r*y(3)-a_y*y(4); ...
    a_y*y(1)-WAIFW_B(2,:)*[y(3);y(7);y(11)]*y(5)-a_m*y(5); ...
    WAIFW_B(2,:)*[y(3);y(7);y(11)]*y(5)+a_y*y(2)-f*y(6)-a_m*y(6); ...
    f*y(6)+a_y*y(3)-r*y(7)-a_m*y(7); ...
    r*y(7)+a_y*y(4)-a_m*y(8); ...
    a_m*y(5)-WAIFW_B(3,:)*[y(3);y(7);y(11)]*y(9)-d*y(9); ...
    WAIFW_B(3,:)*[y(3);y(7);y(11)]*y(9)+a_m*y(6)-f*y(10)-d*y(10); ...
    f*y(10)+a_m*y(7)-r*y(11)-d*y(11); ...
    r*y(11)+a_m*y(8)-d*y(12)];

[t,y] = ode45(odefunB,tspan,y0);

subplot(2,2,3)
plot(t/365, y(:,1)/N_y, t/365, y(:,5)/N_m, t/365, y(:,9)/N_o)
ylim([0 1])
title('Proportion of susceptible(WAIFW A)')
legend('young','middle','old','Location', 'best')

subplot(2,2,4)
plot(t/365, y(:,7)/100000, t/365, y(:,11)/100000)
ylim([0 10E-5])
title('Daily number of new infections / 100,000')
legend('middle','old','Location', 'best')