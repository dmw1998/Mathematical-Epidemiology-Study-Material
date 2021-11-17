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
data_coll_A = data_coll;
for t = 2 : trans_period
    II = [data_coll_A{1}(t-1,3); data_coll_A{2}(t-1,3); data_coll_A{3}(t-1,3)];
    lambda = WAIFW_A  * II;
    diff_mat_y = [1-lambda(1)-a_y, 0, 0, 0; ...
        lambda(1), 1-f-a_y, 0, 0;...
        0, f, 1-r-a_y, 0; ...
        0, 0, r, 1-a_y];
    data_coll_A{1}(t,:) = (diff_mat_y*data_coll_A{1}(t-1,:)' )' + [N0*b, 0, 0, 0];

    diff_mat_m = [1-lambda(2)-a_m, 0, 0, 0; ...
        lambda(2), 1-f-a_m, 0, 0;...
        0, f, 1-r-a_m, 0; ...
        0, 0, r, 1-a_m];
    data_coll_A{2}(t,:) = (diff_mat_m*data_coll_A{2}(t-1,:)')' + a_y * data_coll_A{1}(t-1,:);

    diff_mat_o = [1-lambda(3)-d, 0, 0, 0; ...
        lambda(3), 1-f-d, 0, 0;...
        0, f, 1-r-d, 0; ...
        0, 0, r, 1-d];
    data_coll_A{3}(t,:) = (diff_mat_o*data_coll_A{3}(t-1,:)')' + a_m * data_coll_A{2}(t-1,:);
end

% Picture the proportion of the suspecible
tt = 1:trans_period;
figure('WindowState','maximized')
subplot(1,2,1)
plot(tt/365,data_coll_A{1}(:,1)/N_y, tt/365,data_coll_A{2}(:,1)/N_m, tt/365,data_coll_A{3}(:,1)/N_o)
legend('young','middle','old')
ylim([0 1])

WAIFW_B = [BetaB(1), BetaB(2), BetaB(3); BetaB(2), BetaB(2), BetaB(3); BetaB(3), BetaB(3), BetaB(3)];
data_coll_B = data_coll;
for t = 2 : trans_period
    II = [data_coll_B{1}(t-1,3); data_coll_B{2}(t-1,3); data_coll_B{3}(t-1,3)];
    lambda = WAIFW_B * II;
    diff_mat_y = [1-lambda(1)-a_y, 0, 0, 0; ...
        lambda(1), 1-f-a_y, 0, 0;...
        0, f, 1-r-a_y, 0; ...
        0, 0, r, 1-a_y];
    data_coll_B{1}(t,:) = (diff_mat_y*data_coll_B{1}(t-1,:)' )' + [N0*b, 0, 0, 0];

    diff_mat_m = [1-lambda(2)-a_m, 0, 0, 0; ...
        lambda(2), 1-f-a_m, 0, 0;...
        0, f, 1-r-a_m, 0; ...
        0, 0, r, 1-a_m];
    data_coll_B{2}(t,:) = (diff_mat_m*data_coll_B{2}(t-1,:)')' + a_y * data_coll_B{1}(t-1,:);

    diff_mat_o = [1-lambda(3)-d, 0, 0, 0; ...
        lambda(3), 1-f-d, 0, 0;...
        0, f, 1-r-d, 0; ...
        0, 0, r, 1-d];
    data_coll_B{3}(t,:) = (diff_mat_o*data_coll_B{3}(t-1,:)')' + a_m * data_coll_B{2}(t-1,:);
end

subplot(1,2,2)
plot(tt/365,data_coll_B{1}(:,1)/N_y, tt/365,data_coll_B{2}(:,1)/N_m, tt/365,data_coll_B{3}(:,1)/N_o)
legend('young','middle','old')
ylim([0 1])

beta = 1.0530e-05;
WAIFW_H = beta*ones(3,3);

%% Part II
% Given beta for homogeneous

% Vaccination coverage (86%)
% Use the daily number of new infections per 100,000
c = 0.86;                       % level of vaccination coverage 86%
end_time = 200*365;
% data_coll = {[S(1), I(1), I(1), N_y-S(1)-I(1)-I(1)], [S(2), I(2), I(2), N_m-S(2)-I(2)-I(2)], [S(3), I(3), I(3), N_o-S(3)-I(3)-I(3)]};
data_coll = {[5000, 20, 20, N_y-5000-20-20], [3000, 3, 3, N_m-3000-3-3], [2700, 3, 3, N_o-2700-3-3]};

data_coll_A = data_coll;
for t = 2 : trans_period
    II = [data_coll_A{1}(t-1,3); data_coll_A{2}(t-1,3); data_coll_A{3}(t-1,3)];
    lambda = WAIFW_A  * II;
    diff_mat_y = [1-lambda(1)-a_y, 0, 0, 0; ...
        lambda(1), 1-f-a_y, 0, 0;...
        0, f, 1-r-a_y, 0; ...
        0, 0, r, 1-a_y];
    data_coll_A{1}(t,:) = (diff_mat_y*data_coll_A{1}(t-1,:)' )' + [N0*b, 0, 0, 0];

    diff_mat_m = [1-lambda(2)-a_m, 0, 0, 0; ...
        lambda(2), 1-f-a_m, 0, 0;...
        0, f, 1-r-a_m, 0; ...
        0, 0, r, 1-a_m];
    data_coll_A{2}(t,:) = (diff_mat_m*data_coll_A{2}(t-1,:)')' + a_y * data_coll_A{1}(t-1,:);

    diff_mat_o = [1-lambda(3)-d, 0, 0, 0; ...
        lambda(3), 1-f-d, 0, 0;...
        0, f, 1-r-d, 0; ...
        0, 0, r, 1-d];
    data_coll_A{3}(t,:) = (diff_mat_o*data_coll_A{3}(t-1,:)')' + a_m * data_coll_A{2}(t-1,:);
end

vac = [(1-c)*b*N0, 0, 0, c*b*N0];
for t = trans_period+1 : end_time
    II = [data_coll_A{1}(t-1,3); data_coll_A{2}(t-1,3); data_coll_A{3}(t-1,3)];
    lambda = WAIFW_A * II;
    diff_mat_y = [1-lambda(1)-a_y, 0, 0, 0; ...
        lambda(1), 1-f-a_y, 0, 0;...
        0, f, 1-r-a_y, 0; ...
        0, 0, r, 1-a_y];
    data_coll_A{1}(t,:) = data_coll_A{1}(t-1,:) * diff_mat_y + vac;

    diff_mat_m = [1-lambda(2)-a_m, 0, 0, 0; ...
        lambda(2), 1-f-a_m, 0, 0;...
        0, f, 1-r-a_m, 0; ...
        0, 0, r, 1-a_m];
    data_coll_A{2}(t,:) = data_coll_A{2}(t-1,:) * diff_mat_m + a_y * data_coll_A{1}(t-1,:);

    diff_mat_o = [1-lambda(3)-d, 0, 0, 0; ...
        lambda(3), 1-f-d, 0, 0;...
        0, f, 1-r-d, 0; ...
        0, 0, r, 1-d];
    data_coll_A{3}(t,:) = data_coll_A{3}(t-1,:) * diff_mat_o + a_m * data_coll_A{2}(t-1,:);
end

figure
tt = 1:end_time;
plot(tt/365, data_coll_A{2}(:,3)/N_m)
ylim([0 1])