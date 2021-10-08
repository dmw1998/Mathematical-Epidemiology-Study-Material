clear; close all; clc;

%% Basic settings
% Population of compartments
N = 1e5;
I0 = 1;
S0 = N - I0;

% Parameters
unit_time_val = [2,3,4,5];
f = 1/8;                    % 1/f = pre-infectious period
r = 1/7;                    % 1/r = infectious period
R0 = 13;                    % basic reproduction number
beta = R0/N * r;
life_expect = 70*365;       % life expectancy (days)
d = 1/life_expect;          % death rate
lambda = 1/life_expect;     % birth rate

%% Compute S(t), E(t), I(t), and R(t)
S = cell(1,length(unit_time_val));
E = cell(1,length(unit_time_val));
I = cell(1,length(unit_time_val));
R = cell(1,length(unit_time_val));

for j = 1:length(unit_time_val)
    % Define time interval
    time_stamp = 1:unit_time_val(j):(50*365+1);
    
    % Memory allocation
    S_temp = zeros(length(time_stamp), 1);
    E_temp = zeros(length(time_stamp), 1);
    I_temp = zeros(length(time_stamp), 1);
    R_temp = zeros(length(time_stamp), 1);
    
    % Initial states
    S_temp(1) = S0;
    I_temp(1) = I0;
    
    % Adjust unit time to parameters
    f_temp = f .* unit_time_val(j);
    r_temp = r .* unit_time_val(j);
    beta_temp = beta .* unit_time_val(j);
    d_temp = d .* unit_time_val(j);
    lambda_temp = lambda .* unit_time_val(j);
    for i = 1:length(time_stamp)-1
        S_temp(i+1) = S_temp(i) - beta_temp .* I_temp(i) .* S_temp(i) + lambda_temp .* N - d_temp .* S_temp(i);
        E_temp(i+1) = E_temp(i) + beta_temp .* I_temp(i) .* S_temp(i) - f_temp .* E_temp(i) - d_temp .* E_temp(i);
        I_temp(i+1) = I_temp(i) + f_temp .* E_temp(i) - r_temp .* I_temp(i) - d_temp .* I_temp(i);
        R_temp(i+1) = R_temp(i) + r_temp .* I_temp(i) - d_temp .* R_temp(i);
    end
    S{1,j} = S_temp;
    E{1,j} = E_temp;
    I{1,j} = I_temp;
    R{1,j} = R_temp;
    
end

%% Plot
figure1 = figure('pos', [10 10 1200 600]);
subplot(2,2,1)
hold on;
for i = 1:length(unit_time_val)
    time_stamp = 1:unit_time_val(i):(50*365+1);
    plot(time_stamp/365, S{1,i}, 'LineWidth', 2)
end
hold off;
xlabel('time (years)')
ylabel('the number of people')
legend({'unit time : 2days', 'unit time : 3days', 'unit time : 4days', 'unit time : 5days'})  
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Susceptible (S)')

subplot(2,2,2)
hold on;
for i = 1:length(unit_time_val)
    time_stamp = 1:unit_time_val(i):(50*365+1);
    plot(time_stamp/365, E{1,i}, 'LineWidth', 2)
end
hold off;
xlabel('time (years)')
ylabel('the number of people')
legend({'unit time : 2days', 'unit time : 3days', 'unit time : 4days', 'unit time : 5days'})  
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Pre-infectious (E)')

subplot(2,2,3)
hold on;
for i = 1:length(unit_time_val)
    time_stamp = 1:unit_time_val(i):(50*365+1);
    plot(time_stamp/365, I{1,i}, 'LineWidth', 2)
end
hold off;
xlabel('time (years)')
ylabel('the number of people')
legend({'unit time : 2days', 'unit time : 3days', 'unit time : 4days', 'unit time : 5days'})  
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Infectious (I)')

subplot(2,2,4)
hold on;
for i = 1:length(unit_time_val)
    time_stamp = 1:unit_time_val(i):(50*365+1);
    plot(time_stamp/365, R{1,i}, 'LineWidth', 2)
end
hold off;
xlabel('time (years)')
ylabel('the number of people')
legend({'unit time : 2days', 'unit time : 3days', 'unit time : 4days', 'unit time : 5days'})  
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Recovery (R)')

saveas(gca, 'demographic_measles_change_unit_time.eps', 'epsc')
