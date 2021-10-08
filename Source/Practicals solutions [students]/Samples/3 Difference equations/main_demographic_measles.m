clear; close all; clc;

%% Basic settings
% Population of compartments
N = 1e5;
I0 = 1;
S0 = N - I0;

% Parameters
f = 1/8;                    % 1/f = pre-infectious period
r = 1/7;                    % 1/r = infectious period
R0 = 13;                    % basic reproduction number
beta = R0/N * r;
life_expect = 70*365;       % life expectancy (days)
d = 1/life_expect;          % death rate
lambda = 1/life_expect;     % birth rate

% Time
time_stamp = 1:(50*365+1);

%% Compute S(t), E(t), I(t), and R(t)
S = zeros(length(time_stamp),1);
E = zeros(length(time_stamp),1);
I = zeros(length(time_stamp),1);
R = zeros(length(time_stamp),1);

S(1) = S0;
I(1) = I0;
for i = 1:length(time_stamp)-1
    S(i+1) = S(i) - beta .* I(i) .* S(i) + lambda .* N - d .* S(i);
    E(i+1) = E(i) + beta .* I(i) .* S(i) - f .* E(i) - d .* E(i);
    I(i+1) = I(i) + f .* E(i) - r .* I(i) - d .* I(i);
    R(i+1) = R(i) + r .* I(i) - d .* R(i);
end


%% Plot
figure1 = figure('pos', [10 10 1200 600]);
subplot(2,2,1)
plot(time_stamp/365, S, 'LineWidth', 2)
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Susceptible (S)')

subplot(2,2,2)
plot(time_stamp/365, E, 'LineWidth', 2)
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Pre-infectious (E)')

subplot(2,2,3)
plot(time_stamp/365, I, 'LineWidth', 2)
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Infectious (I)')

subplot(2,2,4)
plot(time_stamp/365, R, 'LineWidth', 2)
hold off;
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Recovery (R)')

saveas(gca, 'demographic_measles.eps', 'epsc')

