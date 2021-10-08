clear; close all; clc;

%% Part 1
%% Basic settings
% Population of compartments
N = 1e5;
I0 = 1;
S0 = N - I0;
E0 = 0;
R0 = 0;

% Parameters
f_val = [1/5, 1/8, 1/20];   % 1/f = pre-infectious period
r = 1/7;                    % 1/r = infectious period
R_0 = 13;                    % basic reproduction number
beta = R_0/N * r;            % transmission rate

% Time
dt = 1;
start_time = 0;
final_time = 365;
time_stamp = start_time:dt:final_time;

%% Solve ODE system
S = zeros(length(time_stamp), length(f_val));
E = zeros(length(time_stamp), length(f_val));
I = zeros(length(time_stamp), length(f_val));
R = zeros(length(time_stamp), length(f_val));
for i = 1:length(f_val)
    f = f_val(i);
    fode = @(t,y) [ - beta .* y(1) .* y(3); ...
        beta .* y(1) .* y(3) - f .* y(2); ...
        f .* y(2) - r .* y(3);
        r .* y(3)];
    
    y0 = [S0; E0; I0; R0];
    [sol_t, sol_y] = ode45(fode, time_stamp, y0); 
    
    S(:,i) = sol_y(:,1);
    E(:,i) = sol_y(:,2);
    I(:,i) = sol_y(:,3);
    R(:,i) = sol_y(:,4);
end

%% Plot
figure1 = figure('pos', [10 10 1200 600]);
subplot(2,2,1)
hold on;
for i = 1:length(f_val)
    plot(time_stamp, S(:,i), 'LineWidth', 2)
end
hold off;
xlabel('time (days)')
ylabel('the number of people')
legend({'f=1/5', 'f=1/8', 'f=1/20'})
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Susceptible (S)')

subplot(2,2,2)
hold on;
for i = 1:length(f_val)
    plot(time_stamp, E(:,i), 'LineWidth', 2)
end
hold off;
xlabel('time (days)')
ylabel('the number of people')
legend({'f=1/5', 'f=1/8', 'f=1/20'})
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Pre-infectious (E)')

subplot(2,2,3)
hold on;
for i = 1:length(f_val)
    plot(time_stamp, I(:,i), 'LineWidth', 2)
end
hold off;
xlabel('time (days)')
ylabel('the number of people')
legend({'f=1/5', 'f=1/8', 'f=1/20'})
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Infectious (I)')

subplot(2,2,4)
hold on;
for i = 1:length(f_val)
    plot(time_stamp, R(:,i), 'LineWidth', 2)
end
hold off;
xlabel('time (days)')
ylabel('the number of people')
legend({'f=1/5', 'f=1/8', 'f=1/20'},'Location','best')
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Recovery (R)')

saveas(gca, 'basic_measles_1y.eps', 'epsc')

%% Part 2
%% Additional Initial Conditions
life_expect = 70*365;       % life expectancy (days)
d = 1/life_expect;          % death rate
b = 1/life_expect;     % birth rate
time_stamp = 0:100*365;
f = 1/8;
%% Compute S(t), E(t), I(t), and R(t) using difference equation
S_de = zeros(length(time_stamp),1);
E_de = zeros(length(time_stamp),1);
I_de = zeros(length(time_stamp),1);
R_de = zeros(length(time_stamp),1);

S_de(1) = S0;
I_de(1) = I0;
for i = 1:length(time_stamp)-1
    S_de(i+1) = S_de(i) - beta .* I_de(i) .* S_de(i) + b .* N - d .* S_de(i);
    E_de(i+1) = E_de(i) + beta .* I_de(i) .* S_de(i) - f .* E_de(i) - d .* E_de(i);
    I_de(i+1) = I_de(i) + f .* E_de(i) - r .* I_de(i) - d .* I_de(i);
    R_de(i+1) = R_de(i) + r .* I_de(i) - d .* R_de(i);
end
%% Compute S(t), E(t), I(t), and R(t) using differential equation
fode = @(t,y) [ b .* N - beta .* y(1) .* y(3) - d .* y(1); ...
    beta .* y(1) .* y(3) - f .* y(2) - d .* y(2); ...
    f .* y(2) - r .* y(3) - d .* y(3); ...
    r .* y(3) - d .* y(4)];

y0 = [S0, E0, I0, R0];
[sol_t, sol_y] = ode45(fode, time_stamp, y0);

S_dle = sol_y(:,1);
E_dle = sol_y(:,2);
I_dle = sol_y(:,3);
R_dle = sol_y(:,4);

%% Plot
figure1 = figure('pos', [10 10 1200 400]);
% Difference equation
subplot(1,2,1)
plot(time_stamp/365, S_de, 'LineWidth', 2)
hold on;
plot(time_stamp/365, R_de, 'LineWidth', 2)
hold off;
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on;
grid minor;
legend('Susceptible', 'Recovery','Location','best')
title('Difference equation')
set(gca, 'FontSize', 12)
% Differential equation
subplot(1,2,2)
plot(time_stamp/365, S_dle, 'LineWidth', 2)
hold on;
plot(time_stamp/365, R_dle, 'LineWidth', 2)
hold off;
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on;
grid minor;
legend('Susceptible', 'Recovery','Location','best')
title('Differential equation')
set(gca, 'FontSize', 12)

saveas(gca, 'demographic_measles_100y_1.eps', 'epsc')

figure2 = figure('pos', [10 10 1200 400]);
% Susceptible
subplot(1,2,1)
plot(time_stamp/365, S_de, 'LineWidth', 2)
hold on;
plot(time_stamp/365, S_dle, 'LineWidth', 2)
hold off;
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on;
grid minor;
legend('Difference', 'Differential','Location','best')
title('Susceptible (S)')
set(gca, 'FontSize', 12)
% Recovery
subplot(1,2,2)
plot(time_stamp/365, R_de, 'LineWidth', 2)
hold on;
plot(time_stamp/365, R_dle, 'LineWidth', 2)
hold off;
xlabel('time (years)')
ylabel('the number of people')
ylim([0 1.1*1e5])
grid on;
grid minor;
legend('Difference', 'Differential','Location','best')
title('Recovery (R)')
set(gca, 'FontSize', 12)

saveas(gca, 'demographic_measles_100y_2.eps', 'epsc')
