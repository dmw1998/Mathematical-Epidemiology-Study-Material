clear; close all; clc;

%% Basic settings
% Population of compartments
N = 1e5;
I0 = 1;
S0 = N - I0;

% Parameters
f_val = [1/5, 1/8, 1/20];   % 1/f = pre-infectious period
r = 1/7;                    % 1/r = infectious period
R0 = 13;                    % basic reproduction number
beta = R0/N * r;            % transmission rate

% Time
time_stamp = 1:201;

%% Compute S(t), E(t), I(t), and R(t)
S = zeros(length(time_stamp),length(f_val));
E = zeros(length(time_stamp),length(f_val));
I = zeros(length(time_stamp),length(f_val));
R = zeros(length(time_stamp),length(f_val));
for j = 1:length(f_val)
    f = f_val(j);
    S(1,:) = S0;
    I(1,:) = I0;
    for i = 1:length(time_stamp)-1
        S(i+1,j) = S(i,j) - beta .* I(i,j) .* S(i,j);
        E(i+1,j) = E(i,j) + beta .* I(i,j) .* S(i,j) - f .* E(i,j);
        I(i+1,j) = I(i,j) + f .* E(i,j) - r .* I(i,j);
        R(i+1,j) = R(i,j) + r .* I(i,j);
    end
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
legend({'f=1/5', 'f=1/8', 'f=1/20'})
ylim([0 1.1*1e5])
grid on; grid minor;
set(gca, 'FontSize', 12)
title('Recovery (R)')

saveas(gca, 'basic_measles.eps', 'epsc')
