close all
clear

%% Parameter Setting

N = 1e+5;% Total population size
ave_preinf = 8;% Average duration of pre-infectious(days)
D = 7;% Average duration of infectious(days)
R_0 = 13;% Basic reproduction number
lifeexp = 70;% Life expectancy(years)

% Parameters
beta = R_0 / N / D;
f = 1 / ave_preinf;
r = 1 / D;
d = 1 / lifeexp / 365;
b = d;
H = 1 - 1 / R_0;

% Initial values
S0 = 99999;
E0 = 0;
I0 = N - S0;
R0 = 0;

% Time
dt = 1;% 1 day
t_start = 0;
t_end = 365 * 150;
tspan = t_start:dt:t_end;
nt = length(tspan);

%% Differential Equation

% ODE
fode = @(t,y) [- beta * y(1) * y(3) + b * N - d * y(1);
    beta * y(1) * y(3) - f * y(2) - d * y(2);
    f * y(2) - r * y(3) - d * y(3);
    r * y(3) - d * y(4)];

[~, sol] = ode45(fode, tspan, [S0; E0; I0; R0]);

S = sol(:,1);
E = sol(:,2);
I = sol(:,3);
R = sol(:,4);
prop_sus = S/ N;
prop_imm = R / N;
R_n = R_0 * prop_sus;

% Get new incidence
new_preinf = beta * (S(1:end-1) .* I(1:end-1) + S(2:end) .* I(2:end)) / 2;
new_inf = f * (E(1:end-1) + E(2:end)) / 2;
new_recov = r * (I(1:end-1) + I(2:end)) / 2;

%% Result

figure(1)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlim([tspan(1), tspan(end)])
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, R_n)
xlim([tspan(1), tspan(end)])
xlabel('day')
title('Net reproduction number')
%saveas(gcf, 'result/Q11_1.eps', 'epsc')

figure(2)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlim([3400, 4000])
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, R_n)
xlim([3400, 4000])
xlabel('day')
title('Net reproduction number')
%saveas(gcf, 'result/Q11_2.eps', 'epsc')

figure(3)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, prop_sus)
hold on
plot(tspan, ones(nt,1) / R_0)
hold off
xlabel('day')
title('Proportion of susceptible')
%saveas(gcf, 'result/Q13_1.eps', 'epsc')

figure(4)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlim([3400, 4000])
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, prop_sus)
hold on
plot(tspan, ones(nt,1) / R_0)
hold off
xlim([3400, 4000])
ylim([0, 0.15])
xlabel('day')
title('Proportion of susceptible')

figure(5)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, prop_imm)
hold on
plot(tspan, ones(nt,1) * H)
hold off
xlabel('day')
title('Proportion of immune')
%saveas(gcf, 'result/Q15_1.eps', 'epsc')

figure(6)
subplot(2,1,1)
plot(tspan(2:end), new_inf)
xlim([3400, 4000])
xlabel('day')
title('New infectious')
subplot(2,1,2)
plot(tspan, prop_imm)
hold on
plot(tspan, ones(nt,1) * H)
hold off
xlim([3400, 4000])
ylim([0.8, 1])
xlabel('day')
title('Proportion of immune')
%saveas(gcf, 'result/Q15_2.eps', 'epsc')