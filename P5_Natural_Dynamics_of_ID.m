clear all; clc
%% Setting up difference equations
% We assume that individuals mix randomly and parameter values are given as
% follows:
N = 100000;						% Population 100,000 people
[S0, E0, I0, R0] = deal(99999, 0, 1, 0);		% Initial values
kappa = 1/8; 					% Pre-infectious period 8 days
alpha = 1/7;					% Infectious period 7 days
R_0 = 13;						% Basic reproduction number 13
beta = R_0 * alpha / N;
LE = 70 * 365;                  % Life expectancy 70 years
mu = 1/LE;                      % Birth rate = death rate

trans_period = 50*365;				% model of 50 years

tspan = [0, trans_period];               % Dutring the first 50 years

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

y0 = [S0, E0, I0, R0];

[t,y] = ode45(odefun,tspan,y0);

% Plot a graph for number of susceptible, pre-Infectious, infectious,
% and recovered populations from the 40th year to the 50th year.
tt = t(1:40*365);
plot(tt./365, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 8 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');