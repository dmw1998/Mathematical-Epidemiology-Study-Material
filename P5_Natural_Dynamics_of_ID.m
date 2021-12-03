close all; clear all; clc
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

tspan = 1:trans_period;               % Dutring the first 50 years

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

y0 = [S0, E0, I0, R0];

[t,y] = ode45(odefun,tspan,y0);

% Plot a graph for number of susceptible, pre-Infectious, infectious,
% and recovered populations from the 40th year to the 50th year.
figure('WindowState','maximized')
subplot(2,2,1)
R_n = R_0*y(:,1)/N;         % net reproduction number
new_inf = alpha * y(:,2);   % new infectious on daily base f*E(t)
yyaxis left
plot(t(40*365:50*365), R_n(40*365:50*365));
hold on
plot(t(40*365:50*365), ones(10*365+1,1))
yyaxis right
plot(t(40*365:50*365), new_inf(40*365:50*365));
legend('Net Reproduction number','R_n = 1','new infectious','Location','best');
title('SEIR model for 40^{th} to 50^{th} years for net reproduction number');
xlabel('Time(Days)');

subplot(2,2,2)
S_prop_I = y(:,1)./y(:,3);       % proportion of the population is susceptible to infection
new_inf = alpha * y(:,2);
yyaxis left
plot(t(40*365:50*365), S_prop_I(40*365:50*365));
yyaxis right
plot(t(40*365:50*365), new_inf(40*365:50*365));
legend('S : I','new infectious','Location','best');
title('SEIR model for 40^{th} to 50^{th} years for S:I');
xlabel('Time(Days)');

% Herd immunity threshold = 1 - 1/R_0

subplot(2,2,3)
prop_R = y(:,4)/N;
new_inf = beta * y(:,1).* y(:,3);
yyaxis left
plot(t(40*365:50*365), prop_R(40*365:50*365));
hold on
HIT = 1 - 1/R_0;
plot(t(40*365:50*365), HIT * ones(10*365+1,1))
yyaxis right
plot(t(40*365:50*365), new_inf(40*365:50*365));
legend('proportion of immune','Herd immunity threshold','new infectious','Location','best');
title('SEIR model for 40^{th} to 50^{th} years for immune');
xlabel('Time(Days)');

subplot(2,2,4)
prop_S = y(:,1)/N;
prop_R = y(:,4)/N;
plot(t(40*365:50*365), prop_S(40*365:50*365), t(40*365:50*365), prop_R(40*365:50*365),t(40*365:50*365), HIT * ones(10*365+1,1));
legend('proportion of susceptible','proportion of immune','Herd immunity threshold','Location','best');
title('SEIR model for 40^{th} to 50^{th} years for new infectious');
xlabel('Time(Days)');


% Modify the model to include the vaccination which is introduced 50 years
% after the infection has been circulating in the population so that a
% proportion of newborn individuals are effectively vaccinated. Run the
% model for 80 years and plot the proportion of immune and the number of
% new infectious on each side of y-axis.

% 7.	What happens to the number of new infectious persons per day if the
% proportion of the population which is effectively vaccinated is below
% (60%, 90% coverage) the herd immunity threshold? What happens to the
% number of new infectious persons per day if this proportion is above (93%
% coverage) the herd immunity threshold?

c1 = 0.6;       % vaccine 60% coverage
c2 = 0.9;       % vaccine 60% coverage

odefun60 = @(t,y) [-beta*y(1)*y(3) + (1-c1)*mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - c1*mu*y(4)];

odefun90 = @(t,y) [-beta*y(1)*y(3) + (1-c2)*mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - c2*mu*y(4)];

y0_vac = y(trans_period,:);

tspan_vac = 50*365+1:80*365;

[t_novac,y_novac] = ode45(odefun,tspan_vac,y0_vac);
[t_vac60,y_vac60] = ode45(odefun60,tspan_vac,y0_vac);
[t_vac90,y_vac90] = ode45(odefun90,tspan_vac,y0_vac);

figure('WindowState','maximized')
subplot(1,2,1)
tt = [t; t_vac60];
yynovac = [y; y_novac];
yy60 = [y; y_vac60];
plot(tt./365,yynovac(:,3),tt./365,yy60(:,3))
xlim([30 80])
legend('without vaccine','vacine 60% coverage','Location','best');
title('SEIR model wihe 60% coverage');
xlabel('Time(Years)');

subplot(1,2,2)
yy90 = [y; y_vac90];
plot(tt./365,yynovac(:,3),tt./365,yy90(:,3))
xlim([30 80])
legend('without vaccine','vacine 90% coverage','Location','best');
title('SEIR model wihe 90% coverage');
xlabel('Time(Years)');