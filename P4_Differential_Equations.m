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

% Set up the SEIR model of the transmission dynamics of measles in a closed
% population using differential equations:

trans_period = 150;				% model of 150 days

odefun = @(t,y) [-beta*y(1)*y(3);...
    beta*y(1)*y(3) - kappa*y(2);...
    kappa*y(2) - alpha*y(3);...
    alpha*y(3)];

tspan = [0, trans_period];               % Dutring the first 150 days

y0 = [S0, E0, I0, R0];

[t,y] = ode45(odefun,tspan,y0);

% 1. Plot a graph for number of susceptible, pre-Infectious, infectious,
% and recovered populations during 150 days.
figure('WindowState','maximized')
subplot(3,1,2)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 8 days pre-infectious period');
xlabel('Time(Days)');
ylabel('Population');

% 2. How does the graph change if you change the pre-infectious period to
% be 5 days and 20 days, respectively?

kappa = 1/5; 					% Pre-infectious period 5 days

odefun = @(t,y) [-beta*y(1)*y(3);...
    beta*y(1)*y(3) - kappa*y(2);...
    kappa*y(2) - alpha*y(3);...
    alpha*y(3)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,1)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 5 days pre-infectious period');
xlabel('Time(Days)');
ylabel('Population');

kappa = 1/20; 					% Pre-infectious period 20 days

odefun = @(t,y) [-beta*y(1)*y(3);...
    beta*y(1)*y(3) - kappa*y(2);...
    kappa*y(2) - alpha*y(3);...
    alpha*y(3)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,3)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 20 days pre-infectious period');
xlabel('Time(Days)');
ylabel('Population');

% 3. What happens to the size of the epidemic (as reflected in the number
% of people who are immune at the end) as the pre-infectious period is
% increased? What happens as it is decreased?
% As the pre-infectious period increasing, the size of the epidemic do not
% change but the period of it becomes longer.

%%  Incorporating births and deaths
% Modify the model to include the births and deaths in each time step and
% change the parameter values back to original:
LE = 70 * 365;
mu = 1/LE;          % Birth rate = death rate

kappa = 1/8;

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

% 1. Plot a graph for number of susceptible, pre-Infectious, infectious,
% and recovered populations during 150 days.
figure('WindowState','maximized')
subplot(3,1,2)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 8 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

% 2. How does the graph change if you change the pre-infectious period to
% be 5 days and 20 days, respectively?

kappa = 1/5; 					% Pre-infectious period 5 days

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,1)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 5 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

kappa = 1/20; 					% Pre-infectious period 20 days

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,3)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 150 days with 20 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

% 4. How do the general patterns in the numbers of susceptible and immune
% individuals differ from those predicted using the difference equations in
% the last practical?

% There is no big difference between the two results. One possible reason
% is that the life expectancy is too long comparing with the observed time.

% Hence we change to a longer transmit period, like 50 years.

trans_period = 50 * 365;				% model of 50 years
tspan = [0, trans_period];

kappa = 1/8;

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

figure('WindowState','maximized')
subplot(3,1,2)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 50 years with 8 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

kappa = 1/5; 					% Pre-infectious period 5 days

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,1)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 50 years with 5 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

kappa = 1/20; 					% Pre-infectious period 20 days

odefun = @(t,y) [-beta*y(1)*y(3) + mu*(N-y(1));...
    beta*y(1)*y(3) - kappa*y(2) - mu*y(2);...
    kappa*y(2) - alpha*y(3) - mu*y(3);...
    alpha*y(3) - mu*y(4)];

[t,y] = ode45(odefun,tspan,y0);

subplot(3,1,3)
plot(t, y);
legend('susceptible','pre-infectious','infectious','recover','Location', 'east');
title('SEIR model for 50 years with 20 days pre-infectious period with birth and death');
xlabel('Time(Days)');
ylabel('Population');

% Then we can find that the longer pre-infectious period, the longer
% cycle of the breaking out.

% 5. You should notice that, although the daily number of new infectious
% persons and the numbers of susceptible and immune individuals oscillates
% over time, these oscillations become weaker, and they seem to disappear
% entirely. This pattern is inconsistent with what happens in many
% populations in which measles vaccination has not been introduced, measles
% epidemics occur every two years, which suggests that other factors help
% to sustain the epidemic cycles. Suggest possible factors which might help
% to determine the regular patterns in measles cycles.

% The more new born get vaccinated, the less susceptible increase. This
% will help to control the eipdemic.