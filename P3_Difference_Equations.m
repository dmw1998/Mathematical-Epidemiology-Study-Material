clear all; clc
%% Setting up difference equations
% We assume that individuals mix randomly and parameter values are given as
% follows:
N = 100000;						% Population 100,000 people
[S0, E0, I0] = deal(99999, 0, 1);		% Initial values
kappa = 1/8; 					% Pre-infectious period 8 days
alpha = 1/7;					% Infectious period 7 days
R0 = 13;						% Basic reproduction number 13
beta = R0 * alpha / N;

% Set up the SEIR model of the transmission dynamics of measles in a closed
% population using difference equations:

trans_period = 200;				% model of 200 days
dt = 1;							% for each day

data_coll = zeros(4,trans_period+1);
data_coll(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period+1
    diff_mat = [-beta*data_coll(3,i-1), 0, 0, 0; ...
        beta*data_coll(3,i-1), -kappa, 0, 0;...
        0, kappa, -alpha, 0; ...
        0, 0, alpha, 0];
	data_coll(:,i) =data_coll(:,i-1) + dt*(diff_mat*data_coll(:,i-1));
end

% 1. Plot a graph for number of susceptible, pre-Infectious, infectious,
% and recovered populations during 200 days.
days = 0:dt:trans_period;
plot(days, data_coll');
ax = gca;           % current axes
ax.YLim = [-1e4 11e4];
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 200 days');

% 2. How long does it take before there are no infectious persons in the
% population? Why do no further new infectious persons occur in this
% population after a certain time?
% Find the first day t that I(t) is less than 1, which means that the
% infectious is less than one person. We can take it as on infectious
% persons.
Dday = find(data_coll(3,8:trans_period+1) < 1, 1) + 7;      % 137
% The first time that there is no infectious person must later than a week,
% since the infectious period is 7 days. Hence, we search after 7 days and
% add 7 days back for giving the day in trans_period.

% 3. How does the graph change if you change the pre-infectious period to
% be 5 days and 20 days, respectively?
kappa1 = 1/5; 					% Pre-infectious period 5 days
kappa2 = 1/20; 					% Pre-infectious period 20 days

data_coll1 = zeros(4,trans_period+1);
data_coll1(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period+1
    diff_mat = [-beta*data_coll1(3,i-1), 0, 0, 0; ...
        beta*data_coll1(3,i-1), -kappa1, 0, 0;...
        0, kappa1, -alpha, 0; ...
        0, 0, alpha, 0];
	data_coll1(:,i) =data_coll1(:,i-1) + dt*(diff_mat*data_coll1(:,i-1));
end

data_coll2 = zeros(4,trans_period+1);
data_coll2(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period+1
    diff_mat = [-beta*data_coll2(3,i-1), 0, 0, 0; ...
        beta*data_coll2(3,i-1), -kappa2, 0, 0;...
        0, kappa2, -alpha, 0; ...
        0, 0, alpha, 0];
	data_coll2(:,i) =data_coll2(:,i-1) + dt*(diff_mat*data_coll2(:,i-1));
end

figure('WindowState','maximized');
subplot(2,1,1);
plot(days, data_coll1');
ax = gca;           % current axes
ax.YLim = [0 10e4];
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 200 days (Pre-infectious period 5 days)');

subplot(2,1,2); 
plot(days, data_coll2');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 200 days (Pre-infectious period 20 days)');
% As we can see in the figure window, the longer pre-infectious period the
% later infectious person back to zero.

% 4. Which assumptions would you alter or add to the model to describe the
% transmission of measles in a large population over a period of years?
% The natural birth and death should take into account. And we have several
% assuumptions, such as the natural birth rate and the natural death rate
% are the same.

%% Incorporating births and deaths
