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
xlabel('Time(Days)');
ylabel('Population');

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
xlabel('Time(Days)');
ylabel('Population');

subplot(2,1,2); 
plot(days, data_coll2');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 200 days (Pre-infectious period 20 days)');
xlabel('Time(Days)');
ylabel('Population');
% As we can see in the figure window, the longer pre-infectious period the
% later infectious person back to zero.

% 4. Which assumptions would you alter or add to the model to describe the
% transmission of measles in a large population over a period of years?
% The natural birth and death should take into account. And we have several
% assuumptions, such as the natural birth rate and the natural death rate
% are the same.

%% Incorporating births and deaths
% 5. Assuming that the population size doesn't change over time, what would
% be an appropriate expression for the per capital birth rate and death
% rate? 
% We can set a constant mu as the rate of birth and death. And all the
% newborns are susceptible. Note that mu is a factor. As our assumption,
% the life expectancy is 70 years. Then we set
LE = 70 * 365;
mu = 1/LE;          % Birth rate = death rate

% 6. Is it realistic that all individuals are born into the susceptible
% population? Is this a reasonable assumption to make? What alternative
% assumptions might be appropriate to make the model more realistic?
% It depends on the properties of diseases. If the immunity is maternal,
% then not all the individuals are born into the susceptible population. If
% it is adaptive immunity, then we can assume this in our model. The birth
% rate and death rate are independent, i.e. they are unnecessarily the
% same. And if the immunity is not adaptive, then we need another factor to
% simulate the portioin of those have natural immunity.


% 7. With the birth and death rates included in the equations, run the
% equations for 200 days. Is there a point where no more infections occur?

data_coll3 = zeros(4,trans_period+1);
data_coll3(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period+1
    diff_mat = [-beta*data_coll3(3,i-1), mu, mu, mu; ...
        beta*data_coll3(3,i-1), -kappa-mu, 0, 0;...
        0, kappa, -alpha-mu, 0; ...
        0, 0, alpha, -mu];
	data_coll3(:,i) =data_coll3(:,i-1) + dt*(diff_mat*data_coll3(:,i-1));
end

figure
days = 0:dt:trans_period;
plot(days, data_coll3');
ax = gca;           % current axes
ax.YLim = [-1e4 11e4];
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 200 days with LE = 70 yrs');
xlabel('Time(Days)');
ylabel('Population');

Dday_200dys = find(data_coll3(3,8:trans_period+1) < 1, 1) + 7;     % 140
% Later than the case without new birth and death. But basicly the same.

% What about 10 years? 50 years?
trans_period1 = 10 * 365;
data_coll4 = zeros(4,trans_period1+1);
data_coll4(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period1+1
    diff_mat = [-beta*data_coll4(3,i-1), mu, mu, mu; ...
        beta*data_coll4(3,i-1), -kappa-mu, 0, 0;...
        0, kappa, -alpha-mu, 0; ...
        0, 0, alpha, -mu];
	data_coll4(:,i) =data_coll4(:,i-1) + dt*(diff_mat*data_coll4(:,i-1));
end

Dday_10yrs = find(data_coll4(3,8:trans_period1+1) < 1, 1) + 7;        % 140

trans_period2 = 50 * 365;
data_coll5 = zeros(4,trans_period2+1);
data_coll5(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period2+1
    diff_mat = [-beta*data_coll5(3,i-1), mu, mu, mu; ...
        beta*data_coll5(3,i-1), -kappa-mu, 0, 0;...
        0, kappa, -alpha-mu, 0; ...
        0, 0, alpha, -mu];
	data_coll5(:,i) =data_coll5(:,i-1) + dt*(diff_mat*data_coll5(:,i-1));
end

Dday_50yrs = find(data_coll5(3,8:trans_period2+1) < 1, 1) + 7;        % 140

figure('WindowState','maximized');
subplot(2,1,1);
days = 0:dt:trans_period1;
plot(days, data_coll4');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 10 years with LE = 70 yrs');
xlabel('Time(Days)');
ylabel('Population');

subplot(2,1,2); 
days = 0:dt:trans_period2;
plot(days, data_coll5');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 50 years with LE = 70 yrs');
xlabel('Time(Days)');
ylabel('Population');

% 8. How are the changes in the number of susceptible and immune related if
% you were to simulate the dynamics of measles for 50 years and for 100
% years?

trans_period3 = 100 * 365;
data_coll6 = zeros(4,trans_period3+1);
data_coll6(:,1) = [S0,E0,I0,R0]';
for i = 2:dt:trans_period3+1
    diff_mat = [-beta*data_coll6(3,i-1), mu, mu, mu; ...
        beta*data_coll6(3,i-1), -kappa-mu, 0, 0;...
        0, kappa, -alpha-mu, 0; ...
        0, 0, alpha, -mu];
	data_coll6(:,i) =data_coll6(:,i-1) + dt*(diff_mat*data_coll6(:,i-1));
end

figure('WindowState','maximized');
subplot(2,1,1);
days = 0:dt:trans_period2;
plot(days, data_coll5');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 50 years with LE = 70 yrs');
xlabel('Time(Days)');
ylabel('Population');

subplot(2,1,2); 
days = 0:dt:trans_period3;
plot(days, data_coll6');
legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
title('SEIR model for 100 years with LE = 70 yrs');
xlabel('Time(Days)');
ylabel('Population');

% As the trans_period increases, we can see that the model reach to the
% equilibrium.

% 9. What happens if you change the time step to 2, 3, 4, and 5 days?

% Compute for 200 days
trans_period = 200;				% model of 200 days
dt_opt = [2, 3, 4, 5];
figure('WindowState','maximized');

for n = 1 : length(dt_opt)
    dt = dt_opt(n);							% for every dt days
    data_coll7 = [];
    data_coll7(:,1) = [S0,E0,I0,R0]';
    i = 1;
    for t = dt+1 : dt : trans_period+1
        i = i + 1;
        diff_mat = [-beta*data_coll7(3,i-1), mu, mu, mu; ...
            beta*data_coll7(3,i-1), -kappa-mu, 0, 0;...
            0, kappa, -alpha-mu, 0; ...
            0, 0, alpha, -mu];
        data_coll7(:,i) =data_coll7(:,i-1) + dt*(diff_mat*data_coll7(:,i-1));
    end

    subplot(2,2,n);
    days = 0:dt:trans_period;
    plot(days, data_coll7');
    legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
    title(sprintf('SEIR model for 200 days (dt = %d)', dt));
    xlabel('Time(Days)');
    ylabel('Population');
end

% Compute for 50 years
trans_period = 50 * 365;				% model of 200 days
figure('WindowState','maximized');

for n = 1 : length(dt_opt)
    dt = dt_opt(n);							% for every dt days
    data_coll7 = [];
    data_coll7(:,1) = [S0,E0,I0,R0]';
    i = 1;
    for t = dt+1 : dt : trans_period+1
        i = i + 1;
        diff_mat = [-beta*data_coll7(3,i-1), mu, mu, mu; ...
            beta*data_coll7(3,i-1), -kappa-mu, 0, 0;...
            0, kappa, -alpha-mu, 0; ...
            0, 0, alpha, -mu];
        data_coll7(:,i) =data_coll7(:,i-1) + dt*(diff_mat*data_coll7(:,i-1));
    end

    subplot(2,2,n);
    days = 0:dt:trans_period;
    plot(days, data_coll7');
    legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
    title(sprintf('SEIR model for 50 years (dt = %d)', dt));
    xlabel('Time(Days)');
    ylabel('Population');
end

% Compute for 200 days and 50 years with dt = 10 days
trans_period = [200, 50 * 365];				% model of 200 days and 50 years
dt = 10;                                                % for every 10 days
figure('WindowState','maximized');

for n = 1 : length(trans_period)
    tsp = trans_period(n);
    data_coll7 = [];
    data_coll7(:,1) = [S0,E0,I0,R0]';
    i = 1;
    for t = dt+1 : dt : tsp+1
        i = i + 1;
        diff_mat = [-beta*data_coll7(3,i-1), mu, mu, mu; ...
            beta*data_coll7(3,i-1), -kappa-mu, 0, 0;...
            0, kappa, -alpha-mu, 0; ...
            0, 0, alpha, -mu];
        data_coll7(:,i) =data_coll7(:,i-1) + dt*(diff_mat*data_coll7(:,i-1));
    end

    subplot(2,1,n);
    days = 0 : dt : tsp;
    plot(days, data_coll7');
    legend('susceptible','pre-infectious','infectious','recover','Location', 'best');
    title(sprintf('SEIR model for 50 years (dt = %d)', dt));
    xlabel('Time(Days)');
    ylabel('Population');
end

 % Would it be reasonable to take a time step of 10 days?
 % No. The information is totally damaged.