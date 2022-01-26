close all; clear all; clc
%% Part I
% Fitting the catalytic model to seroprevalence data to estimate the force
% of infection using MLE

% Ignoring the maternal antibodies, fit the catalytic model to the UK data
% (seroprevalence_uk) to estimate the force of infection using MLE in
% which Binomial distribution is assumed for the data.

% Read data
T = readmatrix('seroprevalence_uk.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence

% 1.	Assume the initial value for the force infection in the UK to be
% 0.1. Do you think the true value for the force infection in the UK was
% greater or smaller than that currently assumed? What is the current value
% for negative log-likelihood?

% MLE procedure
theta_0 = 0.1;

fprintf('start MLE procedure \n')
custnloglf = @(theta) -sum(log(binopdf(Pa,Na,1-exp(-theta*time_stamp))));
theta_MLE = fminsearch(custnloglf,theta_0);         % Wanted best lambda

z = 1 - exp(-theta_0*time_stampc);      % z(lambda,a) = 1 - exp(-lambda*a)

fprintf('Negative log-likelihood value of 0.1 : %f \n\n',custnloglf(theta_0));

% 2.	What is the best-fitting value for the force of infection and the
% current value for negative log-likelihood? Plot a graph of negative
% log-likelihood and estimate the 95% confidence interval.

% compute CI using Chi-square dist.
theta_MLE_ci = zeros(1,2);
chi2 = @(x) chi2cdf(x,1) - 0.95;            % Give the cdf of chi2 with degree of freedom = 1, then -95% 
                                            % to find the zero which is the length of the (vertical) inteval
chival_95 = fzero(chi2,2);                  % around the best value with the given confidence interval
                                            % If the wanted confidence interval changed, then the transformation of function will be changed too.
best_val = custnloglf(theta_MLE);           % likelihood function value evaluated by the minimum = wanted lambda
fprintf('Best negative log-likelihood value of %f : %f \n\n',theta_MLE,best_val);

nln = @(theta) custnloglf(theta)-(chival_95/2+best_val);  % Translate the nloglf to find the confidence interval
                                                          % move chival_95/2+best_val downward (the -2log(lambda) = chi2(k), k=1)
theta_MLE_ci(1) = fzero(nln,theta_MLE*0.9);     % lower bdd of CI
theta_MLE_ci(2) = fzero(nln,theta_MLE*1.1);     % upper bdd of CI

y = 1 - exp(-theta_MLE*time_stampc);      % z(lambda,a) = 1 - exp(-lambda*a)
                                          % Use best value lambda

figure
hold on
scatter(time_stamp, data, '.')
plot(time_stampc,z)
plot(time_stampc,y)
legend('Data','Prediction(\lambda = 0.1)','Prediction(\lambda = best)','Location','best')
xlabel('Age (yrs)')
ylabel('Proportion positive')

% 3.	Calculate the average age at infection, the R0 (assuming that the
% life expectancy is 60 years) and herd immunity threshold with 95%
% confidence interval.
L = 60;
A = 1/theta_MLE;        % Average age at infection
R0 = L/A;               % R0
H = 1 - 1/R0;           % herd immunity thershold

A_ci = 1./theta_MLE_ci;
R0_ci = L./A_ci;
H_ci = 1-1./R0_ci;

% Read data
T = readmatrix('seroprevalence_china.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence

% 4.	Calculate the best-fitting force of infection, the average age at
% infection, the R0 (assuming that the life expectancy is the same as that
% in the UK) and herd immunity threshold with 95% confidence interval.

% MLE procedure
theta_0 = 0.1;

fprintf('start MLE procedure \n')
custnloglf = @(theta) -sum(log(binopdf(Pa,Na,1-exp(-theta*time_stamp))));
theta_MLE = fminsearch(custnloglf,theta_0);         % Wanted best lambda

z = 1 - exp(-theta_0*time_stampc);      % z(lambda,a) = 1 - exp(-lambda*a)

fprintf('Negative log-likelihood value of 0.1 : %f \n\n',custnloglf(theta_0));

%5.     Modify the expression for the prevalence of previous infection at each
%age assuming that individuals are immune for the first 6 months of life
%and are then susceptible. Refit the model to estimate the force of
%infection in the UK and China and plot a graph of model predictions and
%observed data.f negative
% log-likelihood and estimate the 95% confidence interval.

% compute CI using Chi-square dist.
theta_MLE_ci = zeros(1,2);
chi2 = @(x) chi2cdf(x,1) - 0.95;            % Give the cdf of chi2 with degree of freedom = 1, then -95% 
                                            % to find the zero which is the length of the (vertical) inteval
chival_95 = fzero(chi2,2);                  % around the best value with the given confidence interval
                                            % If the wanted confidence interval changed, then the transformation of function will be changed too.
best_val = custnloglf(theta_MLE);           % likelihood function value evaluated by the minimum = wanted lambda
fprintf('Best negative log-likelihood value of 0.1 : %f \n\n',best_val);

nln = @(theta) custnloglf(theta)-(chival_95/2+best_val);  % Translate the nloglf to find the confidence interval
                                                          % move chival_95/2+best_val downward (the -2log(lambda) = chi2(k), k=1)
theta_MLE_ci(1) = fzero(nln,theta_MLE*0.9);     % lower bdd of CI
theta_MLE_ci(2) = fzero(nln,theta_MLE*1.1);     % upper bdd of CI

y = 1 - exp(-theta_MLE*time_stampc);      % z(lambda,a) = 1 - exp(-lambda*a)
                                          % Use best value lambda

figure
hold on
scatter(time_stamp, data, '.')
plot(time_stampc,z)
plot(time_stampc,y)
legend('Data','Prediction(\lambda = 0.1)','Prediction(\lambda = best)','Location','best')
xlabel('Age (yrs)')
ylabel('Proportion positive')

% 6.	Estimate the age-specific forces of infection and 95% confidence
% interval using 2 age groups for the UK and China. Compare the graphs of
% model predictions and observed data using constant force of infection and
% age-specific forces of infection. erd immunity threshold with 95%
% confidence interval.
A = 1/theta_MLE;        % Average age at infection
R0 = L/A;               % R0
H = 1 - 1/R0;           % herd immunity thershold

A_ci = 1./theta_MLE_ci;
R0_ci = L./A_ci;
H_ci = 1-1./R0_ci;