% Set up the SEIR model of the transmission dynamics of measles as follows:
% Population 100000 people
% Pre-infectious period 8 days
% Infectious period 7 days
% Initial values (S,E,I,R)=(99999,0,1,0)

N = 100000;						% Population 100,000 people
[S0, E0, I0, R0] = deal(99999, 0, 1, 0);		% Initial values
kappa = 1/8; 					% Pre-infectious period 8 days
alpha = 1/7;					% Infectious period 7 days

% Fit the SEIR model to prevalence data to estimate the transmission rate
% using MLE in which Poisson distribution is assumed for the data
% ("revalence_measles"):

% Read data
T = readmatrix('incidence_measles.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive

% 7.	What is the best-fitting value for the transmission rate and the
% 95% confidence interval? Plot a graph of model predictions and observed
% data.

% MLE procedure
theta_0 = 0.00001;

fprintf('start MLE procedure \n')
custnloglf = @(theta) -sum(log(poisspdf(Pa,1-exp(-theta*time_stamp))));
theta_MLE = fminsearch(custnloglf,theta_0);         % Wanted best lambda

fprintf('Negative log-likelihood value of %f1 : %f \n\n',theta_0,custnloglf(theta_0));


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

z = 1 - exp(-theta_MLE*time_stampc);      % z(lambda,a) = 1 - exp(-lambda*a)
                                          % Use best value lambda

figure
hold on
scatter(time_stamp, data, '.')
plot(time_stampc,z)
legend('Data','Prediction(\lambda = 0.1)','Prediction(\lambda = best)','Location','best')
xlabel('Age (yrs)')
ylabel('Proportion positive')



% 8.	Calculate R0 and herd immunity threshold with 95% confidence
% interval.

