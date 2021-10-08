clc; clear all; close all;

%% PART 1 - Seroprevalence china

% data
T = readmatrix('seroprevalence_china.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:36;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence


% MLE procedure
fprintf('start MLE procedure \n')
nloglf = @(theta) -sum(log(binopdf(Pa,Na,1-exp(-theta*time_stamp))));
theta_MLE = fminsearch(nloglf,0.2);

% compute CI using Chi-square dist.
theta_MLE_ci = zeros(1,2);
chi2 = @(x) chi2cdf(x,1) - 0.95;
chival_95 = fzero(chi2,2); 
best_val = nloglf(theta_MLE);

nln = @(theta) nloglf(theta)-(chival_95/2+best_val);
theta_MLE_ci(1) = fzero(nln,theta_MLE*0.9);
theta_MLE_ci(2) = fzero(nln,theta_MLE*1.1);

%% Question 4 

% 4-1 best-fitting value and log-likelihood value
A = 1/theta_MLE;        % average age of infection
L = 60;                 % life expentancy
R0 = L/A;               % basic reproduction number
H = 1-1/R0;             % herd immunity thershold

fprintf('Q4. best-fitting value : %f \n',theta_MLE);
fprintf('Q4. average age of infection  : %f \n',A)
fprintf('Q4. basic reproduction number : %f \n',R0)
fprintf('Q4. herd immunity threshold   : %f \n',H)

% 4-2 95% CI
A_ci = 1./theta_MLE_ci;
R0_ci = L./A_ci;
H_ci = 1-1./R0_ci;

fprintf('Q4. CI for force of infection : [%f %f] \n',theta_MLE_ci(1),theta_MLE_ci(2));
fprintf('Q4. CI for average age of infection  : [%f %f] \n',A_ci(2),A_ci(1))
fprintf('Q4. CI for basic reproduction number : [%f %f] \n',R0_ci(1),R0_ci(2))
fprintf('Q4. CI for herd immunity threshold   : [%f %f] \n\n',H_ci(1),H_ci(2))


% 4-3 best-fitting plot
y_theta_MLE = solve_catalytic(theta_MLE,time_stampc);

figure1 = figure('pos',[10 10 600 400]);
plot(time_stamp,data,'.','MarkerSize',20);  % data plot
hold on; 
plot(time_stampc,y_theta_MLE,'LineWidth',2)    % prediction plot
temp_legend = sprintf('Prediction(\\lambda = %f)',theta_MLE);
legend('Data',temp_legend,'Location','east')
xlabel('Age (yrs)')
ylabel('Proportion positive')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q4-1 fitting', 'epsc')
