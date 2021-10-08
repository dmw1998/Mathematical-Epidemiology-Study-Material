clc; clear all; close all;

%% PART 1 - Seroprevalence UK

% data
T = readmatrix('seroprevalence_uk.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence


% MLE procedure
fprintf('start MLE procedure \n')
nloglf = @(theta) -sum(log(binopdf($$$$$)));
theta_MLE = fminsearch(nloglf,0.1);

% compute CI using Chi-square dist.
theta_MLE_ci = zeros(1,2);
chi2 = @(x) chi2cdf(x,1) - 0.95;
chival_95 = fzero(chi2,2); 
best_val = nloglf($$$$$);

nln = @(theta) nloglf(theta)-(chival_95/2+best_val);
theta_MLE_ci(1) = fzero(nln,$$$$$);
theta_MLE_ci(2) = fzero(nln,$$$$$);

%% Question 1
% 1-1 Fill a blank in "solve_catalytic"

% 1-2 negative log-likelihood value
theta0 = 0.1;
y_theta0 = solve_catalytic(theta0,time_stampc);
fprintf('Q1. negative log-likelihood value of 0.1 : %f \n\n',nloglf(theta0));

% 1-3 plot for assumed value 0.1
figure1 = figure('pos',[10 10 600 400]);
plot(time_stamp,data,'.','MarkerSize',20);  % data plot
hold on; 
plot(time_stampc,y_theta0,'LineWidth',2)    % prediction plot
legend('Data','Prediction(\lambda = 0.1)','Location','east')
xlabel('Age (yrs)')
ylabel('Proportion positive')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q1-1 fitting', 'epsc')

%% Question 2
% 2-1 plot : negative log-likelihood

theta = theta_MLE*0.9 : 0.001 : theta_MLE*1.1;
for i = 1 : length(theta)
    y_nloglf(i) = nloglf(theta(i));
end
figure1 = figure('pos',[10 10 600 400]);
plot(theta,y_nloglf,'LineWidth',2')
title('Negative log likelihood')
xlim([min(theta) max(theta)])
xlabel('Force of Infection \lambda')
ylabel('Value of Negative log likelihood')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q2-1 nloglf', 'epsc')

% 2-2 best-fitting value and log-likelihood value
fprintf('Q2. best-fitting value : %f \n',$$$$$);
fprintf('Q2. negative log-likelihood value of %f : %f \n',$$$$$,$$$$$);

% 2-3 95% CI
fprintf('Q2. CI for force of infection : [%f %f] \n\n',...
    $$$$$,$$$$$);

% 2-4 best-fitting plot
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
saveas(gca, 'Q2-2 fitting', 'epsc')


%% Question 3
% 3-1 calculate A, R0 and H 
A = $$$$$;        % average age of infection
L = $$$$$;                 % life expentancy
R0 = $$$$$;               % basic reproduction number
H = $$$$$;             % herd immunity thershold

fprintf('Q3. average age of infection  : %f \n',A)
fprintf('Q3. basic reproduction number : %f \n',R0)
fprintf('Q3. herd immunity threshold   : %f \n',H)


% 3-2 calculate CI for A, R0 and H
A_ci = $$$$$;
R0_ci = $$$$$;
H_ci = $$$$$;

fprintf('Q3. CI for average age of infection  : [%f %f] \n',A_ci(2),A_ci(1))
fprintf('Q3. CI for basic reproduction number : [%f %f] \n',R0_ci(1),R0_ci(2))
fprintf('Q3. CI for herd immunity threshold   : [%f %f] \n\n',H_ci(1),H_ci(2))





