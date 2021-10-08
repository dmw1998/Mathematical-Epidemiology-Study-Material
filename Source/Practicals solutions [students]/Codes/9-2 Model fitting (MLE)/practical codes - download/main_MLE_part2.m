clear all; close all;

%% PART 2 - measles incidence

% data
T = readmatrix('incidence_measles.csv');
time_stamp = T(:,1);        % data time
time_stampc = 0:140;   % continuous time stamp for prediction
data = T(:,2);              % positive

param.time_stamp = time_stamp;
param.time_stampc = time_stampc;
param.data = data;

% fixed parameter
N = $$$$$;
initial = [N-1 0 1 0];
f = $$$$$;
r = $$$$$;

param.initial = initial;
param.f = f;
param.r = r;

%% Question 7 
% Fill "solve_SEIR"

% 7-1 MLE procedure
fprintf('start MLE procedure \n')
nloglf = @(theta) -sum(log(max(eps,poisspdf($$$$$))));
theta_MLE = fminsearch(nloglf,$$$$$);

% 7-2 compute CI using Chi-square dist.
theta_MLE_ci = zeros(1,2);
chi2 = @(x) chi2cdf(x,1) - 0.95;
chival_95 = fzero(chi2,2); 
best_val = nloglf($$$$$);

nln = @(theta) nloglf(theta)-(chival_95/2+best_val);
theta_MLE_ci(1) = fzero(nln,$$$$$);
theta_MLE_ci(2) = 1.857517571551e-5;        % cannot find using fzero

% 7-3 best value and its CI
fprintf('Q7. best-fitting value : %e \n',$$$$$);
fprintf('Q7. negative log-likelihood value of %f : %f \n',$$$$$,$$$$$);
fprintf('Q7. CI for force of infection : [%e %e] \n\n',...
    $$$$$,$$$$$);

% 7-4 best-fitting plot
y_theta_MLE = solve_SEIR(theta_MLE,param);

figure1 = figure('pos',[10 10 600 400]);
plot(time_stamp,data,'.','MarkerSize',20);  % data plot
hold on; 
plot(time_stamp,y_theta_MLE,'LineWidth',2)    % prediction plot
temp_legend = sprintf('Prediction(\\beta = %f)',theta_MLE);
legend('Data',temp_legend)
xlabel('Time (day)')
ylabel('Weekly incidence of measles')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q7-1 fitting', 'epsc')

%% Question 8
% 8-1 R0, H
R0 = $$$$$;
H = $$$$$;

R0_ci = $$$$$;
H_ci = $$$$$;

fprintf('Q8. basic reproduction number : %f \n',R0)
fprintf('Q8. herd immunity threshold   : %f \n',H)
fprintf('Q8. CI for basic reproduction number : [%f %f] \n',R0_ci(1),R0_ci(2))
fprintf('Q8. CI for herd immunity threshold   : [%f %f] \n\n',H_ci(1),H_ci(2))










