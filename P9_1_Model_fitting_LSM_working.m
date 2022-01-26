close all; clear all; clc

%% PART I
% Fitting the catalytic model to seroprevalence data to estimate the force
% of infection

% We will first analyze the data from the UK ("seroprevalence_uk"):
% Ignoring the maternal antibodies, write the formula for the proportion of
% age "a" who have ever been infected in terms of force of infection
% "λ".

% Read data
T_uk = readmatrix('seroprevalence_uk.csv');
time_stamp_uk = T_uk(:,1);        % data age
time_stampc_uk = 0:0.01:44;          % continuous time stamp for prediction
Pa_uk = T_uk(:,2);                % positive
Na_uk = T_uk(:,3);                % N
data_uk = Pa_uk./Na_uk;                 % seroprevalence

% 1. Assume the initial value for the force infection in the UK to be 0.12.
% Do you think the true value for the force infection in the UK was greater
% or smaller than that currently assumed? What is the current value for
% squared error?

% LSM procedure
theta_0 = 0.12;

z_uk = 1 - exp(-theta_0*time_stamp_uk);      % z_uk(lambda,a) = 1 - exp(-lambda*a)

LSM_error = sum((z_uk - data_uk).^2);

figure(1)
hold on
plot(time_stamp_uk, data_uk, '.')
plot(time_stamp_uk, z_uk)
ylim([0, 1.05])
title('UK Data : \lambda = 0.12 Estimation')
text(20,0.4,['Least square Error : ',num2str(LSM_error)])
grid on;
hold off

% 2. What is the best-fitting value for the force of infection and the
% current value for squared error? Plot a graph of model predictions and
% observed data.

LSM_func_uk = @(theta) sum(2.*(1-exp(-theta*time_stamp_uk)-data_uk).*exp(-theta*time_stamp_uk).*time_stamp_uk);
theta_LSM_uk = fzero(LSM_func_uk, theta_0);

y_uk = 1 - exp(-theta_LSM_uk*time_stamp_uk);      % z(lambda,a) = 1 - exp(-lambda*a)
LSM_error = sum((y_uk - data_uk).^2);

figure(2)
hold on
plot(time_stamp_uk, data_uk, '.')
plot(time_stamp_uk, y_uk)
ylim([0, 1.05])
title(['UK Data : \lambda = ',num2str(theta_LSM_uk),' Estimation'])
text(20,0.4,['Least square Error : ',num2str(LSM_error)])
grid on;
hold off

% 3. For which age groups does the model underestimate the proportion of
% individuals who are seropositive? For which age groups does it
% overestimate it?

% 15


% 4. According to the formula, what is the average age at infection in the
% UK assuming that the force of infection is independent of age?

A_uk = 1/theta_LSM_uk;        % Average age at infection


% 5. Assuming that the average life expectancy (L) is 60 years, what is the
% R0 for this population according to the expression R0=L/A? What is the
% herd immunity threshold?

L = 60;
R0_uk = L/A_uk;               % R0
H_uk = 1 - 1/R0_uk;           % herd immunity thershold

% Fit the catalytic model to China data ("seroprevalence_china") to
% estimate the force of infection:

% Read data
T_ch = readmatrix('seroprevalence_china.csv');
time_stamp_ch = T_ch(:,1);     % data age
time_stampc_ch = 0:0.01:36;       % continuous time stamp for prediction
Pa_ch = T_ch(:,2);                % positive
Na_ch = T_ch(:,3);                % N
data_ch = Pa_ch./Na_ch;                 % seroprevalence

% 6. Determine the best-fitting force of infection and plot a graph of
% model predictions and observed data. Calculate the average age at
% infection, the R0 (assuming that the life expectancy is the same as that
% in the UK) and herd immunity.

% z = 1 - exp(-theta_0*time_stamp_ch);      % z(lambda,a) = 1 - exp(-lambda*a)

% LSM_error = sum((z - data).^2);

% figure
% hold on
% plot(time_stamp_ch, data, '.')
% plot(time_stamp_ch, z)
% ylim([0, 1.05])
% title('China Data : \lambda = 0.12 Estimation')
% text(20,0.4,['Least square Error : ',num2str(LSM_error)])
% grid on;
% hold off

LSM_func_ch = @(theta) sum(2.*(1-exp(-theta*time_stamp_ch)-data_ch).*exp(-theta*time_stamp_ch).*time_stamp_ch);
theta_LSM_ch = fzero(LSM_func_ch, theta_0);

y_ch = 1 - exp(-theta_LSM_ch*time_stamp_ch);      % z(lambda,a) = 1 - exp(-lambda*a)
LSM_error = sum((y_ch - data_ch).^2);

figure(3)
hold on
plot(time_stamp_ch, data_ch, '.')
plot(time_stamp_ch, y_ch)
ylim([0, 1.05])
title(['UK Data : \lambda = ',num2str(theta_LSM_ch),' Estimation'])
text(20,0.4,['Least square Error : ',num2str(LSM_error)])
grid on;
hold off

A_ch = 1/theta_LSM_ch;        % Average age at infection
R0_ch = L/A_ch;               % R0
H_ch = 1 - 1/R0_ch;           % herd immunity thershold

% 7. How do the values for the force of infection, average age at
% infection, R0 and herd immunity threshold in China compare against those
% for the UK? Suggest possible reasons for these differences.

fprintf('Average age at infection \n')
fprintf('UK: %f \n', A_uk)
fprintf('China: %f \n', A_ch)

fprintf('\n R_0 \n')
fprintf('UK: %f \n', R0_uk)
fprintf('China: %f \n', R0_ch)

fprintf('\n Herd immunity \n')
fprintf('UK: %f \n', H_uk)
fprintf('China: %f \n', H_ch)

% Since R_0 of China is larger than of the UK, the average age at infection
% of China is much smaller than of the UK. For the same reason, herd
% immunity of China is also larger than of the UK.

% 8. Modify the expression for the prevalence of previous infection at each
% age assuming that individuals are immune for the first 6 months of life
% as a result of maternal antibodies.

% 6 months = 0.5 year

% 9. Assuming that individuals are immune for the first 6 months of life
% and are then susceptible, refit the model to estimate the force of
% infection in the UK and China. Plot a graph of model predictions and
% observed data.

LSM_func_uk2 = @(theta) sum(2.*(1-exp(-theta*(time_stamp_uk-0.5))-data_uk).*exp(-theta*(time_stamp_uk-0.5)).*(time_stamp_uk-0.5));
theta_LSM_uk2 = fzero(LSM_func_uk2, theta_0);

y_uk2 = 1 - exp(-theta_LSM_uk2*(time_stamp_uk-0.5));      % z(lambda,a) = 1 - exp(-lambda*a')
LSM_error2 = sum((y_uk2 - data_uk).^2);

figure(4)
hold on
plot(time_stamp_uk, data_uk, '.')
plot(time_stamp_uk, y_uk)
plot(time_stamp_uk, y_uk2)
ylim([0, 1.05])
title(['UK Data : \lambda = ',num2str(theta_LSM_uk),', ',num2str(theta_LSM_uk2),' Estimation'])
legend('Given data', 'Original Est.', 'With Maternal Imm.','Location','best')
text(18,0.4,['Least square Error : ',num2str(LSM_error),' \lambda = ',num2str(theta_LSM_uk)])
text(18,0.3,['Least square Error : ',num2str(LSM_error2),' \lambda = ',num2str(theta_LSM_uk2)])
grid on;
hold off

LSM_func_ch2 = @(theta) sum(2.*(1-exp(-theta*(time_stamp_ch-0.5))-data_ch).*exp(-theta*(time_stamp_ch-0.5)).*(time_stamp_ch-0.5));
theta_LSM_ch2 = fzero(LSM_func_ch2, theta_0);

y_ch2 = 1 - exp(-theta_LSM_ch2*(time_stamp_ch-0.5));      % z(lambda,a) = 1 - exp(-lambda*a')
LSM_error = sum((y_ch2 - data_ch).^2);

figure(5)
hold on
plot(time_stamp_ch, data_ch, '.')
plot(time_stamp_ch, y_ch)
plot(time_stamp_ch, y_ch2)
ylim([0, 1.05])
title(['UK Data : \lambda = ',num2str(theta_LSM_ch),', ',num2str(theta_LSM_ch2),' Estimation'])
legend('Given data', 'Original Est.', 'With Maternal Imm.','Location','best')
text(17,0.4,['Least square Error : ',num2str(LSM_error),' \lambda = ',num2str(theta_LSM_ch)])
text(17,0.3,['Least square Error : ',num2str(LSM_error2),' \lambda = ',num2str(theta_LSM_ch2)])
grid on;
hold off

% Plot the graphs of -ln(Sa/Na) for China and the UK, where Sa is the
% number of susceptible at age "a" and Na is the number of population at
% age "a".

Sa_uk = Na_uk-Pa_uk;
Sa_ch = Na_ch-Pa_ch;

ln_uk = -log(Sa_uk./Na_uk);
ln_ch = -log(Sa_ch./Na_ch);

figure(6)
hold on
plot(time_stamp_uk, ln_uk, '.', 'MarkerSize', 10)
title('-ln(Sa/Na) for the UK')
grid on
hold off

figure(7)
hold on
plot(time_stamp_ch, ln_ch, '.', 'MarkerSize', 10)
title('-ln(Sa/Na) for China')
grid on
hold off

% 10. According to the plots, is the assumption that the force of infection
% is independent of age in these populations justified? At what age does it
% look as though the force of infection changes in these populations?

% It may be difficult to determine whether the FOI is depend on age only
% from those plots. 

% From the plots in Q3, we can know that the FOI may depend on the age. We
% can around the age as 15 years old. Then we can define a piecewise
% funciton to simulate the data.

% 11. Estimate the age-specific forces of infection using 2 age groups for
% the UK and China. Suggest reasons for the differences in the force of
% infection between China and the UK.

% Divide the observe data sets by age (less than 15 and larger than 15)
% Observing the data sets

f1_uk = @(theta) sum(2.*(1-exp(-theta*time_stamp_uk(1:15))-data_uk(1:15)).*exp(-theta*time_stamp_uk(1:15)).*time_stamp_uk(1:15));
theta_LSM_pf_uk_1 = fzero(f1_uk,theta_0);

f2_uk = @(theta) sum(2.*(1-exp(-15*(theta_LSM_pf_uk_1-theta)).*exp(-theta*time_stamp_uk(16:end))-data_uk(16:end)).*((time_stamp_uk(16:end)-15).*exp((-15*theta_LSM_pf_uk_1+theta*15-theta*time_stamp_uk(16:end)))));
theta_LSM_pf_uk_2 = fzero(f2_uk,theta_0);

y_pf_uk = time_stampc_uk;
y_pf_uk(1:1500) = 1 - exp(-theta_LSM_pf_uk_1*time_stampc_uk(1:1500));      % z(lambda1,a) = 1 - exp(-lambda1*a)
y_pf_uk(1501:end) = 1 - exp(-15*(theta_LSM_pf_uk_1-theta_LSM_pf_uk_2))*exp(-theta_LSM_pf_uk_2*time_stampc_uk(1501:end));
% z(lambda2) = 1 - exp(-15(lambda1-lambda2)*exp(-lambda2*a)

LSM_error = sum((y_pf_uk - data_uk).^2);

figure(8)
hold on
plot(time_stamp_uk, data_uk, '.')
plot(time_stamp_uk, y_uk)
plot(time_stampc_uk, y_pf_uk)
ylim([0, 1.05])
title(['UK Data : \lambda = ',num2str(theta_LSM_uk),' Estimation'])
legend('Given data', 'Original Est.', 'Age-dependent Est.','Location','best')
text(20,0.4,['Least square Error : ',num2str(LSM_error)])
grid on;
hold off


f1_ch = @(theta) sum(2.*(1-exp(-theta*time_stamp_ch(1:7))-data_ch(1:7)).*exp(-theta*time_stamp_ch(1:7)).*time_stamp_ch(1:7));
theta_LSM_pf_ch_1 = fzero(f1_ch,theta_0);

f2_ch = @(theta) sum(2.*(1-exp(-15*(theta_LSM_pf_ch_1-theta)).*exp(-theta*time_stamp_ch(8:end))-data_ch(8:end)).*((time_stamp_ch(8:end)-15).*exp((-15*theta_LSM_pf_ch_1+theta*15-theta*time_stamp_ch(8:end)))));
theta_LSM_pf_ch_2 = fzero(f2_ch,theta_0);

y_pf_ch = time_stampc_ch;
y_pf_ch(1:1500) = 1 - exp(-theta_LSM_pf_ch_1*time_stampc_ch(1:1500));      % z(lambda1,a) = 1 - exp(-lambda1*a)
y_pf_ch(1501:end) = 1 - exp(-15*(theta_LSM_pf_ch_1-theta_LSM_pf_ch_2))*exp(-theta_LSM_pf_ch_2*time_stampc_ch(1501:end));
% z(lambda2) = 1 - exp(-15(lambda1-lambda2)*exp(-lambda2*a)

LSM_error = sum((y_pf_ch - data_ch).^2);

figure(9)
hold on
plot(time_stamp_ch, data_ch, '.')
plot(time_stamp_ch, y_ch)
plot(time_stampc_ch, y_pf_ch)
ylim([0, 1.05])
title(['China Data : \lambda = ',num2str(theta_LSM_ch),' Estimation'])
legend('Given data', 'Original Est.', 'Age-dependent Est.','Location','best')
text(20,0.4,['Least square Error : ',num2str(LSM_error)])
grid on;
hold off