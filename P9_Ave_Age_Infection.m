clear all; clc
T_uk = readmatrix('seroprevalence_uk.csv');
time_stamp_uk = T_uk(:,1);        % data age
time_stampc_uk = 0:0.01:44;          % continuous time stamp for prediction
Pa_uk = T_uk(:,2);                % positive
Na_uk = T_uk(:,3);                % N
data_uk = Pa_uk./Na_uk;                 % seroprevalence

theta_0 = 0.12;
LSM_func_uk = @(theta) sum(2.*(1-exp(-theta*time_stamp_uk)-data_uk).*exp(-theta*time_stamp_uk).*time_stamp_uk);
theta_LSM_uk = fzero(LSM_func_uk, theta_0);

yc_uk = 1 - exp(-theta_LSM_uk*time_stampc_uk);      % z(lambda,a) = 1 - exp(-lambda*a)

A1 = sum(time_stamp_uk.*(Na_uk-Pa_uk));
A2 = sum((Na_uk-Pa_uk));
A = A1/A2

T_ch = readmatrix('seroprevalence_china.csv');
time_stamp_ch = T_ch(:,1);        % data age
time_stampc_ch = 0:0.01:44;          % continuous time stamp for prediction
Pa_ch = T_ch(:,2);                % positive
Na_ch = T_ch(:,3);                % N
data_ch = Pa_ch./Na_ch;                 % seroprevalence

theta_0 = 0.12;
LSM_func_ch = @(theta) sum(2.*(1-exp(-theta*time_stamp_ch)-data_ch).*exp(-theta*time_stamp_ch).*time_stamp_ch);
theta_LSM_ch = fzero(LSM_func_ch, theta_0);

yc_ch = 1 - exp(-theta_LSM_ch*time_stampc_ch);      % z(lambda,a) = 1 - exp(-lambda*a)

A1 = sum(time_stamp_ch.*(Na_ch-Pa_ch));
A2 = sum((Na_ch-Pa_ch));
A = A1/A2