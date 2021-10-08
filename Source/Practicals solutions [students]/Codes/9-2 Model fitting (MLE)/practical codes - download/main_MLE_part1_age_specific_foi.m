clc; clear all; close all;

%% PART 1-3 Age-specific force of infection - UK, China

% uk data
T = readmatrix('seroprevalence_uk.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence

% structured data for uk
uk.Pa = Pa;
uk.Na = Na;
uk.data = data;
uk.time_stamp = time_stamp;
uk.time_stampc = time_stampc;

% china data
T = readmatrix('seroprevalence_china.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:36;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence

% structured data for china
china.Pa = Pa;
china.Na = Na;
china.data = data;
china.time_stamp = time_stamp;
china.time_stampc = time_stampc;

%% before Question 6
% plot the graph : -ln(Sa/Na) 
figure1 = figure('pos',[10 10 1200 400]);

subplot(1,2,1); % uk
plot(uk.time_stamp,$$$$$,'.','MarkerSize',20);  % data plot
legend('proportion of susceptible - uk','Location','southeast')
xlabel('Age (yrs)')
ylabel('Proportion positive')
grid on; grid minor;
set(gca, 'FontSize', 15)

subplot(1,2,2); % china
plot(china.time_stamp,$$$$$,'.','MarkerSize',20);  % data plot
legend('proportion of susceptible - china','Location','southeast')
xlabel('Age (yrs)')
ylabel('Proportion positive')
grid on; grid minor;
set(gca, 'FontSize', 15)

saveas(gca, 'Q6-1 -log(props)', 'epsc')


%% Question 6
% before the MLE procedure, check a "nloglf_age"

%% 6-1 MLE procedure
fprintf('start MLE procedure \n')
uk.theta_MLE = fminsearch(@(theta) nloglf_age($$$$$),[0.1;0.1]);
lb = [0;0];
ub = [1;1];
china.theta_MLE = fminsearchbnd(@(theta) nloglf_age($$$$$),[0.2;0.2],lb,ub);

%% 6-2 compute CI using Chi-square dist.
uk.theta_MLE_ci = zeros(2,2);
china.theta_MLE_ci = zeros(2,2);
chi2 = @(x) chi2cdf(x,2) - 0.95;
chival_95 = fzero(chi2,2); 
uk.best_val = nloglf_age($$$$$,uk);
china.best_val = nloglf_age($$$$$,china);

    % 6-2-1 uk
    uk.nloglf_age1 = @(theta,theta2) nloglf_age($$$$$,uk);
    uk.nloglf_age2 = @(theta1,theta) nloglf_age($$$$$,uk);
    uk.nln1 = @(theta) nloglf_age...
        ([theta;fminsearch(@(theta2) $$$$$,uk.theta_MLE(2))],uk)...
       -(chival_95/2+uk.best_val);
    uk.nln2 = @(theta) nloglf_age...
        ([fminsearch(@(theta1) $$$$$,uk.theta_MLE(1)),theta],uk)...
        -(chival_95/2+uk.best_val);
    uk.theta_MLE_ci(1,1) = fzero(uk.nln1,$$$$$);
    uk.theta_MLE_ci(1,2) = fzero(uk.nln1,$$$$$);
    uk.theta_MLE_ci(2,1) = fzero(uk.nln2,$$$$$);
    uk.theta_MLE_ci(2,2) = fzero(uk.nln2,$$$$$);
    
    % 6-2-2 china
    china.nloglf_age1 = @(theta,theta2) nloglf_age($$$$$,china);
    china.nloglf_age2 = @(theta1,theta) nloglf_age($$$$$,china);
    china.nln1 = @(theta) nloglf_age...
        ([theta;fminsearch(@(theta2) $$$$$,china.theta_MLE(2))],china)...
        -(chival_95/2+china.best_val);
    china.nln2 = @(theta) nloglf_age...
        ([fminsearch(@(theta1) $$$$$,china.theta_MLE(1)),theta],china)...
        -(chival_95/2+china.best_val);
    china.theta_MLE_ci(1,1) = fzero(china.nln1,$$$$$);
    china.theta_MLE_ci(1,2) = fzero(china.nln1,$$$$$);
    china.theta_MLE_ci(2,1) = 0;
    china.theta_MLE_ci(2,2) = fzero(china.nln2,$$$$$);
    
%% 6-3 plot the graph : model predection and observed data

    % uk
    uk.y_theta_MLE = solve_catalytic_age(uk.theta_MLE,uk.time_stampc);
    
    figure1 = figure('pos',[10 10 600 400]);
    plot(uk.time_stamp,uk.data,'.','MarkerSize',20);  % data plot
    hold on;
    plot(uk.time_stampc,uk.y_theta_MLE,'LineWidth',2)    % prediction plot
    temp_legend = sprintf('Prediction(\\lambda = [%f %f])',uk.theta_MLE(1),uk.theta_MLE(2));
    legend('Data',temp_legend,'Location','southeast')
    xlabel('Age (yrs)')
    ylabel('Proportion positive')
    grid on; grid minor;
    set(gca, 'FontSize', 15)
    saveas(gca, 'Q6-2 fitting-uk', 'epsc')
    
    % china
    china.y_theta_MLE = solve_catalytic_age(china.theta_MLE,china.time_stampc);
    
    figure1 = figure('pos',[10 10 600 400]);
    plot(china.time_stamp,china.data,'.','MarkerSize',20);  % data plot
    hold on;
    plot(china.time_stampc,china.y_theta_MLE,'LineWidth',2)    % prediction plot
    temp_legend = sprintf('Prediction(\\lambda = [%f %f])',china.theta_MLE(1),china.theta_MLE(2));
    legend('Data',temp_legend,'Location','southeast')
    xlabel('Age (yrs)')
    ylabel('Proportion positive')
    grid on; grid minor;
    set(gca, 'FontSize', 15)
    saveas(gca, 'Q6-2 fitting-china', 'epsc')

