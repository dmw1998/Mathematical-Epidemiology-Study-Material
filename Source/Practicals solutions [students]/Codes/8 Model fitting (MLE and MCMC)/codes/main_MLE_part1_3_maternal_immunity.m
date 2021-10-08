clc; clear all; close all;

%% PART 1 - Seroprevalence UK and china with maternal immunity

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

%% Question 5
% uk
    % 5-1 MLE procedure
    fprintf('start MLE procedure - uk \n') 
    uk.nloglf = @(theta) -sum(log(max(eps,binopdf...
        (uk.Pa,uk.Na,1-exp(-theta*(uk.time_stamp-0.5))))));
    uk.theta_MLE = fminsearch(uk.nloglf,0.2);
    
    % 5-2 compute CI using Chi-square dist.
    uk.theta_MLE_ci = zeros(1,2);
    chi2 = @(x) chi2cdf(x,1) - 0.95;
    chival_95 = fzero(chi2,2);
    uk.best_val = uk.nloglf(uk.theta_MLE);
    
    uk.nln = @(theta) uk.nloglf(theta)-(chival_95/2+uk.best_val);
    uk.theta_MLE_ci(1) = fzero(uk.nln,uk.theta_MLE*0.8);
    uk.theta_MLE_ci(2) = fzero(uk.nln,uk.theta_MLE*1.2);
    
    % 5-3 values and best-fitting plot
    fprintf('Q5. best-fitting value : %f \n',uk.theta_MLE);
    fprintf('Q5. CI for force of infection : [%f %f] \n\n',...
        uk.theta_MLE_ci(1),uk.theta_MLE_ci(2));
    
    uk.y_theta_MLE = solve_catalytic_MI(uk.theta_MLE,uk.time_stampc);
    
    figure1 = figure('pos',[10 10 600 400]);
    plot(uk.time_stamp,uk.data,'.','MarkerSize',20);  % data plot
    hold on;
    plot(uk.time_stampc,uk.y_theta_MLE,'LineWidth',2)    % prediction plot
    temp_legend = sprintf('Prediction(\\lambda = %f)',uk.theta_MLE);
    legend('Data',temp_legend,'Location','east')
    xlabel('Age (yrs)')
    ylabel('Proportion positive')
    grid on; grid minor;
    set(gca, 'FontSize', 15)
    saveas(gca, 'Q5-1 fitting', 'epsc')

% china
    % 5-4 MLE procedure
    fprintf('start MLE procedure - china \n') 
    china.nloglf = @(theta) -sum(log(max(eps,binopdf...
        (china.Pa,china.Na,1-exp(-theta*(china.time_stamp-0.5))))));
    china.theta_MLE = fminsearch(china.nloglf,0.2);
    
    % 5-5 compute CI using Chi-square dist.
    china.theta_MLE_ci = zeros(1,2);
    chi2 = @(x) chi2cdf(x,1) - 0.95;
    chival_95 = fzero(chi2,2);
    china.best_val = china.nloglf(china.theta_MLE);
    
    china.nln = @(theta) china.nloglf(theta)-(chival_95/2+china.best_val);
    china.theta_MLE_ci(1) = fzero(china.nln,china.theta_MLE*0.8);
    china.theta_MLE_ci(2) = fzero(china.nln,china.theta_MLE*1.2);
    
    % 5-6 values and best-fitting plot
    fprintf('Q5. best-fitting value : %f \n',china.theta_MLE);
    fprintf('Q5. CI for force of infection : [%f %f] \n\n',...
        china.theta_MLE_ci(1),china.theta_MLE_ci(2));
    
    china.y_theta_MLE = solve_catalytic_MI(china.theta_MLE,china.time_stampc);
    
    figure1 = figure('pos',[10 10 600 400]);
    plot(china.time_stamp,china.data,'.','MarkerSize',20);  % data plot
    hold on;
    plot(china.time_stampc,china.y_theta_MLE,'LineWidth',2)    % prediction plot
    temp_legend = sprintf('Prediction(\\lambda = %f)',china.theta_MLE);
    legend('Data',temp_legend,'Location','east')
    xlabel('Age (yrs)')
    ylabel('Proportion positive')
    grid on; grid minor;
    set(gca, 'FontSize', 15)
    saveas(gca, 'Q5-2 fitting', 'epsc')
    
