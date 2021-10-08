close all
clear
restoredefaultpath

%% One-Way Sensitivity Analysis

params = {'N_y', 10000, false
    'prop_immune_y', 0.25, false
    'prop_infection_clinical_y', 0.4, false
    'case_fatality_y', 0.0001, true
    
    'N_o', 40000, false
    'prop_immune_o', 0.4, false
    'prop_infection_clinical_o', 0.4, false
    'case_fatality_o', 0.0002, false
    
    'transmission_coeff', 0.55, true
    'mixing_matrix', [1, 0.15; 0.15, 0.8], false
    
    'gamma', 1 / 3, false
    
    'efficacy', 0.4, false
    'coverage_y.no_vacc', 0, false
    'coverage_o.no_vacc', 0, false
    'coverage_y.A', 0.5, false
    'coverage_o.A', 0, false
    
    'dt', 0.2, false
    't_init', 0, false
    't_termi', 200, false
    
    'cost_per_vacc', 8, false
    'cost_per_clinical_case', 12, true
    'qalyloss_per_clinical_case', 0.005, false
    'qalyloss_per_death', 30, false};

% Set ranges
test_range = [0.75 * 0.0001, 1.25 * 0.0001% case_fatality_y
    0.75 * 0.55, 1.25 * 0.55% transmission_coeff
    0.75 * 12, 1.25 * 12];% cost_per_clinical_case

isTested = cell2mat(params(:, 3));
index4test = find(isTested);
n_parameter = sum(isTested);

% Initialize
ICERs = zeros(n_parameter, 2);

% Loop over parameters
for i = 1:n_parameter
    
    % Loop over test range
    for j = 1:2
        
        % Replace value
        params_test = params;
        params_test{index4test(i), 2} = test_range(i, j);
        
        % CEA
        [~, ~, ~, ~, ICER] = CEA_vacc(params_test);
        
        ICERs(i, j) = ICER;
        
    end
end

% Generate table
array2table([test_range(:, 1), ICERs(:, 1), test_range(:, 2), ICERs(:, 2)], ...
    'VariableNames', {'Lower value', 'ICER lower', 'Upper value', 'ICER upper'}, ...
    'RowNames', {'Case-fatality risk in children', 'Transmission coefficient child-child', 'Cost per clinical case'})

%% Multi-Way Sensitivity Analysis

params = {'N_y', 10000, false
    'prop_immune_y', 0.25, false
    'prop_infection_clinical_y', 0.4, false
    'case_fatality_y', 0.0001, false
    
    'N_o', 40000, false
    'prop_immune_o', 0.4, false
    'prop_infection_clinical_o', 0.4, false
    'case_fatality_o', 0.0002, false
    
    'transmission_coeff', 0.55, false
    'mixing_matrix', [1, 0.15; 0.15, 0.8], false
    
    'gamma', 1 / 3, false
    
    'efficacy', 0.4, false
    'coverage_y.no_vacc', 0, false
    'coverage_o.no_vacc', 0, false
    'coverage_y.A', 0.5, false
    'coverage_o.A', 0, false
    
    'dt', 0.2, false
    't_init', 0, false
    't_termi', 200, false
    
    'cost_per_vacc', 8, true
    'cost_per_clinical_case', 12, true
    'qalyloss_per_clinical_case', 0.005, true
    'qalyloss_per_death', 30, true};

% Set ranges
test_range = [4, 12%cost_per_vacc
    6, 18% cost_per_clinical_case
    0.0025, 0.0075% qalyloss_per_clinical_case
    20, 40];% qalyloss_per_death

n_sample = 100;
isTested = cell2mat(params(:, 3));
n_parameter = sum(isTested);

% Generate samples
test_sample = zeros(n_sample, n_parameter);

rng(0)
for i = 1:n_parameter
    test_min = test_range(i, 1);
    test_max = test_range(i, 2);
    test_sample(:, i) = (test_max - test_min) * rand(1, n_sample) + test_min;
end

% Initialize
CEAs = zeros(n_sample, 7);
RowNames = {};

% Loop over parameters
for i = 1:n_sample
    
    % Replace value
    params(isTested, 2) = num2cell(test_sample(i, :));
    
    % CEA
    [total_cost, total_qalyloss, incremental_cost, incremental_qaly, ICER] = CEA_vacc(params);
    CEAs(i, :) = [total_cost(1), total_qalyloss(1), total_cost(2), total_qalyloss(2), incremental_cost, incremental_qaly, ICER];
    
    RowNames = [RowNames, sprintf('%d', i)];
end

% Calcuate mean and 95% interval
ICER_mean = mean(CEAs(:, 7));
ICER_left95 = quantile(CEAs(:, 7), 0.025);
ICER_right95 = quantile(CEAs(:, 7), 0.975);
fprintf('ICER mean: %f\nICER 95%% interval: (%f, %f)\n', ICER_mean, ICER_left95, ICER_right95)

% Generate table
parameter_names = params(isTested, 1)';
CEA_table = array2table([test_sample, CEAs], 'VariableNames', [parameter_names, ...
    {'No vaccination: total costs', 'No vaccination: total QALY loss', ...
    'Strategy A: total costs', 'Strategy A: total QALY loss', ...
    'Incremental costs', 'Incremental QALYs', 'ICER'}], 'RowNames', RowNames)

% Plot CEA plane
figure(1)
hold on
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5)
h = plot(linspace(14, 34), 5000 * linspace(14, 34), 'k');
hold off
legend(h, 'Threshold')
xlabel('Incremental QALYs')
ylabel('Incremental costs')
title('Cost-effectiveness plane (Uniform distribution)')
saveas(gcf, 'CEA_plane_SA_uniform', 'png')

figure(4)
hold on
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5)

%% Other Choice of Distribution (Normal)

for i = 1:n_parameter
    test_min = test_range(i, 1);
    test_max = test_range(i, 2);
    test_mean = (test_min + test_max) / 2;
    test_sample(:, i) = randn(1, n_sample) * test_mean / 10 + test_mean;
end

% Initialize
CEAs = zeros(n_sample, 7);

% Loop over parameters
for i = 1:n_sample
    
    % Replace value
    params(isTested, 2) = num2cell(test_sample(i, :));
    
    % CEA
    [total_cost, total_qalyloss, incremental_cost, incremental_qaly, ICER] = CEA_vacc(params);
    CEAs(i, :) = [total_cost(1), total_qalyloss(1), total_cost(2), total_qalyloss(2), incremental_cost, incremental_qaly, ICER];
    
end

% Calcuate mean and 95% interval
ICER_mean = mean(CEAs(:, 7));
ICER_left95 = quantile(CEAs(:, 7), 0.025);
ICER_right95 = quantile(CEAs(:, 7), 0.975);
fprintf('ICER mean: %f\nICER 95%% interval: (%f, %f)\n', ICER_mean, ICER_left95, ICER_right95)

% Generate table
parameter_names = params(isTested, 1)';
CEA_table = array2table([test_sample, CEAs], 'VariableNames', [parameter_names, ...
    {'No vaccination: total costs', 'No vaccination: total QALY loss', ...
    'Strategy A: total costs', 'Strategy A: total QALY loss', ...
    'Incremental costs', 'Incremental QALYs', 'ICER'}], 'RowNames', RowNames)

% Plot CEA plane
color = get(gcf, 'defaultAxesColorOrder');

figure(2)
hold on
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5, 'Color', color(2, :))
h = plot(linspace(20, 28), 5000 * linspace(20, 28), 'k');
hold off
legend(h, 'Threshold')
xlabel('Incremental QALYs')
ylabel('Incremental costs')
title('Cost-effectiveness plane (Normal distribution)')
saveas(gcf, 'CEA_plane_SA_normal', 'png')

figure(4)
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5)

%% Other Choice of Distribution (Lognormal)

for i = 1:n_parameter
    test_min = test_range(i, 1);
    test_max = test_range(i, 2);
    test_mean = (test_min + test_max) / 2;
    test_sample(:, i) = lognrnd(log(test_mean), 1, n_sample, 1);
end

% Initialize
CEAs = zeros(n_sample, 7);

% Loop over parameters
for i = 1:n_sample
    
    % Replace value
    params(isTested, 2) = num2cell(test_sample(i, :));
    
    % CEA
    [total_cost, total_qalyloss, incremental_cost, incremental_qaly, ICER] = CEA_vacc(params);
    CEAs(i, :) = [total_cost(1), total_qalyloss(1), total_cost(2), total_qalyloss(2), incremental_cost, incremental_qaly, ICER];
    
end

% Calcuate mean and 95% interval
ICER_mean = mean(CEAs(:, 7));
ICER_left95 = quantile(CEAs(:, 7), 0.025);
ICER_right95 = quantile(CEAs(:, 7), 0.975);
fprintf('ICER mean: %f\nICER 95%% interval: (%f, %f)\n', ICER_mean, ICER_left95, ICER_right95)

% Generate table
parameter_names = params(isTested, 1)';
CEA_table = array2table([test_sample, CEAs], 'VariableNames', [parameter_names, ...
    {'No vaccination: total costs', 'No vaccination: total QALY loss', ...
    'Strategy A: total costs', 'Strategy A: total QALY loss', ...
    'Incremental costs', 'Incremental QALYs', 'ICER'}], 'RowNames', RowNames)

% Plot CEA plane
figure(3)
hold on
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5, 'Color', color(3, :))
h = plot(linspace(0, 250), 5000 * linspace(0, 250), 'k');
hold off
legend(h, 'Threshold')
xlabel('Incremental QALYs')
ylabel('Incremental costs')
title('Cost-effectiveness plane (Lognormal distribution)')
saveas(gcf, 'CEA_plane_SA_lognormal', 'png')

figure(4)
plot(CEAs(:, 6), CEAs(:, 5), 'o', 'MarkerSize', 5)
plot(linspace(0, 250), 5000 * linspace(0, 250), 'k')
hold off
legend('Uniform', 'Normal', 'Lognormal', 'Threshold')
xlabel('Incremental QALYs')
ylabel('Incremental costs')
title('Cost-effectiveness plane')
saveas(gcf, 'CEA_plane_SA_all', 'png')

%% CEA Function

function [total_cost, total_qalyloss, incremental_cost, incremental_qaly, ICER] = CEA_vacc(params)

% Assign parameters
for i = 1:size(params, 1)
    cmd_string = sprintf('%s_ = params{%d, 2};', params{i, 1}, i);
    eval(cmd_string)
end

% Strategy
strategy = {'no_vacc', 'A'};
n_strategy = length(strategy);

% Time
tspan = t_init_:dt_:t_termi_;
nt = length(tspan);

% Initialize
total_clinical_case_y = zeros(1, n_strategy);
total_clinical_case_o = zeros(1, n_strategy);
total_clinical_case = zeros(1, n_strategy);
total_death_y = zeros(1, n_strategy);
total_death_o = zeros(1, n_strategy);
total_death = zeros(1, n_strategy);
num_vacc = zeros(1, n_strategy);

% Loop over scenario
for i = 1:n_strategy
    
    % beta
    denominator = [N_y_, N_o_; N_y_, N_o_];
    WAIFW = transmission_coeff_ * (mixing_matrix_ ./ denominator);
    beta_yy = WAIFW(1, 1);
    beta_yo = WAIFW(1, 2);
    beta_oy = WAIFW(2, 1);
    beta_oo = WAIFW(2, 2);
    
    % Coverage
    eval(sprintf('cov_y = coverage_y.%s_;', strategy{i}))
    eval(sprintf('cov_o = coverage_o.%s_;', strategy{i}))
    
    % Initial values
    R_y0 = N_y_ * prop_immune_y_ + N_y_ * (1 - prop_immune_y_) * cov_y * efficacy_;
    I_y0 = 5;
    S_y0 = N_y_ - I_y0 - R_y0;
    R_o0 = N_o_ * prop_immune_o_ + N_o_ * (1 - prop_immune_o_) * cov_o * efficacy_;
    I_o0 = 5;
    S_o0 = N_o_ - I_o0 - R_o0;
    
    % Initialize
    S_y = zeros(nt, 1);
    I_y = zeros(nt, 1);
    R_y = zeros(nt, 1);
    S_o = zeros(nt, 1);
    I_o = zeros(nt, 1);
    R_o = zeros(nt, 1);
    
    S_y(1) = S_y0;
    I_y(1) = I_y0;
    R_y(1) = R_y0;
    S_o(1) = S_o0;
    I_o(1) = I_o0;
    R_o(1) = R_o0;
    
    %% Solve the Model
    
    % Loop over time
    for j = 2:nt
        
        S_y(j) = S_y(j - 1) + (-beta_yy * I_y(j - 1) * S_y(j - 1) - beta_yo * I_o(j - 1) * S_y(j - 1)) * dt_;
        I_y(j) = I_y(j - 1) + (beta_yy * I_y(j - 1) * S_y(j - 1) + beta_yo * I_o(j - 1) * S_y(j - 1) - gamma_ * I_y(j - 1)) * dt_;
        R_y(j) = R_y(j - 1) + (gamma_ * I_y(j - 1)) * dt_;
        
        S_o(j) = S_o(j - 1) + (-beta_oy * I_y(j - 1) * S_o(j - 1) - beta_oo * I_o(j - 1) * S_o(j - 1)) * dt_;
        I_o(j) = I_o(j - 1) + (beta_oy * I_y(j - 1) * S_o(j - 1) + beta_oo * I_o(j - 1) * S_o(j - 1) - gamma_ * I_o(j - 1)) * dt_;
        R_o(j) = R_o(j - 1) + (gamma_ * I_o(j - 1)) * dt_;
    end
    
    %% Compute Incidence
    
    % Number of clinical cases
    new_inf_y = S_y(1:end - 1) - S_y(2:end);
    new_inf_o = S_o(1:end - 1) - S_o(2:end);
    clinical_case_y = prop_infection_clinical_y_ * new_inf_y;
    clinical_case_o = prop_infection_clinical_o_ * new_inf_o;
    total_clinical_case_y(i) = sum(clinical_case_y);
    total_clinical_case_o(i) = sum(clinical_case_o);
    total_clinical_case(i) = sum(clinical_case_y) + sum(clinical_case_o);
    
    % Number of deaths
    death_y = case_fatality_y_ * clinical_case_y;
    death_o = case_fatality_o_ * clinical_case_o;
    total_death_y(i) = sum(death_y);
    total_death_o(i) = sum(death_o);
    total_death(i) = sum(death_y) + sum(death_o);
    
    % Number of vaccinated
    num_vacc(i) = N_y_ * (1 - prop_immune_y_) * cov_y + N_o_ * (1 - prop_immune_o_) * cov_o;
    
end

% Total costs
total_vacc_cost = cost_per_vacc_ * num_vacc;
total_treatment_cost = cost_per_clinical_case_ * total_clinical_case;
total_cost = total_vacc_cost + total_treatment_cost;

% Total QALY losses
total_qalyloss_clinical_case = qalyloss_per_clinical_case_ * total_clinical_case;
total_qalyloss_death = qalyloss_per_death_ * total_death;
total_qalyloss = total_qalyloss_clinical_case + total_qalyloss_death;

% ICER
incremental_cost = total_cost(2) - total_cost(1);
incremental_qaly = total_qalyloss(1) - total_qalyloss(2);
ICER = incremental_cost / incremental_qaly;

end