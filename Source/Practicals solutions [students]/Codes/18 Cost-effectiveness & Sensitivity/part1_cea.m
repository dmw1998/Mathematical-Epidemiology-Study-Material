close all
clear
restoredefaultpath

%% Settings

% Strategy
strategy = {'no_vacc', 'A', 'B'};
n_strategy = length(strategy);

% Young
N_y = nan;
prop_immune_y = nan;
prop_infection_clinical_y = nan;
case_fatality_y = nan;

% Old
N_o = nan;
prop_immune_o = nan;
prop_infection_clinical_o = nan;
case_fatality_o = nan;

% Transmission
transmission_coeff = nan;
mixing_matrix = nan(2, 2);
denominator = nan(2, 2);
WAIFW = nan(2, 2);

% Progression
gamma = nan;

% Vaccination
efficacy = nan;
coverage_y.no_vacc = nan;
coverage_o.no_vacc = nan;
coverage_y.A = nan;
coverage_o.A = nan;
coverage_y.B = nan;
coverage_o.B = nan;

% Time
dt =nan;
t_init = nan;
t_termi = nan;
tspan = t_init:dt:t_termi;
nt = length(tspan);

%% Model Simulation

% Initialize
total_clinical_case_y = zeros(1, n_strategy);
total_clinical_case_o = zeros(1, n_strategy);
total_clinical_case = zeros(1, n_strategy);
total_death_y = zeros(1, n_strategy);
total_death_o = zeros(1, n_strategy);
total_death = zeros(1, n_strategy);
num_vacc = zeros(1, n_strategy);

strategy_names = {'No vaccination', 'Strategy A', 'Strategy B'};

% Loop over scenario
for i = 1:n_strategy
    
    % beta
    beta_yy = nan;
    beta_yo = nan;
    beta_oy = nan;
    beta_oo = nan;
    
    % Coverage
    eval(sprintf('cov_y = coverage_y.%s;', strategy{i}))
    eval(sprintf('cov_o = coverage_o.%s;', strategy{i}))
    
    % Initial values
    R_y0 = nan;
    I_y0 = nan;
    S_y0 = nan;
    R_o0 = nan;
    I_o0 = nan;
    S_o0 = nan;
    
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
        
        S_y(j) = S_y(j - 1) + nan;
        I_y(j) = I_y(j - 1) + nan;
        R_y(j) = R_y(j - 1) + nan;
        
        S_o(j) = S_o(j - 1) + nan;
        I_o(j) = I_o(j - 1) + nan;
        R_o(j) = R_o(j - 1) + nan;
    end
    
    %% Compute Incidence
    
    % Number of clinical cases
    new_inf_y = nan(nt - 1, 1);
    new_inf_o = nan(nt - 1, 1);
    clinical_case_y = nan(nt - 1, 1);
    clinical_case_o = nan(nt - 1, 1);
    total_clinical_case_y(i) = nan;
    total_clinical_case_o(i) = nan;
    total_clinical_case(i) = nan;
    
    % Number of deaths
    death_y = nan(nt - 1, 1);
    death_o = nan(nt - 1, 1);
    total_death_y(i) = nan;
    total_death_o(i) = nan;
    total_death(i) = nan;
    
    % Number of vaccinated
    num_vacc(i) = nan;
    
    %% Plot
    
    figure(i)
    hold on
    plot(tspan, I_y)
    plot(tspan, I_o)
    hold off
    legend('Child', 'Adult')
    ylim([0, 400])
    xlabel('Days')
    ylabel('Number of infectious')
    title(sprintf('%s', strategy_names{i}))
    saveas(gcf, sprintf('strategy_%s_infectious', strategy{i}), 'png')
    
end

%% Total Costs and QALYs

% Costs and QALYs
cost_per_vacc = nan;
cost_per_clinical_case = nan;
qalyloss_per_clinical_case = nan;
qalyloss_per_death = nan;

% Total costs
total_vacc_cost = nan(1, n_strategy);
total_treatment_cost = nan(1, n_strategy);
total_cost = nan(1, n_strategy);

% Total QALY losses
total_qalyloss_clinical_case = nan(1, n_strategy);
total_qalyloss_death = nan(1, n_strategy);
total_qalyloss = nan(1, n_strategy);

% Generate total table
RowNames = {'Clinical cases', 'Deaths', 'Number vaccinated', ...
    'Vaccine costs', 'Treatment costs', 'Total costs', ...
    'QALY losses for clinical cases', 'QALY losses for deaths', 'Total QALY losses'};
total_table = array2table([total_clinical_case; total_death; num_vacc;
    total_vacc_cost; total_treatment_cost; total_cost;
    total_qalyloss_clinical_case; total_qalyloss_death; total_qalyloss], ...
    'VariableNames', strategy_names, 'RowNames', RowNames)

%% Cost-Effectiveness Analysis

clinical_case_prevented = nan(1, 2);
death_prevented = nan(1, 2);
number_vaccinated = nan(1, 2);
vaccine_cost = nan(1, 2);
treatment_cost = nan(1, 2);
total_cost_spent = nan(1, 2);
qaly_gained_clinical_case = nan(1, 2);
qaly_gained_death = nan(1, 2);
total_qaly_gained = nan(1, 2);
incremental_cost_per_case_prevented = total_cost_spent ./ clinical_case_prevented;
incremental_cost_per_death_prevented = total_cost_spent ./ death_prevented;
incremental_cost_per_qaly_gained = total_cost_spent ./ total_qaly_gained;

% Generate difference table
VariableNames = {'A compared to no vaccination', 'B compared to no vaccination'};
RowNames = {'Clinical cases prevented', 'Deaths prevented', 'Number vaccinated', ...
    'Vaccine costs', 'Treatment costs', 'Total costs spent', ...
    'QALY gained by preventing clinical cases', 'QALY gained by preventing deaths', 'Total QALY gains', ...
    'Incremental cost per case prevented', 'Incremental cost per daeth prevented', 'Incremental cost per QALY gained'};
difference_table = array2table([clinical_case_prevented; death_prevented; number_vaccinated;
    vaccine_cost; treatment_cost; total_cost_spent;
    qaly_gained_clinical_case; qaly_gained_death; total_qaly_gained;
    incremental_cost_per_case_prevented; incremental_cost_per_death_prevented; incremental_cost_per_qaly_gained], ...
    'VariableName', VariableNames, 'RowNames', RowNames)

% Plot cost-effectiveness plane
figure(4)
hold on
plot(total_qaly_gained(1), total_cost_spent(1), '.', 'MarkerSize', 25)
plot(total_qaly_gained(2), total_cost_spent(2), '.', 'MarkerSize', 25)
plot(0:0.1:30, 5000 * (0:0.1:30), 'k')
hold off
xlim([0, 30])
ylim([-1e+4, 5e+4])
legend('Strategy A compared to no vaccination', 'Strategy B compared to no vaccination', ...
    'Threshold', 'Location', 'southwest')
xlabel('Incremental QALYs')
ylabel('Incremental costs')
title('Cost-Effectiveness Plane')
saveas(gcf, 'CEA_plane', 'png')