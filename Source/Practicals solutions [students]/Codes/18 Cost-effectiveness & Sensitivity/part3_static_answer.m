close all
clear
restoredefaultpath

%% Settings

% Strategy
strategy = {'no_vacc', 'A', 'B'};
n_strategy = length(strategy);

% Young
N_y = 10000;
prop_immune_y = 0.25;
prop_infection_clinical_y = 0.4;
case_fatality_y = 0.0001;

% Old
N_o = 40000;
prop_immune_o = 0.4;
prop_infection_clinical_o = 0.4;
case_fatality_o = 0.0002;

% Transmission
transmission_coeff = 0.55;
mixing_matrix = [1, 0.15; 0.15, 0.8];
denominator = [N_y, N_o; N_y, N_o];
WAIFW = transmission_coeff * (mixing_matrix ./ denominator);

% Progression
gamma = 1 / 3;

% Vaccination
efficacy = 0.4;
coverage_y.no_vacc = 0;
coverage_o.no_vacc = 0;
coverage_y.A = 0.5;
coverage_o.A = 0;
coverage_y.B = 0.25;
coverage_o.B = 0.25;

% Time
dt = 0.2;
t_init = 0;
t_termi = 200;
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
    beta_yy = WAIFW(1, 1);
    beta_yo = WAIFW(1, 2);
    beta_oy = WAIFW(2, 1);
    beta_oo = WAIFW(2, 2);
    
    % Coverage
    eval(sprintf('cov_y = coverage_y.%s;', strategy{i}))
    eval(sprintf('cov_o = coverage_o.%s;', strategy{i}))
    
    % Initial values
    R_y0 = N_y * prop_immune_y + N_y * (1 - prop_immune_y) * cov_y * efficacy;
    I_y0 = 5;
    S_y0 = N_y - I_y0 - R_y0;
    R_o0 = N_o * prop_immune_o + N_o * (1 - prop_immune_o) * cov_o * efficacy;
    I_o0 = 5;
    S_o0 = N_o - I_o0 - R_o0;
    
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
        
        S_y(j) = S_y(j - 1) + (-beta_yy * I_y(j - 1) * S_y(j - 1) - beta_yo * I_o(j - 1) * S_y(j - 1)) * dt;
        I_y(j) = I_y(j - 1) + (beta_yy * I_y(j - 1) * S_y(j - 1) + beta_yo * I_o(j - 1) * S_y(j - 1) - gamma * I_y(j - 1)) * dt;
        R_y(j) = R_y(j - 1) + (gamma * I_y(j - 1)) * dt;
        
        S_o(j) = S_o(j - 1) + (-beta_oy * I_y(j - 1) * S_o(j - 1) - beta_oo * I_o(j - 1) * S_o(j - 1)) * dt;
        I_o(j) = I_o(j - 1) + (beta_oy * I_y(j - 1) * S_o(j - 1) + beta_oo * I_o(j - 1) * S_o(j - 1) - gamma * I_o(j - 1)) * dt;
        R_o(j) = R_o(j - 1) + (gamma * I_o(j - 1)) * dt;
    end
    
    %% Compute Incidence
    
    % Number of clinical cases
    new_inf_y = S_y(1:end - 1) - S_y(2:end);
    new_inf_o = S_o(1:end - 1) - S_o(2:end);
    clinical_case_y = prop_infection_clinical_y * new_inf_y;
    clinical_case_o = prop_infection_clinical_o * new_inf_o;
    total_clinical_case_y(i) = sum(clinical_case_y);
    total_clinical_case_o(i) = sum(clinical_case_o);
    total_clinical_case(i) = sum(clinical_case_y) + sum(clinical_case_o);
    
    % Number of deaths
    death_y = case_fatality_y * clinical_case_y;
    death_o = case_fatality_o * clinical_case_o;
    total_death_y(i) = sum(death_y);
    total_death_o(i) = sum(death_o);
    total_death(i) = sum(death_y) + sum(death_o);
    
    % Number of vaccinated
    num_vacc(i) = N_y * (1 - prop_immune_y) * cov_y + N_o * (1 - prop_immune_o) * cov_o;
    
end

%% Total Costs and QALYs

% Costs and QALYs
cost_per_vacc = 8;
cost_per_clinical_case = 12;
qalyloss_per_clinical_case = 0.005;
qalyloss_per_death = 30;

% Total costs
total_vacc_cost = cost_per_vacc * num_vacc;
total_treatment_cost = cost_per_clinical_case * total_clinical_case;
total_cost = total_vacc_cost + total_treatment_cost;

% Total QALY losses
total_qalyloss_clinical_case = qalyloss_per_clinical_case * total_clinical_case;
total_qalyloss_death = qalyloss_per_death * total_death;
total_qalyloss = total_qalyloss_clinical_case + total_qalyloss_death;

% Generate total table
RowNames = {'Clinical cases', 'Deaths', 'Number vaccinated', ...
    'Vaccine costs', 'Treatment costs', 'Total costs', ...
    'QALY losses for clinical cases', 'QALY losses for deaths', 'Total QALY losses'};
total_table = array2table([total_clinical_case; total_death; num_vacc;
    total_vacc_cost; total_treatment_cost; total_cost;
    total_qalyloss_clinical_case; total_qalyloss_death; total_qalyloss], ...
    'VariableNames', strategy_names, 'RowNames', RowNames)

%% Cost-Effectiveness Analysis: Dynamic Model

clinical_case_prevented = [total_clinical_case(1) - total_clinical_case(2),...
    total_clinical_case(1) - total_clinical_case(3)];

death_prevented = [total_death(1) - total_death(2), ...
    total_death(1) - total_death(3)];

number_vaccinated = [num_vacc(2) - num_vacc(1), ...
    num_vacc(3) - num_vacc(1)];

vaccine_cost = [total_vacc_cost(2) - total_vacc_cost(1), ...
    total_vacc_cost(3) - total_vacc_cost(1)];

treatment_cost = [total_treatment_cost(2) - total_treatment_cost(1), ...
    total_treatment_cost(3) - total_treatment_cost(1)];

total_cost_spent = vaccine_cost + treatment_cost;

qaly_gained_clinical_case = [total_qalyloss_clinical_case(1) - total_qalyloss_clinical_case(2), ...
    total_qalyloss_clinical_case(1) - total_qalyloss_clinical_case(3)];

qaly_gained_death = [total_qalyloss_death(1) - total_qalyloss_death(2), ...
    total_qalyloss_death(1) - total_qalyloss_death(3)];

total_qaly_gained = qaly_gained_clinical_case + qaly_gained_death;

incremental_cost_per_case_prevented = total_cost_spent ./ clinical_case_prevented;
incremental_cost_per_death_prevented = total_cost_spent ./ death_prevented;
incremental_cost_per_qaly_gained = total_cost_spent ./ total_qaly_gained;

% Generate CEA table
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

%% Cost-Effectiveness Analysis: Static Model

total_clinical_case_novac = total_clinical_case(1);
total_clinical_case_A = total_clinical_case_y(1) * (1 - efficacy * coverage_y.A) ...
    + total_clinical_case_o(1) * (1 - efficacy * coverage_o.A);
total_clinical_case_static = [total_clinical_case_novac, total_clinical_case_A, total_clinical_case_novac - total_clinical_case_A];

total_death_novac = total_death(1);
total_death_A = total_death_y(1) * (1 - efficacy * coverage_y.A) ...
    + total_death_o(1) * (1 - efficacy * coverage_o.A);
total_death_static = [total_death_novac, total_death_A, total_death_novac - total_death_A];

number_vaccinated_static = [0, number_vaccinated(1), number_vaccinated(1)];

vaccine_cost_static = number_vaccinated_static * cost_per_vacc;

treatment_cost_static = total_clinical_case_static * cost_per_clinical_case;
treatment_cost_static(3) = -treatment_cost_static(3);

total_cost_static = vaccine_cost_static + treatment_cost_static;

qaly_loss_clinical_case_static = total_clinical_case_static * qalyloss_per_clinical_case;

qaly_loss_death_static = total_death_static * qalyloss_per_death;

total_qaly_loss_static = qaly_loss_clinical_case_static + qaly_loss_death_static;

incremental_cost_per_case_prevented_static = [0, 0, total_cost_static(3) ./ total_clinical_case_static(3)];
incremental_cost_per_death_prevented_static = [0, 0, total_cost_static(3) ./ total_death_static(3)];
incremental_cost_per_qaly_gained_static = [0, 0, total_cost_static(3) ./ total_qaly_loss_static(3)];

% Generate CEA table
VariableNames = {'No vaccination', 'Strategy A', 'A compared to no vaccination'};
RowNames = {'Clinical cases', 'Deaths', 'Number vaccinated', ...
    'Vaccine costs', 'Treatment costs', 'Total costs', ...
    'QALY loss for clinical cases', 'QALY loss for deaths', 'Total QALY loss', ...
    'Incremental cost per case prevented', 'Incremental cost per daeth prevented', 'Incremental cost per QALY gained'};
static_table = array2table([total_clinical_case_static; total_death_static; number_vaccinated_static;
    vaccine_cost_static; treatment_cost_static; total_cost_static;
    qaly_loss_clinical_case_static; qaly_loss_death_static; total_qaly_loss_static;
    incremental_cost_per_case_prevented_static; incremental_cost_per_death_prevented_static; incremental_cost_per_qaly_gained_static],...
    'VariableName', VariableNames, 'RowNames', RowNames)