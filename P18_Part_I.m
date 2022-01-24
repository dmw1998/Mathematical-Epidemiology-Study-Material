clear all; clc
%% Parameters
% Young
N_y = 10000;            % Young population
p_yr = 0.25;            % Proportion of immune children
p_yc = 0.4;             % Proportion of clinical case among children
                        % Only part of the population go to clinic, rest recover themselves
d_y = 0.0001;           % case fatality rate for children in clinic
I_y0 = 5;               % Initial young infectious individuals

% Old
N_o = 40000;            % Young population
p_or = 0.4;             % Proportion of immune children
p_oc = 0.4;             % Proportion of clinical case among children
                        % Only part of the population go to clinic, rest recover themselves
d_o = 0.0002;           % case fatality rate for children in clinic
I_o0 = 5;               % Initial old infectious individuals

% Common coefficients
trans_coeff = 0.55;     % Transmission coefficient
mixing_mat = [1, 0.15; 0.15, 0.8];
popu_mat = [N_y, N_o; N_y, N_o];
WAIFW = trans_coeff*(mixing_mat./popu_mat);
beta_yy = WAIFW(1,1);
beta_yo = WAIFW(1,2);
beta_oy = WAIFW(2,1);
beta_oo = WAIFW(2,2);

gamma = 1/3;            % Recover rate (for all person)

e = 0.4;                % Vaccine efficacy
% For Plan A
c_yA = 0.5;             % Vaccine coverage for children in Plan A
c_oA = 0;               % Vaccine coverage for old in Plan A
% For Plan B
c_yB = 0.25;            % Vaccine coverage for children in Plan B
c_oB = 0.25;            % Vaccine coverage for old in Plan B

% Define ODE
dS_y = @(s_y,i_y,i_o,dt) s_y + (-beta_yy*i_y*s_y - beta_yo*i_o*s_y)*dt;
dI_y = @(s_y,i_y,i_o,dt) i_y + (beta_yy*i_y*s_y + beta_yo*i_o*s_y - gamma*i_y)*dt;
% dR_y = @(i_y,dt) (gamma*i_y)*dt;

dS_o = @(s_o,i_y,i_o,dt) s_o + (-beta_oy*i_y*s_o - beta_oo*i_o*s_o)*dt;
dI_o = @(s_o,i_y,i_o,dt) i_o + (beta_oy*i_y*s_o + beta_oo*i_o*s_o - gamma*i_o)*dt;
% dR_o = @(i_o,dt) (gamma*i_o)*dt;

% For using forward Euler method
start_time = 0;         % Initial time, start from 0
end_time = 200;         % Terminal time, end at 200 days
dt = 0.2;               % Time step size

% Cost and OALY
vac_cost = 8;           % Cost for per vaccination
cli_cost = 12;          % Cost for clinical case
cli_loss = 0.005;       % QALY loss per clinical case
death_loss = 30;        % QALY loss per death

%% No vaccine
S_y0 = (1-p_yr)*(N_y-I_y0);
S_o0 = (1-p_or)*(N_o-I_o0);

% Simulation
n = 1;
S_y(n) = S_y0;
I_y(n) = I_y0;

S_o(n) = S_o0;
I_o(n) = I_o0;

s_y= S_y(n);
i_y = I_y(n);
s_o = S_o(n);
i_o = I_o(n);
new_I_y = [];
new_I_o = [];
for n = start_time+2:end_time/dt
    S_y(n) = dS_y(s_y,i_y,i_o,dt);
    I_y(n) = dI_y(s_y,i_y,i_o,dt);

    S_o(n) = dS_o(s_o,i_y,i_o,dt);
    I_o(n) = dI_o(s_o,i_y,i_o,dt);

    new_I_y(n-1) = s_y - S_y(n);
    new_I_o(n-1) = s_o - S_o(n);

    s_y= S_y(n);
    i_y = I_y(n);
    s_o = S_o(n);
    i_o = I_o(n);

end

total_vac_cost = 0;
total_cli = sum(new_I_y*p_yc + new_I_o*p_oc);
total_cli_cost = total_cli*cli_cost;
total_cost = total_vac_cost + total_cli_cost;

total_QALY_cli_loss = total_cli*cli_loss;
total_death = sum(new_I_y*p_yc*d_y) + sum(new_I_o*p_oc*d_o);
total_QALY_death_loss = total_death*death_loss;
total_QALY_loss = total_QALY_cli_loss + total_QALY_death_loss;

%% Plan A
% Initial situation
R_y0 = p_yr*N_y + (1-p_yr)*(N_y-I_y0)*c_yA*e;      % Recoverd Children
S_y0 = N_y - I_y0 - R_y0;                   % Suscepible Children

R_o0 = p_or*N_o + (1-p_or)*(N_o-I_o0)*c_oA*e;      % Recoverd Adults
S_o0 = N_o - I_o0 - R_o0;                   % Susceptible Adults

% Simulation
n = 1;
S_y(n) = S_y0;
I_y(n) = I_y0;

S_o(n) = S_o0;
I_o(n) = I_o0;

s_y= S_y(n);
i_y = I_y(n);
s_o = S_o(n);
i_o = I_o(n);
new_I_y = [];
new_I_o = [];
for n = start_time+2:end_time/dt
    S_y(n) = dS_y(s_y,i_y,i_o,dt);
    I_y(n) = dI_y(s_y,i_y,i_o,dt);

    S_o(n) = dS_o(s_o,i_y,i_o,dt);
    I_o(n) = dI_o(s_o,i_y,i_o,dt);

    new_I_y(n-1) = s_y - S_y(n);
    new_I_o(n-1) = s_o - S_o(n);

    s_y= S_y(n);
    i_y = I_y(n);
    s_o = S_o(n);
    i_o = I_o(n);

end

total_vac_A = S_y0*c_yA+ S_o0*c_oA;
total_vac_cost_A = total_vac_A*vac_cost;
total_cli_A = sum(new_I_y*p_yc + new_I_o*p_oc);
total_cli_cost_A = total_cli_A*cli_cost;
total_cost_A = total_vac_cost_A + total_cli_cost_A;

total_QALY_cli_loss_A = total_cli_A*cli_loss;
total_death_A = sum(new_I_y*p_yc*d_y) + sum(new_I_o*p_oc*d_o);
total_QALY_death_loss_A = total_death_A*death_loss;
total_QALY_loss_A = total_QALY_cli_loss_A + total_QALY_death_loss_A;

%% Plan B
% Initial situation
R_y0 = p_yr*N_y + (1-p_yr)*(N_y-I_y0)*c_yB*e;      % Recoverd Children
S_y0 = N_y - I_y0 - R_y0;                   % Suscepible Children

R_o0 = p_or*N_o + (1-p_or)*(N_o-I_o0)*c_oB*e;      % Recoverd Adults
S_o0 = N_o - I_o0 - R_o0;                   % Susceptible Adults

% Simulation
n = 1;
S_y(n) = S_y0;
I_y(n) = I_y0;

S_o(n) = S_o0;
I_o(n) = I_o0;

s_y= S_y(n);
i_y = I_y(n);
s_o = S_o(n);
i_o = I_o(n);
new_I_y = [];
new_I_o = [];
for n = start_time+2:end_time/dt
    S_y(n) = dS_y(s_y,i_y,i_o,dt);
    I_y(n) = dI_y(s_y,i_y,i_o,dt);

    S_o(n) = dS_o(s_o,i_y,i_o,dt);
    I_o(n) = dI_o(s_o,i_y,i_o,dt);

    new_I_y(n-1) = s_y - S_y(n);
    new_I_o(n-1) = s_o - S_o(n);

    s_y= S_y(n);
    i_y = I_y(n);
    s_o = S_o(n);
    i_o = I_o(n);

end

total_vac_B = S_y0*c_yB + S_o0*c_oB;
total_vac_cost_B = total_vac_B*vac_cost;
total_cli_B = sum(new_I_y*p_yc + new_I_o*p_oc);
total_cli_cost_B = total_cli_B*cli_cost;
total_cost_B = total_vac_cost_B + total_cli_cost_B;

total_QALY_cli_loss_B = total_cli_B*cli_loss;
total_death_B = sum(new_I_y*p_yc*d_y) + sum(new_I_o*p_oc*d_o);
total_QALY_death_loss_B = total_death_B*death_loss;
total_QALY_loss_B = total_QALY_cli_loss_B + total_QALY_death_loss_B;

%% Prevented
% Clinical cases prevented
cli_prev_A = total_cli - total_cli_A
cli_prev_B = total_cli - total_cli_B

% Deaths prevented
death_prev_A = total_death - total_death_A
death_prev_B = total_death - total_death_B

% Number vaccinated
total_vac_A
total_vac_B

% Vaccine costs
vac_cost_prev_A = total_vac_cost_A
vac_cost_prev_B = total_vac_cost_B

% Treatment costs
cli_cost_prev_A = total_cli_cost - total_cli_cost_A
cli_cost_prev_B = total_cli_cost - total_cli_cost_B

% Total costs
total_cost_prev_A = total_cost - total_cost_A
total_cost_prev_B = total_cost - total_cost_B

% QALYs gained by preventing cases
QALY_cli_gained_A = total_QALY_cli_loss - total_QALY_cli_loss_A
QALY_cli_gained_B = total_QALY_cli_loss - total_QALY_cli_loss_B

% QALYs gained by preventing deaths
QALY_death_gained_A = total_QALY_death_loss - total_QALY_death_loss_A
QALY_death_gained_B = total_QALY_death_loss - total_QALY_death_loss_B

% Total QALYs gained
total_QALY_gained_A = total_QALY_loss - total_QALY_loss_A
total_QALY_gained_B = total_QALY_loss - total_QALY_loss_B

% Incremental cost per case prevented
inc_cost_case_A = total_cost_prev_A / cli_prev_A
inc_cost_case_B = total_cost_prev_B / cli_prev_B

% Incremental cost per death prevented
inc_cost_death_A = total_cost_prev_A / death_prev_A
inc_cost_death_B = total_cost_prev_B / death_prev_B

% Incremental cost per QALY gained
inc_cost_QALY_A = total_cost_prev_A / total_QALY_gained_A
inc_cost_QALY_B = total_cost_prev_B / total_QALY_gained_B