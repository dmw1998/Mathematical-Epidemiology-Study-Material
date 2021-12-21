close all; clear all; clc

%% 3. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 1;               % Transmission Rate
b = 0;                  % Birth Rate
gamma = 0.5;            % Recover Rate
I0 = 2;                 % Initial Infectious Inidvidual

%% Set other parameter
total_time = 14;        % To 14 for each run
number_simul = 500;     % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time
S0 = N - I0;            % Initial susceptible individuals

%% Define transition probility
p_11 = @(s,i) beta/N*i*s*dt;            % (-1,1) one from S to I
p0_1 = @(s,i) gamma*i*dt;               % (0,-1) one from I to R
p1_1 = @(s,i) b*i*dt;                   % (1,-1) one born in S while one death in I
p10 = @(s,i) b*(N-s-i)*dt;              % (1,0)  one born in S while one death in R
p00 = @(s,i) 1 - p_11(s,i) - p0_1(s,i) - p1_1(s,i) - p10(s,i);          % Unchanged

%% Simulation
n = total_time/dt;
data_coll_S = zeros(n, number_simul);
data_coll_I = zeros(n, number_simul);
data_coll_S(1,:) = N - I0;
data_coll_I(1,:) = I0;

for k = 1 : number_simul
    for j = 2 : n
        s = data_coll_S(j-1,k);         % # of S
        i = data_coll_I(j-1,k);         % # of I

        if i == 0                       % Infective individual is zero
            data_coll_S(j,k) = s;
            data_coll_I(j,k) = i;
            continue                    % Nothing change
        elseif s == 0                   % Susceptible individual is zero
            if i == N                   % All are infective
                a1 = 0;                 % No S to I
                a2 = p0_1(s,i);         % One I to R
                a3 = a2 + p1_1;         % New born in S and one death in I
                a4 = a3;                % No death in R
            else                        % Some people in R == death in R
                a1 = 0;                 % No S to I
                a2 = p0_1(s,i);         % One I to R
                a3 = a2 + p0_1(s,i);    % New born in S and one death in I
                a4 = a3 + p1_1(s,i);    % New born in S and one death in R
            end
        else                            % Both S and I are not zero
            a1 = p_11(s,i);              % Less than a1, one from S to I
            a2 = a1 + p0_1(s,i);        % Greater than a1 but less than a2, one from I to R
            a3 = a2 + p1_1(s,i);        % Greater than a2 but less than a3, (1,-1)
            a4 = a3 + p10(s,i);         % Greater than a3 but less than a4, (1,0)
        end                             % Greater than a4, not changed

        c = rand;                       % Pick a random number between (0,1)
        if c < a1                       % Less than a1, one from S to I
            data_coll_S(j,k) = s-1;
            data_coll_I(j,k) = i+1;
        elseif c > a1 && c < a2         % Greater than a1 but less than a2, one from I to R
            data_coll_I(j,k) = i-1;
        elseif c > a2 && c < a3         %  Greater than a2 but less than a3, (1,-1)
            data_coll_S(j,k) = s+1;
            data_coll_I(j,k) = i-1;
        elseif c > a3 && c < a4         % Greater than a3 but less than a4, (1,0)
            data_coll_S(j,k) = s+1;
        else                            % Greater than a4, unchanged
            data_coll_S(j,k) = s;
            data_coll_I(j,k) = i;
        end

    end
end

figure
plot(data_coll_I)
hold on

average_rel = mean(data_coll_I,2);
plot(average_rel, 'k', 'LineWidth', 1.5)

%% 4. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 1;               % Transmission Rate
b = 0.25;               % Birth Rate
gamma = 0.25;           % Recover Rate
I0 = 2;                 % Initial Infectious Inidvidual

%% Set other parameter
total_time = 14;        % To 14 for each run
number_simul = 500;     % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time
S0 = N - I0;            % Initial susceptible individuals

%% Define transition probility
p_11 = @(s,i) beta/N*i*s*dt;            % (-1,1) one from S to I
p0_1 = @(s,i) gamma*i*dt;               % (0,-1) one from I to R
p1_1 = @(s,i) b*i*dt;                   % (1,-1) one born in S while one death in I
p10 = @(s,i) b*(N-s-i)*dt;              % (1,0)  one born in S while one death in R
p00 = @(s,i) 1 - p_11(s,i) - p0_1(s,i) - p1_1(s,i) - p10(s,i);          % Unchanged

%% Simulation
n = total_time/dt;
data_coll_S = zeros(n, number_simul);
data_coll_I = zeros(n, number_simul);
data_coll_S(1,:) = N - I0;
data_coll_I(1,:) = I0;

for k = 1 : number_simul
    for j = 2 : n
        s = data_coll_S(j-1,k);         % # of S
        i = data_coll_I(j-1,k);         % # of I

        if i == 0                       % Infective individual is zero
            data_coll_S(j,k) = s;
            data_coll_I(j,k) = i;
            continue                    % Nothing change
        elseif s == 0                   % Susceptible individual is zero
            if i == N                   % All are infective
                a1 = 0;                 % No S to I
                a2 = p0_1(s,i);         % One I to R
                a3 = a2 + p1_1;         % New born in S and one death in I
                a4 = a3;                % No death in R
            else                        % Some people in R == death in R
                a1 = 0;                 % No S to I
                a2 = p0_1(s,i);         % One I to R
                a3 = a2 + p0_1(s,i);    % New born in S and one death in I
                a4 = a3 + p1_1(s,i);    % New born in S and one death in R
            end
        else                            % Both S and I are not zero
            a1 = p_11(s,i);              % Less than a1, one from S to I
            a2 = a1 + p0_1(s,i);        % Greater than a1 but less than a2, one from I to R
            a3 = a2 + p1_1(s,i);        % Greater than a2 but less than a3, (1,-1)
            a4 = a3 + p10(s,i);         % Greater than a3 but less than a4, (1,0)
        end                             % Greater than a4, not changed

        c = rand;                       % Pick a random number between (0,1)
        if c < a1                       % Less than a1, one from S to I
            data_coll_S(j,k) = s-1;
            data_coll_I(j,k) = i+1;
        elseif c > a1 && c < a2         % Greater than a1 but less than a2, one from I to R
            data_coll_I(j,k) = i-1;
        elseif c > a2 && c < a3         %  Greater than a2 but less than a3, (1,-1)
            data_coll_S(j,k) = s+1;
            data_coll_I(j,k) = i-1;
        elseif c > a3 && c < a4         % Greater than a3 but less than a4, (1,0)
            data_coll_S(j,k) = s+1;
        else                            % Greater than a4, unchanged
            data_coll_S(j,k) = s;
            data_coll_I(j,k) = i;
        end

    end
end

figure
plot(data_coll_I)
hold on

average_rel = mean(data_coll_I,2);
plot(average_rel, 'k', 'LineWidth', 1.5)