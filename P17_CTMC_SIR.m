close all; clear all; clc

%% 3. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 1;               % Transmission Rate
b = 0;                  % Birth Rate
gamma = 0.5;            % Recover Rate
I0 = 2;                 % Initial Infectious Inidvidual

%% Set other parameter
total_time = 40;        % To 40 for each run
number_simul = 50;      % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time
S0 = N - I0;            % Initial susceptible individuals

%% Define transition probility
p_11 = @(s,i,dt) beta/N*i*s*dt;            % (-1,1) one from S to I
p0_1 = @(s,i,dt) gamma*i*dt;               % (0,-1) one from I to R
p1_1 = @(s,i,dt) b*i*dt;                   % (1,-1) one born in S while one death in I
p10 = @(s,i,dt) b*(N-s-i)*dt;              % (1,0)  one born in S while one death in R
p00 = @(s,i,dt) 1 - p_11(s,i,dt) - p0_1(s,i,dt) - p1_1(s,i,dt) - p10(s,i,dt);          % Unchanged

%% Simulation
data_cell_I = {}; data_cell_S = {};
index_cell = {};
data_cell_I = cell(1, number_simul);
data_cell_S = cell(1, number_simul);
index_cell = cell(1,number_simul);

for k = 1 : number_simul
    data_coll_I = []; data_coll_S = [];
    data_coll_I(1) = I0;                    % Initialize the data set
    data_coll_S(1) = S0;                    % Initialize the data set
    
    % T_ind keep the time of people changed
    n = 1;
    T_ind = [];                             % Initialize the change time
    T_ind(1) = 0;                           % Start from time = 0
    while T_ind(n) < total_time
        s = data_coll_S(n);
        i = data_coll_I(n);
        % Pick a random number to decide the length of time interval
        % i.e. to give the next time which the change happends
        c = rand;

        n = n+1;                            % Next event
        % By given formular, give M (the total change rate) and T_len (length of time interval)
        M = 1 - p00(s,i,1);
        if M == 0
            T_len = 0.01;                   % Step one step
            data_coll_S(n) = data_coll_S(n-1);
            data_coll_I(n) = data_coll_I(n-1);
            T_ind(n) = T_ind(n-1) + T_len;
            continue
        else
            T_len = -log(c)/M;
            T_ind(n) = T_ind(n-1) + T_len;      % Next change will happen at this time
        end

        % Determine the change interval
        a1 = p_11(s,i,T_len)/(1 - p00(s,i,T_len));
        a2 = a1 + p0_1(s,i,T_len)/(1 - p00(s,i,T_len));
        a3 = a2 + p1_1(s,i,T_len)/(1 - p00(s,i,T_len));

        c = rand;                           % Pick a random number between (0,1)
        if c < a1                           % (-1,1)
            data_coll_S(n) = data_coll_S(n-1)-1;
            data_coll_I(n) = data_coll_I(n-1)+1;
        elseif c >= a1 && c < a2            % (0,-1)
            data_coll_S(n) = data_coll_S(n-1);
            data_coll_I(n) = data_coll_I(n-1)-1;
        elseif c >= a1 && c < a3            % (1,-1)
            data_coll_S(n) = data_coll_S(n-1)+1;
            data_coll_I(n) = data_coll_I(n-1)-1;
        else                                % (1,0)
            data_coll_S(n) = data_coll_S(n-1)+1;
            data_coll_I(n) = data_coll_I(n-1);
        end
    end

    index_cell{1,k} = T_ind;
    data_cell_S{1,k} = data_coll_S;
    data_cell_I{1,k} = data_coll_I;

end

figure
hold on
for k = 1 : number_simul
    index_coll = index_cell{1,k};
    data_coll_I = data_cell_I{1,k};
    plot(index_coll, data_coll_I)
end
xlim([0,40]);

%% 4. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 1;               % Transmission Rate
b = 0.25;               % Birth Rate
gamma = 0.25;           % Recover Rate
I0 = 2;                 % Initial Infectious Inidvidual

%% Define transition probility
p_11 = @(s,i,dt) beta/N*i*s*dt;            % (-1,1) one from S to I
p0_1 = @(s,i,dt) gamma*i*dt;               % (0,-1) one from I to R
p1_1 = @(s,i,dt) b*i*dt;                   % (1,-1) one born in S while one death in I
p10 = @(s,i,dt) b*(N-s-i)*dt;              % (1,0)  one born in S while one death in R
p00 = @(s,i,dt) 1 - p_11(s,i,dt) - p0_1(s,i,dt) - p1_1(s,i,dt) - p10(s,i,dt);          % Unchanged

%% Simulation
data_cell_I = {}; data_cell_S = {};
index_cell = {};
data_cell_I = cell(1, number_simul);
data_cell_S = cell(1, number_simul);
index_cell = cell(1,number_simul);

for k = 1 : number_simul
    data_coll_I = []; data_coll_S = [];
    data_coll_I(1) = I0;                    % Initialize the data set
    data_coll_S(1) = S0;                    % Initialize the data set
    
    % T_ind keep the time of people changed
    n = 1;
    T_ind = [];                             % Initialize the change time
    T_ind(1) = 0;                           % Start from time = 0
    while T_ind(n) < total_time
        s = data_coll_S(n);
        i = data_coll_I(n);
        % Pick a random number to decide the length of time interval
        % i.e. to give the next time which the change happends
        c = rand;

        n = n+1;                            % Next event
        % By given formular, give M (the total change rate) and T_len (length of time interval)
        M = 1 - p00(s,i,1);
        if M == 0
            T_len = 0.01;                   % Step one step
            data_coll_S(n) = data_coll_S(n-1);
            data_coll_I(n) = data_coll_I(n-1);
            T_ind(n) = T_ind(n-1) + T_len;
            continue
        else
            T_len = -log(c)/M;
            T_ind(n) = T_ind(n-1) + T_len;      % Next change will happen at this time
        end

        % Determine the change interval
        a1 = p_11(s,i,T_len)/(1 - p00(s,i,T_len));
        a2 = a1 + p0_1(s,i,T_len)/(1 - p00(s,i,T_len));
        a3 = a2 + p1_1(s,i,T_len)/(1 - p00(s,i,T_len));

        c = rand;                           % Pick a random number between (0,1)
        if c < a1                           % (-1,1)
            data_coll_S(n) = data_coll_S(n-1)-1;
            data_coll_I(n) = data_coll_I(n-1)+1;
        elseif c >= a1 && c < a2            % (0,-1)
            data_coll_S(n) = data_coll_S(n-1);
            data_coll_I(n) = data_coll_I(n-1)-1;
        elseif c >= a1 && c < a3            % (1,-1)
            data_coll_S(n) = data_coll_S(n-1)+1;
            data_coll_I(n) = data_coll_I(n-1)-1;
        else                                % (1,0)
            data_coll_S(n) = data_coll_S(n-1)+1;
            data_coll_I(n) = data_coll_I(n-1);
        end
    end

    index_cell{1,k} = T_ind;
    data_cell_S{1,k} = data_coll_S;
    data_cell_I{1,k} = data_coll_I;

end

figure
hold on
for k = 1 : number_simul
    index_coll = index_cell{1,k};
    data_coll_I = data_cell_I{1,k};
    plot(index_coll, data_coll_I)
end
xlim([0,40]);