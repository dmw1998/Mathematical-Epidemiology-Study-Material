close all; clear all; clc

%% 1. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 1;               % Transmission Rate
b = 0.25;               % Birth Rate
gamma = 0.25;           % Recover Rate
I0 = 2;                 % Initial Infectious Inidvidual

%% Set other parameter
total_time = 50;        % To 50 for each run
number_simul = 10;     % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time

%% Define transition probility
p1 = @(i,dt) beta/N*i*(N-i)*dt;          % j = i+1 One plus
p1_ = @(i,dt) (b+gamma)*i*dt;            % j = i-1 One minus
p0 = @(i,dt) 1 - p1(i,dt) - p1_(i,dt);   % j = i Unchanged

%% Simulation
data_cell = cell(1, number_simul);
index_cell = cell(1,number_simul);

for k = 1 : number_simul
    data_coll = [];
    data_coll(1) = I0;                    % Initialize the data set
    
    % T_ind keep the time of people changed
    n = 1;
    T_ind = [];                             % Initialize the change time
    T_ind(1) = 0;                           % Start from time = 0
    while T_ind(n) < total_time
        i = data_coll(n);
        % Pick a random number to decide the length of time interval
        % i.e. to give the next time which the change happends
        c = rand;

        n = n+1;                            % Next event
        % By given formular, give M (the total change rate) and T_len (length of time interval)
        if i == 0
            T_len = 0.01;           % Step one step
            data_coll(n) = data_coll(n-1);
            T_ind(n) = T_ind(n-1) + T_len;
            continue
        else
            M = 1 - p0(i,1);             % M = p1 + p1_ (the probability of changes)
            T_len = -log(c)/M;
        end

        T_ind(n) = T_ind(n-1) + T_len;      % Next change will happen at this time

        % Determine the change interval
        a1 = p1(i,T_len)/(p1(i,T_len)+p1_(i,T_len));

        c = rand;           % Pick a random number between (0,1)
        if c < a1           % Less than a1, I plus one
            data_coll(n) = data_coll(n-1)+1;
        else                % Larger than a2, I minus one
            data_coll(n) = data_coll(n-1)-1;
        end
    end

    index_cell{1,k} = T_ind;
    data_cell{1,k} = data_coll;

end

figure
hold on
for k = 1 : number_simul
    index_coll = index_cell{1,k};
    data_coll = data_cell{1,k};

    if data_coll(length(data_coll)) == 0
        continue
    else
        plot(index_coll, data_coll)
    end
end

%% 2. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 0.25;               % Transmission Rate
b = 0.25;               % Birth Rate
gamma = 0.25;           % Recover Rate
I0 = 20;                 % Initial Infectious Inidvidual

%% Simulation
data_cell = {}; index_cell = {};
data_cell = cell(1, number_simul);
index_cell = cell(1,number_simul);

for k = 1 : number_simul
    data_coll = [];
    data_coll(1) = I0;                    % Initialize the data set
    
    % T_ind keep the time of people changed
    n = 1;
    T_ind = [];                             % Initialize the change time
    T_ind(1) = 0;                           % Start from time = 0
    while T_ind(n) < total_time
        i = data_coll(n);
        % Pick a random number to decide the length of time interval
        % i.e. to give the next time which the change happends
        c = rand;

        n = n+1;                            % Next event
        % By given formular, give M (the total change rate) and T_len (length of time interval)
        if i == 0
            T_len = dt;           % Step one step
            data_coll(n) = data_coll(n-1);
            T_ind(n) = T_ind(n-1) + T_len;
            continue
        else
            M = 1 - p0(i,1);             % M = p1 + p1_ (the probability of changes)
            T_len = -log(c)/M;
        end

        T_ind(n) = T_ind(n-1) + T_len;      % Next change will happen at this time

        % Determine the change interval
        a1 = p1(i,T_len)/(p1(i,T_len)+p1_(i,T_len));

        c = rand;           % Pick a random number between (0,1)
        if c < a1           % Less than a1, I plus one
            data_coll(n) = data_coll(n-1)+1;
        else                % Larger than a2, I minus one
            data_coll(n) = data_coll(n-1)-1;
        end
    end

    index_cell{1,k} = T_ind;
    data_cell{1,k} = data_coll;

end

figure
hold on
for k = 1 : number_simul
    index_coll = index_cell{1,k};
    data_coll = data_cell{1,k};
    plot(index_coll, data_coll)
end