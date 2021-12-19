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
number_simul = 500;       % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time

%% Define transition probility
p1 = @(i) beta/N*i*(N-i)*dt;          % j = i+1 One plus
p1_ = @(i) (b+gamma)*i*dt;            % j = i-1 One minus
p0 = @(i) 1 - p1(i) - p1_(i);         % j = i Unchanged

%% Make transition matrix
T(1,1) = 1;
for i = 1 : N
    T(i,i+1) = p1_(i);
    T(i+1,i+1) = p0(i);
    T(i+2,i+1) = p1(i);
end
T(N+2,:) = [];

%% Simulation
n = total_time/dt;
data_coll = zeros(n, number_simul);
data_coll(1,:) = I0;

for i = 1 : number_simul
    for j = 2 : n
        if data_coll(j-1,i) == 0            % Infectious individual is zero
            continue                        % Won't change
        elseif data_coll(j-1,i) == N        % Infectious individual is all
            a1 = T(N+1,N+1);                % Only one minus or unchange
            a2 = 1;
        else
            m = data_coll(j-1,i)+1;         % Find next level
            a1 = T(m,m);                    % Less than a1, unchange
            a2 = T(m-1,m)+T(m,m);           % Larger than a2, plus one
        end

        c = rand;           % Pick a random number between (0,1)
        if c < a1           % Less than a1, unnchange
            data_coll(j,i) = data_coll(j-1,i);
        elseif c > a2       % Larger than a2, I plus one
            data_coll(j,i) = data_coll(j-1,i)+1;
        else                % Between a1 and a2, I minus one
            data_coll(j,i) = data_coll(j-1,i) - 1;
        end
    end
end

figure
plot(data_coll)
hold on

final_I = data_coll(n,:);
indx = find(final_I==0);
final_I(indx) = [];

data_coll(:, indx) = [];

average_rel = mean(data_coll,2);
plot(average_rel, 'k', 'LineWidth', 1.5)

figure
histogram(final_I)

%% 2. Set given parameter
N = 100;                % Population
dt = 0.01;              % Delta t, Time Step
beta = 0.25;               % Transmission Rate
b = 0.25;               % Birth Rate
gamma = 0.25;           % Recover Rate
I0 = 20;                 % Initial Infectious Inidvidual

%% Set other parameter
total_time = 50;        % To 50 for each run
number_simul = 500;       % The number of simulations
infuctive_unit = 1;     % Only one person trans for each time

%% Define transition probility
p1 = @(i) beta/N*i*(N-i)*dt;          % j = i+1 One plus
p1_ = @(i) (b+gamma)*i*dt;            % j = i-1 One minus
p0 = @(i) 1 - p1(i) - p1_(i);         % j = i Unchanged

%% Make transition matrix
T(1,1) = 1;
for i = 1 : N
    T(i,i+1) = p1_(i);
    T(i+1,i+1) = p0(i);
    T(i+2,i+1) = p1(i);
end
T(N+2,:) = [];

%% Simulation
n = total_time/dt;
data_coll = zeros(n, number_simul);
data_coll(1,:) = I0;

for i = 1 : number_simul
    for j = 2 : n
        if data_coll(j-1,i) == 0            % Infectious individual is zero
            continue                        % Won't change
        elseif data_coll(j-1,i) == N        % Infectious individual is all
            a1 = T(N+1,N+1);                % Only one minus or unchange
            a2 = 1;
        else
            m = data_coll(j-1,i)+1;         % Find next level
            a1 = T(m,m);                    % Less than a1, unchange
            a2 = T(m-1,m)+T(m,m);           % Larger than a2, plus one
        end

        c = rand;           % Pick a random number between (0,1)
        if c < a1           % Less than a1, unnchange
            data_coll(j,i) = data_coll(j-1,i);
        elseif c > a2       % Larger than a2, I plus one
            data_coll(j,i) = data_coll(j-1,i)+1;
        else                % Between a1 and a2, I minus one
            data_coll(j,i) = data_coll(j-1,i) - 1;
        end
    end
end

figure
plot(data_coll)
hold on

average_rel = mean(data_coll,2);
plot(average_rel, 'k', 'LineWidth', 1.5)