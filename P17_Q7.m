close all; clear all; clc
%% Set given parameter
N = 20;                 % Population
dt = 0.01;              % Delta t, Time Step
beta_set = [0.5, 2.5];  % Transmission Rate
b = 0;                  % Birth Rate
gamma = 0.5;            % Recover Rate
I0 = 1;                 % Initial infectioin
S0 = N - I0;

total_time = 50;
number_simul = 50;      % The number of simulations
prob_outbreak = [];
cat = {};

% Define transition probility
p_11 = @(beta,s,i,dt) beta/N*i*s*dt;            % (-1,1) one from S to I
p0_1 = @(beta,s,i,dt) gamma*i*dt;               % (0,-1) one from I to R
p1_1 = @(beta,s,i,dt) b*i*dt;                   % (1,-1) one born in S while one death in I
p10 = @(beta,s,i,dt) b*(N-s-i)*dt;              % (1,0)  one born in S while one death in R
p00 = @(beta,s,i,dt) 1 - p_11(beta,s,i,dt) - p0_1(beta,s,i,dt) - p1_1(beta,s,i,dt) - p10(beta,s,i,dt);          % Unchanged

for beta = beta_set
    cases = 0;

    for k = 1 : number_simul
        i = I0;                    % Initialize the data set
        s = S0;                    % Initialize the data set
        
        % T_ind keep the time of people changed
        T_ind = 0;                              % Initialize the change time
                                                % Start from time = 0
        while T_ind < total_time
            % Pick a random number to decide the length of time interval
            % i.e. to give the next time which the change happends
            c = rand;

            % By given formular, give M (the total change rate) and T_len (length of time interval)
            M = 1 - p00(beta,s,i,1);
            if M == 0
                T_len = 0.01;                   % Step one step
                T_ind = T_ind + T_len;
                continue
            else
                T_len = -log(c)/M;
                T_ind = T_ind + T_len;          % Next change will happen at this time
            end

            % Determine the change interval
            a1 = p_11(beta,s,i,T_len)/(1 - p00(beta,s,i,T_len));
            a2 = a1 + p0_1(beta,s,i,T_len)/(1 - p00(beta,s,i,T_len));
            a3 = a2 + p1_1(beta,s,i,T_len)/(1 - p00(beta,s,i,T_len));

            c = rand;                           % Pick a random number between (0,1)
            if c < a1                           % (-1,1)
                s = s-1;
                i = i+1;
            elseif c >= a1 && c < a2            % (0,-1)
                i = i-1;
            elseif c >= a1 && c < a3            % (1,-1)
                s = s+1;
                i = i-1;
            else                                % (1,0)
                s = s+1;
            end
        end

        if i > 0
            cases = cases + 1;            % Collect the case of one outbreak
        end

    end

    prob_outbreak = [prob_outbreak; cases];
end

figure
bar(prob_outbreak);
xlabel("R_0"); ylabel("probability of an outbreak");
legend(["Simulation", "Theory"],"Location","best");