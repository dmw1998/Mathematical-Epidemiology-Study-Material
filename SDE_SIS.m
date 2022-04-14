clear all; clc; close all

N = 100;                % Population
beta = 0.25;            % Transmissibility
birth = 0.25;           % Birth rate
gamma = 0.25;           % Recovery rate
delta_time = 0.01;      % Interevent time
total_time = 50;        % Total time
simulate = 500;         % Number of simulations
I0 = 20;                % Initial infection

t = 0:delta_time:total_time;
n = length(t);

b = @(i,dt) beta*i*(N-i)*dt/N;
d = @(i,dt) (birth+gamma)*i*dt;

simulations = zeros(simulate,n);

for i = 1 : simulate
    simulations(i,1) = I0;
    for j = 2 : n 
        I = simulations(i,j-1);
        E = b(I,1) - d(I,1);
        V = b(I,1) + d(I,1);
        simulations(i,j) = I + E*delta_time + sqrt(V)*sqrt(delta_time)*normrnd(0,1);

        if simulations(i,j) < 0
            simulations(i,j) = 0;
        end
    end
end

ave_result = mean(simulations);

figure('WindowState','maximized')
plot(t, simulations)
hold on
plot(t, ave_result, 'k', 'LineWidth', 1.5)
xlabel('Time')
ylabel('Number of Infctives')