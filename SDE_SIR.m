clear all; clc; close all

N = 100;                % Population
beta = 1;               % Transmissibility
birth = 0;              % Birth rate
gamma = 0.5;            % Recovery rate
delta_time = 0.01;      % Interevent time
total_time = 50;        % Total time
simulate = 500;         % Number of simulations
I0 = 2;                 % Initial infection
S0 = N - I0;            % Initial susceptible

t = 0:delta_time:total_time;
n = length(t);

p1 = @(s,i) beta*i*s*delta_time/N;
p2 = @(s,i) gamma*i*delta_time;
p3 = @(s,i) birth*i*delta_time;
p4 = @(s,i) birth*(N-s-i)*delta_time;
p5 = @(s,i) 1-(p1(s,i)+p2(s,i)+p3(s,i)+p4(s,i));
E = @(s,i) [(-1)*p1(s,i)+(1)*p3(s,i)+(1)*p4(s,i); (1)*p1(s,i)+(-1)*p2(s,i)+(-1)*p3(s,i)];
V = @(s,i) [(-1)^2*p1(s,i)+1^2*p3(s,i)+1^2*p4(s,i) (-1)*1*p1(s,i)+1*(-1)*p3(s,i); ...
    (-1)*1*p1(s,i)+1*(-1)*p3(s,i) 1^2*p1(s,i)+(-1)^2*p2(s,i)+(-1)^2+p3(s,i)];

simulations = zeros(simulate,n);

for i = 1 : simulate
    D = [S0; I0];
    simulations(i,1) = I0;
    for j = 2 : n
        [b, d] = eig(V(D(1),D(2)));
        sqrtV = b*sqrt(d)*b^(-1);
        normal = [normrnd(0,sqrt(delta_time)); normrnd(0,sqrt(delta_time))];
        D = D + E(D(1),D(2)) + sqrtV*normal;

        if (D(2) < 0 || simulations(i,j-1) == 0)
            D(2) = 0;
        end

        simulations(i,j) = D(2);

    end
end

ave_result = mean(simulations);

figure('WindowState','maximized')
plot(t, simulations)
hold on
plot(t, ave_result, 'k', 'LineWidth', 1.5)
xlabel('Time')
ylabel('Number of Infctives')