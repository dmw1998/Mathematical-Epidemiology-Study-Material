%% PART 2 Incorporating births and deaths
clear all; clc;

%% parameter setting
N = 1e5;
LE = 70 * 365; % life expactancy

pre_inf = 8;
f = 1/pre_inf;

inf = 7;
gamma = 1/inf;

R_0 = 13;
beta = R_0 / (N * inf);

% initial
S0 = N - 1;
E0 = 0 ;
I0 = 1 ;
R0 = 0 ;

%% 7. How does this change your answer to question Q2?
tspan = 18250;
dt = 1;

% Allocate memories
S = zeros(tspan,1);
E = zeros(size(S));
I = zeros(size(S));
R = zeros(size(S));

% initial 
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;
b = 1/LE; d = b;

% Solve difference equation (DFE)
for i = 1:tspan-1
    S(i+1) = S(i) + (b * N -beta * S(i) * I(i) - d * S(i)) * dt;
    E(i+1) = E(i) + (beta * S(i) * I(i) - f * E(i) - d * E(i)) * dt;
    I(i+1) = I(i) + (f * E(i) - gamma * I(i) - d * I(i)) * dt;
    R(i+1) = R(i) + (gamma * I(i) - d * R(i)) * dt;
end

%% Q7. plot
time_stamp = 0:dt:tspan-1;

hold on
plot(time_stamp,S(:,1));
plot(time_stamp,E(:,1));
plot(time_stamp,I(:,1));
plot(time_stamp,R(:,1));

xlabel('time (days)');
ylabel('the number of people');
legend('S(susceptible)','E(pre-infectious)','I(infectious)','R(recover)');
title('SEIR model for 200 days')
set(gca, 'FontSize', 12)
hold off
saveas(gca, 'Q7', 'epsc')

%% find eradication date
for i = 1:tspan-1
    erd = find(I<1);
end
erd;
fprintf('Eradication occurs on the %d th day\n', erd(3))
%% 8. How are the changes in the number of susceptible and immune related if you were to simulate the dynamics of measles for 50 years and for 100 years?

tspan_val = [50*365, 100*365];
dt = 1;

% Allocate memories
S = cell(1,length(tspan_val));
E = cell(1,length(tspan_val));
I = cell(1,length(tspan_val));
R = cell(1,length(tspan_val));

for j = 1:length(tspan_val)
    
    % Define time interval
    time_stamp = 1:dt:tspan_val(j)+1;
    
    % Allocate memories
    S_temp = zeros(length(time_stamp), 1);
    E_temp = zeros(length(time_stamp), 1);
    I_temp = zeros(length(time_stamp), 1);
    R_temp = zeros(length(time_stamp), 1);
    
    % Initial states
    S_temp(1) = S0;
    E_temp(1) = E0;
    I_temp(1) = I0;
    R_temp(1) = R0;
    
    % Solve difference equation (DFE)
    for i = 1:dt:tspan_val(j)+1
        S_temp(i+1) = S_temp(i) + (- beta .* I_temp(i) .* S_temp(i) + b .* N - d.* S_temp(i)) * dt;
        E_temp(i+1) = E_temp(i) + (beta .* I_temp(i) .* S_temp(i) - f .* E_temp(i) - d .* E_temp(i)) * dt;
        I_temp(i+1) = I_temp(i) + (f .* E_temp(i) - gamma .* I_temp(i) - d .* I_temp(i)) * dt;
        R_temp(i+1) = R_temp(i) + (gamma .* I_temp(i) - d .* R_temp(i)) * dt;
    end
    
    S{1,j} = S_temp;
    E{1,j} = E_temp;
    I{1,j} = I_temp;
    R{1,j} = R_temp;
end

%% Q8. plot

figure('pos', [10 10 1600 900]);

for i = 1:length(tspan_val)
    subplot(2,1,i) % 2*2 plot
    time_stamp = 0:dt:tspan_val(i)+1;
    hold on;
    plot(time_stamp./365, S{1,i});
    plot(time_stamp./365, E{1,i});
    plot(time_stamp./365, I{1,i});
    plot(time_stamp./365, R{1,i});
    hold off;
    
    xlabel('time (years)')
    ylabel('the number of people')
    xlim([0 100])
    legend('S(susceptible)','E(pre-infectious)','I(infectious)','R(recover)')  
    grid on; grid minor;
    set(gca, 'FontSize', 12)
    title(sprintf('SEIR model for %d years', tspan_val(i)/365))
end
saveas(gca, 'Q8', 'epsc')


%% 9-1. What happens if you change the time step to 2, 3, 4, and 5 days?
dt = [2,3,4,5];

S = cell(1,length(dt));
E = cell(1,length(dt));
I = cell(1,length(dt));
R = cell(1,length(dt));

for j = 1:length(dt)
    % Define time interval
    time_stamp = 0:dt(j):tspan;
    
    % Memory allocation
    S_temp = zeros(length(time_stamp), 1);
    E_temp = zeros(length(time_stamp), 1);
    I_temp = zeros(length(time_stamp), 1);
    R_temp = zeros(length(time_stamp), 1);
    
    % Initial states
    S_temp(1) = S0;
    E_temp(1) = E0;
    I_temp(1) = I0;
    R_temp(1) = R0;
    
    % Adjust unit time to parameters
    f_temp = f .* dt(j);
    gamma_temp = gamma .* dt(j);
    beta_temp = beta .* dt(j);
    d_temp = d .* dt(j);
    b_temp = b .* dt(j);
    for i = 1:length(time_stamp)-1
        S_temp(i+1) = S_temp(i) - beta_temp .* I_temp(i) .* S_temp(i) + b_temp .* N - d_temp.* S_temp(i);
        E_temp(i+1) = E_temp(i) + beta_temp .* I_temp(i) .* S_temp(i) - f_temp .* E_temp(i) - d_temp .* E_temp(i);
        I_temp(i+1) = I_temp(i) + f_temp .* E_temp(i) - gamma_temp .* I_temp(i) - d_temp .* I_temp(i);
        R_temp(i+1) = R_temp(i) + gamma_temp .* I_temp(i) - d_temp .* R_temp(i);
    end
    S{1,j} = S_temp;
    E{1,j} = E_temp;
    I{1,j} = I_temp;
    R{1,j} = R_temp;
    
end

%% Q8. Plot
figure('pos', [10 10 1600 900]);

for i = 1:length(dt)
    subplot(2,2,i) % 2*2 plot
    time_stamp = 1:dt(i):tspan+1;
    hold on;
    plot(time_stamp, S{1,i});
    plot(time_stamp, E{1,i});
    plot(time_stamp, I{1,i});
    plot(time_stamp, R{1,i});
    hold off;
    
    xlabel('time (days)')
    ylabel('the number of people')
    xlim([0 18250])
    legend('S(susceptible)','E(pre_infectious)','I(infectious)','R(recover)')  
    grid on; grid minor;
    set(gca, 'FontSize', 12)
    title(sprintf('unit time = %d days :SEIR model for 200 days', dt(i)))
end
saveas(gca, 'Q9.eps', 'epsc')





