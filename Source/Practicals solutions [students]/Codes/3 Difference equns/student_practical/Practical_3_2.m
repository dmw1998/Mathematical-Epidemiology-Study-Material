%% PART 2 Incorporating births and deaths
clear all; clc;

%% parameter setting
N = 1e5;
LE = 70*365;%--% % life expactancy

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
tspan = 200;
dt = 1;

% Allocate memories
S = zeros(tspan,1);
E = zeros(size(S));
I = zeros(size(S));
R = zeros(size(S));

% initial 
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;
b = %--%; 
d = %--%; % birth and death

% Solve difference equation (DFE)
for i = %--%
    S(i+1) = %--%
    E(i+1) = %--%
    I(i+1) = %--%
    R(i+1) = %--%
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
grid on; grid minor;
hold off

saveas(gca, 'Q7', 'epsc')

%% Q7. find eradication date
for i = 1:tspan-1
    erd = find(I<1);
end
erd;
fprintf('Eradication occurs on the %d th day\n', erd(3))

%% 8. How are the changes in the number of susceptible and immune related if you were to simulate the dynamics of measles for 50 years and for 100 years?

tspan_val = [%--%];
dt = 1;

% Allocate memories
S = cell(1,%--%);
E = cell(1,%--%);
I = cell(1,%--%);
R = cell(1,%--%);

for j = %--%
    
    % Define time interval
    time_stamp = %--%
    
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
    for i = %--%
        S_temp(i+1) = %--%
        E_temp(i+1) = %--%
        I_temp(i+1) = %--%
        R_temp(i+1) = %--%
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
    plot(time_stamp./365, S{1,i}, 'LineWidth', 2);
    plot(time_stamp./365, E{1,i}, 'LineWidth', 2);
    plot(time_stamp./365, I{1,i}, 'LineWidth', 2);
    plot(time_stamp./365, R{1,i}, 'LineWidth', 2);
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

%--% make your self!

%% Q9. Plot

%--% make your self!
