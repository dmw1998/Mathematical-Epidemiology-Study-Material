%% PART 1 Setting up difference equations
clear all; clc;
%% fill the %--%
%% 1. Plot a graph for number of S,E,I,R during 200 days.
% parameter setting
N = 100000;%             % Population

pre_inf = 8;%           % pre-infectious duration
f = 1/8;%

inf = 7;%                % infectious duration
gamma = 1/7;%

R_0 = 13;%               % reproduction number
beta = 13/700000;%

% initial
S0 = 99999;%--%
E0 = 0;%--%
I0 = 1;%--%
R0 = 0;%--%

%% 1. run the model for 200 days
tspan =200; % time--%
dt = 1;

% Allocate memories
S = zeros(tspan,1);
E = zeros(tspan,1);%--%
I = zeros(tspan,1);%--%
R = zeros(tspan,1);%--%

% initial 
S(1) = S0; E(1) = E0; I(1) = I0; R(1) = R0;

% Solve difference equation (DFE)
for i = 1:200%--%
    S(i+1) = S(i)-beta*S(i)*I(i);%--%
    E(i+1) = E(i)+beta*S(i)*I(i)-f*E(i);%--%
    I(i+1) = I(i)+f*E(i)-gamma*I(i);%--%
    R(i+1) = R(i)+gamma*I(i)%--%
end

%% Q1. plot
% Define time interval
time_stamp =1:201; %--%;

hold on;
plot(time_stamp,S(:,1));
plot(time_stamp,E(:,1));
plot(time_stamp,I(:,1));
plot(time_stamp,R(:,1));
legend('susceptible','pre-infectious','infectious','recover');
title('SEIR model for 200 days');
grid on; grid minor;
hold off;

saveas(gca, 'Q1', 'epsc')

%% 2-1. How long does it take before there are no infectious persons in the population?
for i = 1:tspan-1
    erd = find(I<1);
end
erd;
fprintf('Eradication occurs on the %d th day\n', erd)

%% 3. How does the graph change if you change the pre-infectious period to be 5 days and 20 days, respectively?
pre_inf = [5,20];
f_val = 1./pre_inf;

% Allocate memories
S = cell(1,%--%);
E = cell(1,%--%);
I = cell(1,%--%);
R = cell(1,%--%);

for j = %--%
    
    % Define time interval
    time_stamp = %--%
    
    % Allocate memories
    S_temp = zeros(%--%, 1);
    E_temp = zeros(%--%, 1);
    I_temp = zeros(%--%, 1);
    R_temp = zeros(%--%, 1);
    
    % Initial states
    S_temp(1) = %--%
    E_temp(1) = %--%
    I_temp(1) = %--%
    R_temp(1) = %--%
    
    % Solve difference equation (DFE)
    for i = 1:dt:tspan+1
        S_temp(i+1) = %--%
        E_temp(i+1) = %--%
        I_temp(i+1) = %--%
        R_temp(i+1) = %--%
    end
    
    S{1,j} = %--%
    E{1,j} = %--%
    I{1,j} = %--%
    R{1,j} = %--%
end

%% Q3. plot

figure('pos', [10 10 1600 900]);

for i = 1:length(f_val)
    subplot(3,1,i) % 3*1 plot
   % define time interval
    time_stamp = 0:dt:tspan+1;
    hold on;
    plot(time_stamp, S{1,i});
    plot(time_stamp, E{1,i});
    plot(time_stamp, I{1,i});
    plot(time_stamp, R{1,i});
    hold off;
    
    xlabel('time (days)')
    ylabel('the number of people')
    xlim([0 200])
    legend('S(susceptible)','E(pre-infectious)','I(infectious)','R(recover)')  
    grid on; grid minor;
    set(gca, 'FontSize', 12)
    title(sprintf('Pre-infection duration = %d days : SEIR model for 200 days', pre_inf(i)))
end
saveas(gca, 'Q3', 'epsc')
