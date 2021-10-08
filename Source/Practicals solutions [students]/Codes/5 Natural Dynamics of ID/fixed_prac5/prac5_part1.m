clear all ; clc

%% Parameter
N = 100000;
preinf = 8;% days
D = 7;% days
r0 = 13;
lifeexp = 70; %years

% parameter set
beta = r0/(N*D);
f = 1/preinf ;
r = 1/D ;
d = 1/(lifeexp*365);
b = d;
h = 1-1/r0;

% Initial values
S0 = 99999;
E0 = 0;
I0 = N - S0;
R0 = 0;

% Time
dt = 1;% 1 day
t_start = 0;
t_end = 365*150;
tspan = t_start:dt:t_end;
nt = length(tspan);

%% Differentail equation (1)

% ODE
fode = @(t,y) [-beta*y(1)*y(3) + b*N - d*y(1);
                beta*y(1)*y(3) - f*y(2) - d*y(2);
                f*y(2) - r*y(3) - d*y(3);
                r*y(3) - d*y(4) ] ;
            
[~,sol] = ode45 (fode, tspan, [S0; E0; I0; R0]); 

S = sol(:,1);   % sol 1열
E = sol(:,2);   % sol 2열
I = sol(:,3);   % sol 3열
R = sol(:,4);   % sol 4열

propsus = sol(:,1) / N;   % S/N
propimm = sol(:,4) / N;   % R/N
Rn = r0 * propsus;        % Rn = R0 * S/N


% Get new incidence
new_preinf = beta * (S(1:end-1) .* I(1:end-1) + S(2:end) .* I(2:end)) / 2;
new_inf = f * (E(1:end-1) + E(2:end)) / 2;
new_recov = r * (I(1:end-1) + I(2:end)) / 2;

                
%% Result Q1

figure(1)
yyaxis left % y축 왼쪽 
plot (tspan,Rn);
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('Net reproduct number')
ylim ([0,1.5]);
hold on;

yyaxis right % y축 오른쪽
plot (tspan(2:end),new_inf);
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('Number of new infectious/day');
ylim ([0,8]);
hold on ;
grid on; grid minor;

title ('Comparison between Rn and the daily number of new infectious');
legend ('Rn','the number of new infectious');




 
%% RESULT Q2

yyaxis left % y축 왼쪽 
plot (tspan,propsus,'-b');
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylim ([0,0.15]);
ylabel('Proportion of susceptible')
hold on;


yyaxis right % y축 오른쪽
plot (tspan(2:end),new_inf);
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylim ([0,8]);
title ('Comparison between the proportion of susceptible and the daily number of new infectious');
legend('proportion of susceptible','New infectious');


grid on; grid minor;

%% RESULT  Q3

yyaxis left % y축 왼쪽 
plot (tspan,propimm,'-b');
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('Proportion of immune')
ylim ([0.8,1.0]);
hold on;


yyaxis right % y축 오른쪽
plot (tspan(2:end),new_inf);
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('The number of new infectious/day')
ylim ([0,8]);
hold on ;

title ('Comparison between the proportion immune and the daily number of new infectious');
legend ('Proportion of immnune', 'the number of New infectious');

grid on; grid minor;


%% all plot

yyaxis left % y축 왼쪽
hold on
plot (tspan,propsus, '-r');
plot (tspan,propimm, '-c');
plot (tspan,Rn, '-b');
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('Proportion of susceptible or immune or Rn')
ylim ([0,1.5]);
hold off;


yyaxis right % y축 오른쪽
plot (tspan(2:end),new_inf);
xlim ([14650, 18250]);
xlabel ('Time(days)');
ylabel ('The number of new infectious/day')
ylim ([0,8]);

title ('Comparison between the proportion or Rn and the daily number of new infectious');
legend ('Proportion of susceptible', 'Proportion of immune', 'Rn', 'the number of New infectious');

grid on; grid minor;
%% Parameter
clear all; clc;

N = 100000;
preinf = 8;% days
D = 7;% days
r0 = 13;
lifeexp = 70; %years

% parameter set
beta = r0/(N*D);
f = 1/preinf ;
r = 1/D ;
d = 1/(lifeexp*365);
b = d;
h = 1-1/r0;

% Initial values
S0 = 99999;
E0 = 0;
I0 = N - S0;
R0 = 0;

dt = 1;% 1 day
t_start = 0;
t_end = 365*150;
tspan = t_start:dt:t_end;
nt = length(tspan);

%% ODE 2 insert vac60
vacc_val = [0 60 75 90 95];

S = zeros(length(tspan), length(vacc_val));
E = zeros(length(tspan), length(vacc_val));
I = zeros(length(tspan), length(vacc_val));
R = zeros(length(tspan), length(vacc_val));
propsus = zeros(size(S));
propimm = zeros(size(S));
Rn = zeros(size(S));
new_preinf = zeros(length(tspan)-1, length(vacc_val));
new_inf = zeros(size(new_preinf));
new_recov = zeros(size(new_preinf));

%ODE
for i = 1:length(vacc_val)
fode = @(t,y) [-beta*y(1)*y(3) + (1-vac(vacc_val(i),t))*b*N - d*y(1);
                beta*y(1)*y(3) - f*y(2) - d*y(2);
                f*y(2) - r*y(3) - d*y(3);
                r*y(3) + vac(vacc_val(i),t)*b*N - d*y(4) ] ;
            
[~,sol] = ode45 (fode, tspan, [S0; E0; I0; R0]); 

S(:,i) = sol(:,1);
E(:,i) = sol(:,2);
I(:,i) = sol(:,3);
R(:,i) = sol(:,4);

propsus(:,i) = S(:,i) / N;
propimm(:,i) = R(:,i) / N;
Rn(:,i) = r0 * propsus(:,i);

end
%% Get new incidence
for i = 1:length(vacc_val)
new_preinf(:,i) = beta * (S(1:end-1,i) .* I(1:end-1,i) + S(2:end,i) .* I(2:end,i)) / 2;
new_inf(:,i) = f * (E(1:end-1,i) + E(2:end,i)) / 2;
new_recov(:,i) = r * (I(1:end-1,i) + I(2:end,i)) / 2;
end

%% Result Q4
hold on;
plot (tspan(2:end)/365,new_inf(:,1),'-.');
plot (tspan(2:end)/365,new_inf(:,2));
plot (tspan(2:end)/365,new_inf(:,4));
xlim ([10950/365, 150]);
xlabel ('Time(years)');
ylabel ('the number of new infectious')
title ('The long-term daily number of new infectious')
legend ('vac. coverage 0%', 'vac.coverage 60%', 'vac.coverage 90%')

grid on; grid minor;


%% Result Q5

hold on;
plot (tspan(2:end)/365,new_inf(:,1),'-.');
plot (tspan(2:end)/365,new_inf(:,5));
xlim ([10950/365, 150]);
xlabel ('Time(years)');
ylabel ('the number of new infectious')
title ('The long-term daily number of new infectious')
legend ('vac. coverage 0%', 'vac.coverage 95%')

grid on; grid minor;