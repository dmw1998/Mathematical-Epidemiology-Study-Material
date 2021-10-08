clear all ; clc

%% R0 change
% Parameter
N = 100000;
preinf = 8;% days
D = 7;% days
r0_val = [13 5 18];
lifeexp = 70; %years

% parameter set
beta = zeros(1, length(r0_val));
h = zeros(1, length(r0_val));
for i = 1:length(r0_val)
beta(i) = r0_val(i)/(N*D);
h(i) = 1-1/r0_val(i);
end
f = 1/preinf ;
r = 1/D ;
d = 1/(lifeexp*365);
b = d;


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

S = zeros(length(tspan),length(r0_val));
E = zeros(size(S));
I = zeros(size(S));
R = zeros(size(S));
propsus = zeros(size(S));
propimm = zeros(size(S));
Rn = zeros(size(S));

% ODE
for i = 1:length(r0_val)
fode = @(t,y) [-beta(i)*y(1)*y(3) + b*N - d*y(1);
                beta(i)*y(1)*y(3) - f*y(2) - d*y(2);
                f*y(2) - r*y(3) - d*y(3);
                r*y(3) - d*y(4) ] ;
            
[~,sol] = ode45 (fode, tspan, [S0; E0; I0; R0]); 

S(:,i) = sol(:,1);
E(:,i) = sol(:,2);
I(:,i) = sol(:,3);
R(:,i) = sol(:,4);

propsus(:,i) = S(:,i) / N;
propimm(:,i) = R(:,i) / N;
Rn(:,i) = r0_val(i) * propsus(:,i);
end

%% Get new incidence
new_preinf = zeros(length(tspan)-1, length(r0_val));
new_inf = zeros(size(new_preinf));
new_recov = zeros(size(new_preinf));

for i = 1:length(r0_val)
new_preinf(:,i) = beta(i) * (S(1:end-1,i) .* I(1:end-1,i) + S(2:end,i) .* I(2:end,i)) / 2;
new_inf(:,i) = f * (E(1:end-1,i) + E(2:end,i)) / 2;
new_recov(:,i) = r * (I(1:end-1,i) + I(2:end,i)) / 2;
end
%% Result Q1 (50~100) R0=13 vs. 5 vs. 18

hold on;
plot (tspan(2:end)/365,new_inf(:,1));
plot (tspan(2:end)/365,new_inf(:,2));
plot (tspan(2:end)/365,new_inf(:,3));
xlim ([50, 100]);
xlabel ('Time(years)');
ylabel ('The number of new infectious')
hold off;

title ('the effect of the R0 on the cycles in the daily number of new infectious');
legend ('R0 = 13', 'R0 = 5', 'R0 = 18')
grid on; grid minor;
 
%% life exp change 
clear all; clc;

N = 100000;
preinf = 8;% days
D = 7;% days-inf
r0 = 13;
lifeexp_val = [50 70 90]; %years

% parameter set
beta = r0/(N*D);
f = 1/preinf ;
r = 1/D ;
d = zeros(length(lifeexp_val),1);
b = zeros(length(lifeexp_val),1);
for i = 1:length(lifeexp_val)
d(i) = 1/(lifeexp_val(i)*365);
b(i) = d(i);
end
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


% ODE
for i = 1:length(lifeexp_val)
fode = @(t,y) [-beta*y(1)*y(3) + b(i)*N - d(i)*y(1);
                beta*y(1)*y(3) - f*y(2) - d(i)*y(2);
                f*y(2) - r*y(3) - d(i)*y(3);
                r*y(3) - d(i)*y(4) ] ;
            
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
new_preinf = zeros(length(tspan)-1, length(lifeexp_val));
new_inf = zeros(size(new_preinf));
new_recov = zeros(size(new_preinf));

for i = 1:length(lifeexp_val)
new_preinf(:,i) = beta * (S(1:end-1,i) .* I(1:end-1,i) + S(2:end,i) .* I(2:end,i)) / 2;
new_inf(:,i) = f * (E(1:end-1,i) + E(2:end,i)) / 2;
new_recov(:,i) = r * (I(1:end-1,i) + I(2:end,i)) / 2;
end

%% plot life 50 vs 70 vs 90
hold on;
plot (tspan(2:end)/365,new_inf(:,1),'-.');
plot (tspan(2:end)/365,new_inf(:,2));
plot (tspan(2:end)/365,new_inf(:,3));
xlim ([0, 150]);
ylim ([0, 15]);
xlabel ('Time(years)');
ylabel ('The number of new infectious')
hold off ;
legend ('birth rate(1/70*365)', 'birth rate(1/50*365)', 'birth rate(1/90*365)');
title ('the effect of the birth rate on the inter-epidemic period')

%% Influenza pre =2, D =2, RO=2
clear all ; clc

N = 100000;
preinf = 2;% days
D = 2;% days-inf
r0 = 2;
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
dt = 0.01;% 1 day
t_start = 0;
t_end = 1365*100;
tspan = t_start:dt:t_end;
nt = length(tspan);


% ODE
fode = @(t,y) [-beta*y(1)*y(3) + b*N - d*y(1);
                beta*y(1)*y(3) - f*y(2) - d*y(2);
                f*y(2) - r*y(3) - d*y(3);
                r*y(3) - d*y(4) ] ;
            
[~,sol] = ode45 (fode, tspan, [S0; E0; I0; R0]); 

S = sol(:,1);
E = sol(:,2);
I = sol(:,3);
R = sol(:,4);

propsus = sol(:,1) / N;
propimm = sol(:,4) / N;
Rn = r0 * propsus;


% Get new incidence
new_preinf = beta * (S(1:end-1) .* I(1:end-1) + S(2:end) .* I(2:end)) / 2;
new_inf = f * (E(1:end-1) + E(2:end)) / 2;
new_recov = r * (I(1:end-1) + I(2:end)) / 2;

%% plot
plot (tspan(2:end),new_inf);
xlim ([18250, 36500]);
xlabel ('day');
hold on ;

legend ('new-inf');
title ('inflienza')