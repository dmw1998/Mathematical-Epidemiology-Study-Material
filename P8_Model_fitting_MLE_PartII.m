% Set up the SEIR model of the transmission dynamics of measles as follows:
% Population 100000 people
% Pre-infectious period 8 days
% Infectious period 7 days
% Initial values (S,E,I,R)=(99999,0,1,0)

N = 100000;						% Population 100,000 people
[S0, E0, I0, R0] = deal(99999, 0, 1, 0);		% Initial values
kappa = 1/8; 					% Pre-infectious period 8 days
alpha = 1/7;					% Infectious period 7 days
