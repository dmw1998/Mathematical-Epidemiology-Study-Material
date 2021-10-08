function params = set_params()
params.N = 6000;
params.f = 1/8; % /days
params.r = 1/7; % /days
params.initial = [params.N-1; 0; 1; 0];
end