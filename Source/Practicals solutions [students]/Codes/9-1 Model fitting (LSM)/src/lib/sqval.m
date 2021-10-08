function sq = sqval(sero_prev, simulation)
%negLL - calculate square root of squared sum 
%
% Syntax: sq = sqval(sero_prev, simulation)
%
% Input: seroprevelence data and simulation. They must have the same size of elements.
% Output: Least square value. 

%% Check size error
str_data_size = sprintf("%g ",size(sero_prev));
str_simul_size = sprintf("%g ",size(simulation));
errMsg = "sqval: data has the size of ["+str_data_size ...
    +"], but simulation result has the size of ["+str_simul_size+"]";
check_simul_size = prod(size(sero_prev) == size(simulation));
assert(logical(check_simul_size), errMsg);

%% Calculate the objective function

sq = sum((sero_prev - simulation).^2);
end