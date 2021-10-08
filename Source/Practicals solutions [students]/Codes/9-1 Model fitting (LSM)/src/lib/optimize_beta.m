function optimal_beta = optimize_beta(initial_guess, params, inc_data)
time_points = inc_data(:,1);
inc_points = inc_data(:,2);

optimal_beta = fminsearchbnd(...
                @(x) objective(x, params, time_points, inc_points), ...
                initial_guess, ...
                zeros(size(initial_guess)), ...
                []);
end

function obj = objective(bet, params, time_points, inc_points)
    sol = solve_SEIR(bet,params,time_points(end));
    simul = simul_inc(sol,params.f);
    obj = sum((inc_points - simul).^2);
end