function optimal_beta = estimate_beta(inc_data, params, initial_guess)
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
    simul = get_weekly_incidence_from_sol(sol,params.f);
    obj = sqrt(sum((inc_points - simul).^2));
end