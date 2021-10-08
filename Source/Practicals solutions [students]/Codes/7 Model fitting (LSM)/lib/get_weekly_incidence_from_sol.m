function simul = get_weekly_incidence_from_sol(sol,f)
dt = sol.x(2)-sol.x(1);
    
E = sol.y(2,:);
daily_inc = f*(E(1:end-1) + E(2:end))/2*dt;

simul = zeros(sol.x(end)/7,1);
for i = 1:length(simul)-1
    for j = 1:7
        simul(i) = simul(i) + daily_inc(7*(i-1)+j);
    end
end
end