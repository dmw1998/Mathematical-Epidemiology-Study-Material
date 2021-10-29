function weekly_inc = solve_SEIR(theta,param)

para = [theta;param.f;param.r];
time_stampc = param.time_stampc;

sol = ode45(@(t,y) ode_SEIR(t,y,para),time_stampc,param.initial);
E = deval(time_stampc,sol,2);
daily_inc = param.f*(E(1:end-1)+E(2:end))/2;
for i = 1 : 20
    weekly_inc(i) = sum(daily_inc(7*(i-1)+1:7*i));
end
weekly_inc = weekly_inc';


