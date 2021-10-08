function sol = solve_SEIR(bet,params,end_time)
f = params.f;
r = params.r;
time = 1:end_time;

sol = ode45(@SEIR_eq, time, params.initial, odeset("NonNegative",1:4));
sol.y = deval(time, sol);
sol.x = time;

    function dydt = SEIR_eq(t,y)
        S = y(1);
        E = y(2);
        I = y(3);
        
        dSdt = -bet*S.*I;
        dEdt = bet*S.*I - f*E;
        dIdt = f*E - r*I;
        dRdt = r*I;
        
        dydt = [dSdt; dEdt; dIdt; dRdt];
    end
end