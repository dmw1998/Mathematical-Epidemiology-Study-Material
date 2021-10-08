function dy = ode_SEIR(t,y,para)
    
S = y(1);
E = y(2);
I = y(3);

beta = para(1);
f = para(2);
r = para(3);

dS = -beta*S.*I;
dE = beta*S.*I -f*E;
dI = f*E -r*I;
dR = r*I;

dy = [dS;dE;dI;dR];
