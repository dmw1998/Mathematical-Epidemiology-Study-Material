from cProfile import label
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 100000						# Population 100,000 people
[S0, E0, I0, R0] = [99999, 0, 1, 0]		# Initial values
kappa = 1/8 					# Pre-infectious period 8 days
alpha = 1/7					    # Infectious period 7 days
R_0 = 13						# Basic reproduction number 13
beta = R_0 * alpha / N
LE = 70 * 365;                  # Life expectancy 70 years
mu = 1/LE;                      # Birth rate = death rate

c1 = 0.6                        # Vaccine 60% coverage
c2 = 0.9                        # Vaccine 90% coverage

trans_period = 50*365           # No vaccine applied for first 50 years
trans_period_vac = 80*365       # Obsever another 30 years with vaccine

# Define ODE function
def model(y,t,beta,kappa,alpha,mu,N):
    dydt = [-beta*y[0]*y[2] + mu*(N-y[0]),
            beta*y[0]*y[2] - kappa*y[1] - mu*y[1],
            kappa*y[1] - alpha*y[2] - mu*y[2],
            alpha*y[2] - mu*y[3]]
    return dydt

# Initial guess
y0 = [S0, E0, I0, R0]

# Time points
t1 = np.linspace(0,trans_period_vac,trans_period_vac+1)

# Solve ODE
data_coll = odeint(model,y0,t1,args=(beta,kappa,alpha,mu,N))

I = data_coll[:,2]

plt.figure()

plt.plot(t1/365,I,label='No vaccine')

# Define ODE function
def model_vac(y,t,beta,kappa,alpha,mu,N,c):
    dydt = [-beta*y[0]*y[2] + (1-c)*mu*(N-y[0]),
            beta*y[0]*y[2] - kappa*y[1] - mu*y[1],
            kappa*y[1] - alpha*y[2] - mu*y[2],
            alpha*y[2] - c*mu*y[3]]
    return dydt

# Initial guess
y0 = data_coll[trans_period,:]

# Time points
t2 = np.linspace(trans_period+1,trans_period_vac,trans_period_vac-trans_period+1)

# Solve ODE
data_coll1 = odeint(model_vac,y0,t1,args=(beta,kappa,alpha,mu,N,c1))

plt.plot(t1/365,data_coll1[:,2],label='60 coverage')

# Solve ODE
data_coll2 = odeint(model_vac,y0,t1,args=(beta,kappa,alpha,mu,N,c2))

plt.plot(t1/365,data_coll2[:,2],label='90 coverage')

plt.legend()
plt.xlabel('Time(Years)')
plt.xlim((30,80))
plt.ylim((0,120))

plt.show()

"""
ivp_solver: RK45 https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
odeint: LSODA https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html#scipy.integrate.odeint

lsoda: https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html
"""