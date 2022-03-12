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

trans_period = 50*365           # Model of 50 years

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
t = np.linspace(0,trans_period,trans_period+1)

# Solve ODE
data_coll = odeint(model,y0,t,args=(beta,kappa,alpha,mu,N))

S = data_coll[:,0]
E = data_coll[:,1]
I = data_coll[:,2]
R = data_coll[:,3]

R_n = R_0*S/N               # Net reproduction number
new_infection = kappa*E     # New infections
S_I = S/I                   # Susceptible to infection
prop_S = S/N                # Proportion of susceptible
prop_R = R/N                # Proportion of immune

t = t/365

# Observe net reproduction number and new infectious
fig1, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time(Years)')
ax1.set_ylabel('Net reproduction number',color=color)
ax1.plot(t[35*365::],R_n[35*365::],color=color)
ax1.axhline(y=1,color='r',linestyle='--')
ax1.tick_params(axis='y',labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('New Infectious',color=color)
ax2.plot(t[35*365::],new_infection[35*365::],color=color)
ax2.tick_params(axis='y',labelcolor=color)

fig1.tight_layout()

# Observe proportion of the population
fig2, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time(Years)')
ax1.set_ylabel('Proportion of population',color=color)
ax1.plot(t[35*365::],prop_R[35*365::],color=color)
ax1.tick_params(axis='y',labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('New Infectious',color=color)
ax2.plot(t[35*365::],new_infection[35*365::],color=color)
ax2.tick_params(axis='y',labelcolor=color)

fig2.tight_layout()

# Observe proportion of R, herd immunity and new infectious
fig3, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time(Years)')
ax1.set_ylabel('Proportion of R',color=color)
ax1.plot(t[35*365::],R_n[35*365::],color=color)
ax1.axhline(y=1-1/R_0,color='r',linestyle='--')
ax1.tick_params(axis='y',labelcolor=color)

ax2 = ax1.twinx()

color = 'tab:blue'
ax2.set_ylabel('New Infectious',color=color)
ax2.plot(t[35*365::],new_infection[35*365::],color=color)
ax2.tick_params(axis='y',labelcolor=color)

fig3.tight_layout()

# Observe 
plt.figure()

plt.plot(t[35*365::],prop_S[35*365::],label='Proportion of S')
plt.plot(t[35*365::],prop_R[35*365::],label='Proportion of R')
plt.axhline(y=1-1/R_0,linestyle='--',label='Herd immunity threshold')
plt.legend()
plt.xlabel('Time(Years)')

plt.show()