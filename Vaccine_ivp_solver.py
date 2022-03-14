from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

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

def ode_sys(t,y):
    N = 100000						# Population 100,000 people
    kappa = 1/8 					# Pre-infectious period 8 days
    alpha = 1/7					    # Infectious period 7 days
    R_0 = 13						# Basic reproduction number 13
    beta = R_0 * alpha / N
    LE = 70 * 365;                  # Life expectancy 70 years
    mu = 1/LE;                      # Birth rate = death rate

    dydt = [-beta*y[0]*y[2] + mu*(N-y[0]),
            beta*y[0]*y[2] - kappa*y[1] - mu*y[1],
            kappa*y[1] - alpha*y[2] - mu*y[2],
            alpha*y[2] - mu*y[3]]
    return dydt

t_begin = 0.
t_end = trans_period_vac
t_nsamples = trans_period_vac+1
t_space = np.linspace(t_begin, t_end, t_nsamples)
y_init = [S0, E0, I0, R0]

method = 'RK45'
result = solve_ivp(ode_sys,[t_begin,t_end],y_init,method=method,dense_output=True)

y_init = result.sol(t_space)[:,trans_period]

# S = result.sol(t_space)[0].T
# E = result.sol(t_space)[1].T
# I = result.sol(t_space)[2].T
# R = result.sol(t_space)[3].T

# R_n = R_0*S/N               # Net reproduction number
# new_infection = kappa*E     # New infections
# S_I = S/I                   # Susceptible to infection
# prop_S = S/N                # Proportion of susceptible
# prop_R = R/N                # Proportion of immune

# t = t_space/365

# # Observe net reproduction number and new infectious
# fig1, ax1 = plt.subplots()

# color = 'tab:red'
# ax1.set_xlabel('Time(Years)')
# ax1.set_ylabel('Net reproduction number',color=color)
# ax1.plot(t[35*365:trans_period],R_n[35*365:trans_period],color=color)
# ax1.axhline(y=1,color='r',linestyle='--')
# ax1.tick_params(axis='y',labelcolor=color)

# ax2 = ax1.twinx()

# color = 'tab:blue'
# ax2.set_ylabel('New Infectious',color=color)
# ax2.plot(t[35*365:trans_period],new_infection[35*365:trans_period],color=color)
# ax2.tick_params(axis='y',labelcolor=color)

# fig1.tight_layout()

# # Observe proportion of the population
# fig2, ax1 = plt.subplots()

# color = 'tab:red'
# ax1.set_xlabel('Time(Years)')
# ax1.set_ylabel('Proportion of population',color=color)
# ax1.plot(t[35*365:trans_period],prop_R[35*365:trans_period],color=color)
# ax1.tick_params(axis='y',labelcolor=color)

# ax2 = ax1.twinx()

# color = 'tab:blue'
# ax2.set_ylabel('New Infectious',color=color)
# ax2.plot(t[35*365:trans_period],new_infection[35*365:trans_period],color=color)
# ax2.tick_params(axis='y',labelcolor=color)

# fig2.tight_layout()

# # Observe proportion of R, herd immunity and new infectious
# fig3, ax1 = plt.subplots()

# color = 'tab:red'
# ax1.set_xlabel('Time(Years)')
# ax1.set_ylabel('Proportion of R',color=color)
# ax1.plot(t[35*365:trans_period],R_n[35*365:trans_period],color=color)
# ax1.axhline(y=1-1/R_0,color='r',linestyle='--')
# ax1.tick_params(axis='y',labelcolor=color)

# ax2 = ax1.twinx()

# color = 'tab:blue'
# ax2.set_ylabel('New Infectious',color=color)
# ax2.plot(t[35*365:trans_period],new_infection[35*365:trans_period],color=color)
# ax2.tick_params(axis='y',labelcolor=color)

# fig3.tight_layout()

# # Observe 
# plt.figure()

# plt.plot(t[35*365:trans_period],prop_S[35*365:trans_period],label='Proportion of S')
# plt.plot(t[35*365:trans_period],prop_R[35*365:trans_period],label='Proportion of R')
# plt.axhline(y=1-1/R_0,linestyle='--',label='Herd immunity threshold')
# plt.legend()
# plt.xlabel('Time(Years)')

I_no_vac = result.sol(t_space)[2].T

plt.figure()
plt.plot(t_space[30*365::]/365, I_no_vac[30*365::], label='No vaccine')

t_begin = trans_period
t_space = np.linspace(t_begin, t_end, t_end-t_begin+1)

# Define ODE function
def model_vac_60(t,y):
    N = 100000						# Population 100,000 people
    kappa = 1/8 					# Pre-infectious period 8 days
    alpha = 1/7					    # Infectious period 7 days
    R_0 = 13						# Basic reproduction number 13
    beta = R_0 * alpha / N
    LE = 70 * 365                   # Life expectancy 70 years
    mu = 1/LE                       # Birth rate = death rate
    c = 0.6

    dydt = [-beta*y[0]*y[2] + (1-c)*mu*(N-y[0]),
            beta*y[0]*y[2] - kappa*y[1] - mu*y[1],
            kappa*y[1] - alpha*y[2] - mu*y[2],
            alpha*y[2] - c*mu*y[3]]
    return dydt

def model_vac_90(t,y):
    N = 100000						# Population 100,000 people
    kappa = 1/8 					# Pre-infectious period 8 days
    alpha = 1/7					    # Infectious period 7 days
    R_0 = 13						# Basic reproduction number 13
    beta = R_0 * alpha / N
    LE = 70 * 365                   # Life expectancy 70 years
    mu = 1/LE                       # Birth rate = death rate
    c = 0.9

    dydt = [-beta*y[0]*y[2] + (1-c)*mu*(N-y[0]),
            beta*y[0]*y[2] - kappa*y[1] - mu*y[1],
            kappa*y[1] - alpha*y[2] - mu*y[2],
            alpha*y[2] - c*mu*y[3]]
    return dydt

result_60 = solve_ivp(model_vac_60,[t_begin,t_end],y_init,method=method,dense_output=True)
result_90 = solve_ivp(model_vac_90,[t_begin,t_end],y_init,method=method,dense_output=True)

I_vac_60 = result_60.sol(t_space)[2].T
I_vac_90 = result_90.sol(t_space)[2].T

plt.plot(t_space/365, I_vac_60, label='60% vaccine')
plt.plot(t_space/365, I_vac_90, label='90% vaccine')
plt.legend()
plt.show()