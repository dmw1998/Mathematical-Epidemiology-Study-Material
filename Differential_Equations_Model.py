import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 100000						# Population 100,000 people
[S0, E0, I0, R0] = [99999, 0, 1, 0]		# Initial values
kappa = 1/8 					# Pre-infectious period 8 days
alpha = 1/7					# Infectious period 7 days
R_0 = 13						# Basic reproduction number 13
beta = R_0 * alpha / N
trans_period = 150              # Model of 150 days

# Define ODE function
def model(y,t,beta,kappa,alpha):
    dydt = [-beta*y[0]*y[2],
            beta*y[0]*y[2] - kappa*y[1],
            kappa*y[1] - alpha*y[2],
            alpha*y[2]]
    return dydt

# Initial guess
y0 = [S0, E0, I0, R0]

# Time points
t = np.linspace(0,trans_period,301)

# Solve ODE
data_coll = odeint(model,y0,t,args=(beta,kappa,alpha))

plt.figure()
plt.plot(t, data_coll[:,0], label='susceptible')
plt.plot(t, data_coll[:,1], label='infectious')
plt.plot(t, data_coll[:,2], label='pre-infectious')
plt.plot(t, data_coll[:,3], label='recovered')
plt.legend()
plt.title('SEIR model for %d days pre-inf %d days' %(trans_period,1/kappa))
plt.xlabel('Time(Days)')
plt.ylabel('Population')

## Consider natural death
LE = 70*365         # Life expectancy
mu = 1/LE           # Birth rate = death rate

trans_period = 100*365              # Model of 100 years

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

t = t/365

plt.figure()
plt.plot(t, data_coll[:,0], label='susceptible')
plt.plot(t, data_coll[:,1], label='infectious')
plt.plot(t, data_coll[:,2], label='pre-infectious')
plt.plot(t, data_coll[:,3], label='recovered')
plt.legend()
plt.title('SEIR model for %d years pre-inf %d days' %(trans_period/365,1/kappa))
plt.xlabel('Time(Years)')
plt.ylabel('Population')

plt.show()