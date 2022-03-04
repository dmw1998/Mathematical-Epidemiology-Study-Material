from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

N = 100000          # Population 100.000 people
[S0, E0, I0, R0] = [99999, 0, 1, 0]     # Initial values
kappa = 1/8         # Pre-infectious period 8 days
alpha = 1/7         # Infectious period 7 days
R_0 = 13            # Basic reproduction number 13
beta = R_0*alpha/N  # Contact rate

trans_period = 200  # Simulate the firt 200 days
dt = 1              # Time step 1 day

data_coll = []
I = [I0]
data_coll.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll[i-1][2], 0, 0, 0],
                        [beta*data_coll[i-1][2], -kappa, 0, 0],
                        [0, kappa, -alpha, 0],
                        [0, 0, alpha, 0]])
    data_coll.append(data_coll[i-1] + dt*np.matmul(diff_mat,data_coll[i-1]))
    I.append(data_coll[i][2])

days = range(0,trans_period+1)
plt.figure()
fig = plt.plot(days, data_coll)
plt.legend(fig,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 200 days pre-inf %d days' %(1/kappa))
plt.xlabel('Time(Days)')
plt.ylabel('Population')

kappa1 = 1/5; 					# Pre-infectious period 5 days
kappa2 = 1/20; 					# Pre-infectious period 20 days

data_coll1 = []
I1 = [I0]
data_coll1.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll1[i-1][2], 0, 0, 0],
                        [beta*data_coll1[i-1][2], -kappa1, 0, 0],
                        [0, kappa1, -alpha, 0],
                        [0, 0, alpha, 0]])
    data_coll1.append(data_coll1[i-1] + dt*np.matmul(diff_mat,data_coll1[i-1]))
    I1.append(data_coll1[i][2])

plt.figure()
fig1 = plt.plot(days, data_coll1)
plt.legend(fig1,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 200 days pre-inf %d days' %(1/kappa1))
plt.xlabel('Time(Days)')
plt.ylabel('Population')

data_coll2 = []
I2 = [I0]
data_coll2.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll2[i-1][2], 0, 0, 0],
                        [beta*data_coll2[i-1][2], -kappa2, 0, 0],
                        [0, kappa2, -alpha, 0],
                        [0, 0, alpha, 0]])
    data_coll2.append(data_coll2[i-1] + dt*np.matmul(diff_mat,data_coll2[i-1]))
    I2.append(data_coll2[i][2])

plt.figure()
fig2 = plt.plot(days, data_coll2)
plt.legend(fig2,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 200 days pre-inf %d days' %(1/kappa2))
plt.xlabel('Time(Days)')
plt.ylabel('Population')

fig3, ax = plt.subplots()
ax.plot(days, I, label='%d days' %(1/kappa))
ax.plot(days, I1, label='%d days' %(1/kappa1))
ax.plot(days, I2, label='%d days' %(1/kappa2))
ax.set_title('SEIR model with different pre-infectious periods')
ax.set_xlabel('Time(Days)')
ax.set_ylabel('Population')
ax.legend()

## Consider natural death
LE = 70*365         # Life expectancy
mu = 1/LE           # Birth rate = death rate

data_coll3 = []
I3 = [I0]
data_coll3.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll3[i-1][2], mu, mu, mu],
                        [beta*data_coll3[i-1][2], -kappa-mu, 0, 0],
                        [0, kappa, -alpha-mu, 0],
                        [0, 0, alpha, -mu]])
    data_coll3.append(data_coll3[i-1] + dt*np.matmul(diff_mat,data_coll3[i-1]))
    # I3.append(data_coll3[i][2])

fig4, ax = plt.subplots()
ax.plot(days, data_coll,'--',label=['susceptible','pre-infectious','infectious','recover'])
ax.plot(days, data_coll3,label=['susceptible','pre-infectious','infectious','recover'])
ax.set_title('SEIR model for 200 days pre-inf %d days' %(1/kappa))
ax.set_xlabel('Time(Days)')
ax.set_ylabel('Population')
ax.legend()

trans_period = 10*365
days = range(0,trans_period+1)
days = np.array(days)/365
data_coll4 = []
I4 = [I0]
data_coll4.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll4[i-1][2], mu, mu, mu],
                        [beta*data_coll4[i-1][2], -kappa-mu, 0, 0],
                        [0, kappa, -alpha-mu, 0],
                        [0, 0, alpha, -mu]])
    data_coll4.append(data_coll4[i-1] + dt*np.matmul(diff_mat,data_coll4[i-1]))
    # I4.append(data_coll4[i][2])

plt.figure()
fig5 = plt.plot(days, data_coll4)
plt.legend(fig5,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 10 years pre-inf %d days' %(1/kappa))
plt.xlabel('Time(years)')

trans_period = 50*365
days = range(0,trans_period+1)
days = np.array(days)/365
data_coll5 = []
I5 = [I0]
data_coll5.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll5[i-1][2], mu, mu, mu],
                        [beta*data_coll5[i-1][2], -kappa-mu, 0, 0],
                        [0, kappa, -alpha-mu, 0],
                        [0, 0, alpha, -mu]])
    data_coll5.append(data_coll5[i-1] + dt*np.matmul(diff_mat,data_coll5[i-1]))
    # I5.append(data_coll5[i][2])

plt.figure()
fig6 = plt.plot(days, data_coll5)
plt.legend(fig6,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 50 years pre-inf %d days' %(1/kappa))
plt.xlabel('Time(years)')

trans_period = 100*365
days = range(0,trans_period+1)
days = np.array(days)/365
data_coll6 = []
I6 = [I0]
data_coll6.append(np.array([S0, E0, I0, R0]))
for i in range(1,trans_period+1):
    diff_mat = np.array([[-beta*data_coll6[i-1][2], mu, mu, mu],
                        [beta*data_coll6[i-1][2], -kappa-mu, 0, 0],
                        [0, kappa, -alpha-mu, 0],
                        [0, 0, alpha, -mu]])
    data_coll6.append(data_coll6[i-1] + dt*np.matmul(diff_mat,data_coll6[i-1]))
    # I6.append(data_coll6[i][2])

plt.figure()
fig6 = plt.plot(days, data_coll6)
plt.legend(fig6,['susceptible','pre-infectious','infectious','recover'])
plt.title('SEIR model for 100 years pre-inf %d days' %(1/kappa))
plt.xlabel('Time(years)')

plt.show()