
#----------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import figstyle

#----------------------------------------------------------------------------------

n_blocks = 20
nstates = 21
dP_mean = np.zeros(nstates)
dV_mean = np.zeros(nstates)
P_mean = np.zeros(nstates)
V_mean = np.zeros(nstates)
P = np.zeros((nstates,n_blocks))
V = np.zeros((nstates,n_blocks))

n_samples = 1000
P_new = np.zeros((nstates,n_samples))
V_new = np.zeros((nstates,n_samples))
integral = np.zeros(n_samples)

for state in range(nstates):
            
    data = pd.read_csv(f'{state}_output_{state}.txt', sep=" ")
    
    for i in range(n_blocks):
        P[state,i] = np.mean(data['Press'].values[250*i:250*(i+1)])
        V[state,i] = np.mean(data['Volume'].values[250*i:250*(i+1)])
    
    P_mean[state] = np.mean(P[state,:])
    V_mean[state] = np.mean(V[state,:])
    dP_mean[state] = 2*np.std(P[state,:])   
    dV_mean[state] = 2*np.std(V[state,:])  

    P_new[state,:] = np.random.normal(P_mean[state], dP_mean[state], n_samples)
    V_new[state,:] = np.random.normal(V_mean[state], dV_mean[state], n_samples)         
    
#----------------------------------------------------------------------------------
conv = 1e-27*24.21726e-3 #atm A3 --> atm L --> kcal
Navogadro = 6.02214076e23 #atomos/mol
Natomos = 500.0 #atomos

from scipy.integrate import simps

for i in range(n_samples):
    integral[i] = simps(P_new[:,i],V_new[:,i])*conv*Navogadro/Natomos*4.184

final = np.mean(integral)
dfinal = 2*np.std(integral) 
print(final,dfinal)
#----------------------------------------------------------------------------------

fig = plt.figure(figsize=(2.5,2.5), dpi=300)
plt.plot(V_mean,P_mean,'k-')
plt.errorbar(V_mean, P_mean, yerr = dP_mean, fmt="o")
plt.ylabel(r'$P$ (atm)')
plt.xlabel(r'$V$ (A$^3$)')
plt.savefig('figura.png', dpi=300)
plt.show()


#----------------------------------------------------------------------------------



