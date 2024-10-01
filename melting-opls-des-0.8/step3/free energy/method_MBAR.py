
#-----------------------------------------------------------------------------------------

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import figstyle

#-----------------------------------------------------------------------------------------

def h1(x):
    return eta + x*(1-eta)

def h1p(x): 
    return (1-eta)

def h2(x):
    return (eta + x*(1-eta))**2

def h2p(x): 
    return 2*(eta + x*(1-eta))*(1-eta)

#-----------------------------------------------------------------------------------------

npoints = 201 
T = 660.0 #K
kb = 1.9872e-3 #kcal/mol/K
eta = 0.1
lambda_v = np.array([0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
nstates = len(lambda_v)

kT = T*kb
U = 'h1*U_vdwl + h2*U_coulk'
Up = 'h1p*U_vdwl + h2p*U_coulk' 

#-----------------------------------------------------------------------------------------

dUv = np.zeros(nstates)
dUc = np.zeros(nstates)
samples = mx.pooledsample()
mx.verbose=True
for state in range(nstates):
    print(f'state {state}')
    lbda = lambda_v[state]
    data = pd.read_csv(f'../{state}/output.out', sep=" ")

    data['dU_vdwl'] = data['c_lj']/h1(lbda)*h1p(lbda) + data['E_tail']/h1(lbda)*h1p(lbda)
    data['U_vdwl'] = data['c_lj']/h1(lbda) + data['E_tail']/h1(lbda)
    data['dU_coulk'] = data['c_coul']/h2(lbda)*h2p(lbda) + data['E_long']/h2(lbda)*h2p(lbda)
    data['U_coulk'] = data['c_coul']/h2(lbda) + data['E_long']/h2(lbda)
    
    data.drop(index=range(5000), inplace=True)
    data.index = np.arange(0, len(data))
    kwargs = {'lbda': lbda, 'beta': 1.0/kT}
    kwargs['h1'] = h1(lbda) 
    kwargs['h2'] = h2(lbda) 
    samples.append(mx.sample(data, f'beta*({U})', acfun='Temp', **kwargs))

    dUv[state] = np.mean(data['dU_vdwl'].values)
    dUc[state] = np.mean(data['dU_coulk'].values)

samples.subsampling(integratedACF=True)

mixture = mx.mixture(samples, engine=mx.MBAR())

results = mixture.free_energies()
results['F'] = results['f']*kT
results['dF'] = results['df']*kT

#-----------------------------------------------------------------------------------------

variables = dict(lbda=[], h1=[], h2=[], h1p=[], h2p=[])
for point in range(npoints):
    lbda = point/(npoints - 1)
    variables['lbda'].append(lbda)
    variables['h1'].append(h1(lbda))
    variables['h1p'].append(h1p(lbda))
    variables['h2'].append(h2(lbda))
    variables['h2p'].append(h2p(lbda))

properties = dict(
    Up=f'{Up}'
    )

combinations = dict(
    F='kT*f'   
    )

reweighting = mixture.reweighting(
    potential=f'beta*({U})',
    properties=properties,
    combinations=combinations,
    conditions=pd.DataFrame(variables), 
    beta=1/kT,
    kT=kT 
    )

#-----------------------------------------------------------------------------------------
Nat = 500.0

deltaG1 = results['F'].iloc[-1]/Nat*4.184
d_deltaG1 = results['dF'].iloc[-1]/Nat*4.184

print(f'Delta G = {deltaG1} {d_deltaG1} kJ/mol')

fig, ax = plt.subplots(2, 1, figsize=(2.5, 5.0),  dpi=300)
ax[0].plot(reweighting['lbda'], reweighting['F']/Nat*4.184, 'k-', label='Interpolated')
ax[0].plot(results['lbda'], results['F']/Nat*4.184, 'bo', label='Simulated')
ax[1].plot(reweighting['lbda'], reweighting['Up']/Nat*4.184, 'k-', label='')
ax[1].plot(results['lbda'], dUv/Nat*4.184 + dUc/Nat*4.184, 'bo', label='')
ax[0].set_ylabel('$A$ (kJ/mol)')
ax[1].set_ylabel(r'$\left\langle \frac{\partial U}{\partial \lambda} \right\rangle$ (kJ/mol)')
ax[1].set_xlabel('$\lambda$')
ax[0].set_xticklabels([])
plt.savefig('figure_MBAR.png', dpi=300)
plt.show()






