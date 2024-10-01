
#-----------------------------------------------------------------------------------------

import mics as mx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import figstyle

#-----------------------------------------------------------------------------------------

conv = 1e-27*24.21726e-3 #atm A3 --> atm L --> kcal
Navogadro = 6.02214076e23 #atomos/mol
npoints = 201 
kb = 1.9872e-3 #kcal/mol/K
temp = np.array([590,600,610,620,630,640,650,660])
nstates = len(temp)
deltaG = np.zeros(npoints)
deltaG_upper = np.zeros(npoints)
deltaG_lower = np.zeros(npoints)

deltaf_alch = 0.47/4.184/(kb*660.0)

deltaHp = np.zeros(nstates)

#-----------------------------------------------------------------------------------------
samples_liq = mx.pooledsample()
samples_sol = mx.pooledsample()
for state in range(nstates):

    print(f'state {state}')

    kwargs = {}
    kwargs['T'] = temp[state]
    kwargs['beta'] = 1/(kb*temp[state])  
    
    data = pd.read_csv(f'../liquid/{int(temp[state])}/output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(10000), inplace=True)
    data.index = np.arange(0, len(data))
    datal = data
    samples_liq.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))

    data = pd.read_csv(f'../solid/{int(temp[state])}/output.out', sep=" ")
    data['U_pot'] = data['PotEng']
    data['PV'] = 1.0*data['Volume']*conv*Navogadro
    data.drop(index=range(20000), inplace=True)
    data.index = np.arange(0, len(data))
    datas = data
    samples_sol.append(mx.sample(data, 'beta*(U_pot + PV)', acfun='Temp', **kwargs))
    
    deltaHp[state] = (np.mean(datas['Enthalpy'])- np.mean(datal['Enthalpy']))/500.0

samples_liq.subsampling(integratedACF=True)
samples_sol.subsampling(integratedACF=True)

mixture_liq = mx.mixture(samples_liq, engine=mx.MBAR())
mixture_sol = mx.mixture(samples_sol, engine=mx.MBAR())

results_liq = mixture_liq.free_energies()
results_sol = mixture_sol.free_energies()

#-----------------------------------------------------------------------------------------
# REWEIGHTING

variables = dict(beta=[], T=[])
temp_new = np.linspace(temp.min(), temp.max(), npoints)

for point in range(npoints):

    variables['beta'].append( 1/(kb*temp_new[point])  )
    variables['T'].append( temp_new[point]  )


reweighting_liq = mixture_liq.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

reweighting_sol = mixture_sol.reweighting(
    potential='beta*(U_pot + PV)',
    conditions=pd.DataFrame(variables)
    )

Nat = 500.0

reweighting_liq['f'] = reweighting_liq['f']/Nat
reweighting_sol['f'] = reweighting_sol['f']/Nat
results_liq['f'] = results_liq['f']/Nat
results_sol['f'] = results_sol['f']/Nat

#-----------------------------------------------------------------------------------------

for point in range(npoints):
    deltaG[point] = (reweighting_sol['f'][point] - reweighting_liq['f'][point] - \
                     results_sol['f'][7] + results_liq['f'][7] + \
                         deltaf_alch)*kb*reweighting_liq['T'][point]
#-----------------------------------------------------------------------------------------

indice = np.where(np.diff(np.sign(deltaG)))[0] 
raiz = (reweighting_liq['T'].values[indice[0]] + reweighting_liq['T'].values[indice[0]+1] )/2
print(raiz)

#-----------------------------------------------------------------------------------------
fig = plt.figure(figsize=(3.0, 4.0),  dpi=300)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_ylabel(r'$\beta G - \beta_\text{ref} G_\text{ref}$')
ax2.set_ylabel(r'$\Delta G_\textsc{l,s}$ (kJ/mol)')
ax2.set_xlabel(r'$T$ (K)')
ax1.set_xticklabels([])

p1, = ax1.plot(reweighting_liq['T'], reweighting_liq['f'], 'c-', label='Liquid')
p2, = ax1.plot(results_liq['T'], results_liq['f'], 'co', label='')
p3, = ax1.plot(reweighting_sol['T'], reweighting_sol['f'], 'm-', label='Solid')
p4, = ax1.plot(results_sol['T'], results_sol['f'], 'mv', label='')

from matplotlib.legend_handler import HandlerTuple
ax1.legend([(p1, p2), (p3, p4)], ['Liquid', 'Solid'], 
               handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2)

ax2.plot(reweighting_liq['T'], deltaG*4.184, 'k-', label='')
ax2.plot(np.linspace(590, 660, 100), np.linspace(0, 0, 100), 'k--', label='')
ax2.legend()

fig.align_ylabels([ax1, ax2])
plt.tight_layout()

plt.savefig('figura_mbar.png', dpi=300)
plt.show()



