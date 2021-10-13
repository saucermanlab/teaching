# toyModuleCompute.py
# 10/1/2021 by JS
# Modeling the following reactions:
# 1) L+R<->LR, 2) LR+E<->LRE, 3) LRE catalyzes S->P, 4) P->S
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import toyModelODEfunc

# define parameters
k1f = 1    # [uM^-1 s^-1] react1 forward rate constant
k1r = 1    # [s^-1] react1 reverse rate constant
k2f = 1    # [uM^-1 s^-1] react2 forward rate constant
k2r = 1    # [s^-1] react2 reverse rate constant
kcat = 1   # [s^-1] catalytic rate constant for enzyme
Km = 1     # [uM] Michaelis constant for enzyme
k4 = 1     # [s^-1] react4 rate constant
L = 1      # [uM] concentration of ligand
Rtot = 1   # [uM] total concentration of receptor
Etot = 1   # [uM] total concentration of enzyme
Stot = 1   # [uM] total concentration of substrate
params = [k1f,k1r,k2f,k2r,kcat,Km,k4,L]

y0 = [Rtot, 0, Etot, 0, Stot, 0]
tspan = [0, 10]

# Run a single simulation
sol = scipy.integrate.solve_ivp(toyModelODEfunc.ODE,tspan,y0,args=(params,))

#plt.subplot(221)
plt.plot(sol.t,sol.y.T)
plt.show()

# Run dose response over a range of total ligand concentrations
paramRange = np.logspace(-2,2,20)
yfinal = []
for param in paramRange:
    L = param
    params = [k1f,k1r,k2f,k2r,kcat,Km,k4,L]
    sol = scipy.integrate.solve_ivp(toyModelODEfunc.ODE,tspan,y0,args=(params,))
    yfinal.append(sol.y[-1,-1])
    
#plt.subplot(222)
plt.semilogx(paramRange,yfinal)
plt.show()
# not getting the subplots to work correctly, but simulation is correct.
    