# LRparamEst.py
# Parameter estimation for a ligand-receptor model
# 10/1/2021 by JS
import numpy as np
from scipy import optimize, integrate, interpolate, linalg
import matplotlib.pyplot as plt

# Parameter assignments
L = 1
kon = 1
koff = 0.1
params = [L,kon,koff]

def ODEfunc(t,y,params):
    L,k1f,k1r = params
    R, LR = y    
    react1 = k1f*L*R - k1r*LR       # [uM/s] L+R<->LR    
    dR = -react1            # [uM/s] free receptor
    dLR = react1            # [uM/s] ligand-receptor complex    
    return [dR,dLR]

# Simulation parameters
y0 = [1,0]
tspan = [0,6]
paramsFixed = koff
p0 = [0.9, 1.1]

# Load experimental data
xdata = np.arange(7)
ydata = [0,0.58,0.79,0.93,0.88,0.94,0.80]

# Run parameter estimation
def objectiveFunction(paramsEst,paramsFixed,tspan,xdata,ydata):
    koff = paramsFixed
    L, kon = paramsEst
    params = [L,kon,koff]
    sol = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
    finterp = interpolate.interp1d(sol.t,sol.y[1,:]) # create an interpolation function
    error = ydata - finterp(xdata)       #  compute error at resampled locations
    return error

result = optimize.least_squares(objectiveFunction, p0, args=(paramsFixed,tspan,xdata,ydata))
paramsEst = result.x

# Run test of parameter estimation, shows a good fit
params = [paramsEst[0],paramsEst[1],koff]
sol = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
plt.plot(xdata,ydata,'o',sol.t,sol.y[1,:].T,'-')
plt.show()

# generate parameter ensemble
# note in this particular example the covariance is huge. If the parameters had
# been more identifiable it would be much smaller.
covariance = result.cost*linalg.inv((np.matmul(result.jac.T,result.jac)))/len(ydata)
n = 10;
paramEnsemble = np.random.multivariate_normal(paramsEst,covariance,n)
paramsStd = np.std(paramEnsemble,axis=0)

# plot parameter estimates
plt.bar(['L','kon'],paramsEst)
plt.errorbar(['L','kon'],paramsEst,paramsStd)
plt.show()