# exampleOptimization.py
# JS 10/3/2021
import numpy as np
import scipy.optimize, scipy.linalg
import matplotlib.pyplot as plt

def biexp(p,x):
    a, b, c = p
    return c*np.exp(-a*x)*(1-np.exp(-b*x))

def objectiveFunction(p,xdata,ydata):
    ysim = biexp(p,xdata)
    error = ydata - ysim
    return error

# generate simulated data for fitting
p0 = [1,2,3]
xdata = np.linspace(0,5)
ydata = biexp(p0,xdata) + 0.1*np.random.normal(size=len(xdata)) # add noise std dev 0.1

# run optimization
result = scipy.optimize.least_squares(objectiveFunction, p0, args=(xdata,ydata))
paramEst = result.x
covariance = result.cost*scipy.linalg.inv((np.matmul(result.jac.T,result.jac)))/len(ydata)

# plot the results
plt.plot(xdata,ydata,'o',xdata,biexp(paramEst,xdata),'-')
plt.show()

# generate parameter ensemble
n = 10;
paramEnsemble = np.random.multivariate_normal(paramEst,covariance,n)

paramStd = np.std(paramEnsemble,axis=0)

# plot parameter estimates
plt.bar(['a','b','c'],paramEst)
plt.errorbar(['a','b','c'],paramEst,paramStd)
plt.show()

# plot simulatons with the parameter ensemble
y = biexp(paramEst,xdata)
plt.plot(xdata,y)
for i in range(n):
    y = biexp(paramEnsemble[i],xdata)
    plt.plot(xdata,y)
plt.show()