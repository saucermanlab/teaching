import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
p = 0.5
y0 = [1]
tspan = [0,10]
t = np.linspace(0,10,100)

def f(t,y,p):
    dydt = - p*y
    return dydt

sol = scipy.integrate.solve_ivp(f,tspan,y0,t_eval=t,args=(p,))

plt.plot(sol.t,sol.y.T)
plt.show()