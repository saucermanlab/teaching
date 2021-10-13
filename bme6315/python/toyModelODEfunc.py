# toyModuleODEfunc.py
# 10/1/2021 by JS

def ODE(t,y,params):
    k1f,k1r,k2f,k2r,kcat,Km,k4,L = params
    
    R, LR, E, LRE, S, P = y
    
    # reaction rates
    react1 = k1f*L*R - k1r*LR       # [uM/s] L+R<->LR
    react2 = k2f*LR*E - k2r*LRE     # [uM/s] LR+E<->LRE
    react3 = kcat*LRE*S/(Km+S)      # [uM/s] LRE catalyzes S->P
    react4 = k4*P                   # [uM/s] P->S
    
    # differential equations
    dR = -react1            # [uM/s] free receptor
    dLR = react1-react2     # [uM/s] ligand-receptor complex
    dE = -react2            # [uM/s] free enzyme
    dLRE = react2           # [uM/s] ligand-receptor-enzyme complex
    dS = react4-react3      # [uM/s] substrate
    dP = react3-react4      # [uM/s] product
    
    return [dR,dLR,dE,dLRE,dS,dP]