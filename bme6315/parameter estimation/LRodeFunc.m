function dydt = LRodeFunc(t,y,params)

% Parameter assignments
L = params(1);
kon = params(2);
koff = params(3);

% Variable assignments
R = y(1);
LR = y(2);

% Reaction fluxes
Jon = kon*L*R;
Joff = koff*LR;

% Differential equations
dR = Joff - Jon;
dLR = Jon - Joff;
dydt = [dR; dLR];   % repackage ODE's into a single array