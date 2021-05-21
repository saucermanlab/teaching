function [dydt, algvars] = toyModelODEfunc(t,y,params)
% ODEfunc defines the system of ODEs describing the model
% 1/10/2020 JS: updates to simplify

% Assign names for parameter values and state variables
[k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L] = params{:};

R = y(1);
LR = y(2);
E = y(3);
LRE = y(4);
S = y(5);
P = y(6);
Pnuc = y(7);

% Reaction rates
react1 = k1f*L*R - k1r*LR;      % [uM/s] L+R<->LR
react2 = k2f*LR*E - k2r*LRE;    % [uM/s] LR+E<->LRE
react3 = kcat*LRE*S/(Km+S);     % [uM/s] LRE catalyzes S->P
react4 = k4*P;                  % [uM/s] P->S
react5 = k5*(P-Pnuc);           % [uM/s] Pnuc <-> P

% Differential equations;
dR = -react1;             % [uM/s] free receptor
dLR = react1-react2;            % [uM/s] ligand-receptor complex
dE = -react2;                   % [uM/s] free enzyme
dLRE = react2;                  % [uM/s] ligand-receptor-enzyme complex
dS = react4-react3;             % [uM/s] substrate
dP = react3-react4-react5;      % [uM/s] product
dPnuc = react5/Vratio;          % [uM/s] nuclear product
dydt = [dR;dLR;dE;dLRE;dS;dP;dPnuc]; % Reassemble differential equations
algvars = [react1,react2,react3,react4,react5]; % optional for seeing fluxes