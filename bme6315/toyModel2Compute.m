function toyModelCompute
% simple toy model for testing template
% 1/10/2020 JS: updates
% 6/02/2011 - Jeff Saucerman
%
% Modeling the following reactions:
% 1) L+R<->LR, 2) LR+E<->LRE, 3) LRE catalyzes S->P, 4) P->S, 5) P <-> Pnuc

%% define parameters
k1f = 1;    % [uM^-1 s^-1] react1 forward rate constant
k1r = 1;    % [s^-1] react1 reverse rate constant
k2f = 1;    % [uM^-1 s^-1] react2 forward rate constant
k2r = 1;    % [s^-1] react2 reverse rate constant
kcat = 1;   % [s^-1] catalytic rate constant for enzyme
Km = 1;     % [uM] Michaelis constant for enzyme
k4 = 1;     % [s^-1] react4 rate constant
k5 = 1;
Vratio = 0.1;

L = 1;   % [uM] concentration of ligand
Rtot = 1;   % [uM] total concentration of receptor
Etot = 1;   % [uM] total concentration of enzyme
Stot = 1;   % [uM] total concentration of substrate
params = {k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L};

%% Run single simulation
% 'R','LR','E','LRE','S','P','Pnuc'
y0 = [Rtot; 0; Etot; 0; Stot; 0; 0];   % or you could load saved data like: y0 = load('yfinal.dat');
tspan = [0 10];
options = [];%odeset('MaxStep',5e-3);
[t,y] = ode23(@toyModel2ODEfunc,tspan,y0,options,params);

yfinal = y(end,:)';
% save -ascii 'yfinal.dat' yfinal;    % save final values to a file

% plot timecourse
subplot(1,3,1);
plot(t,y);
xlabel('Time (sec)'); ylabel('y(t) (\muM)'); legend('R','LR','E','LRE','S','P','Pnuc');

% optional: re-evaluate ODE function after solving ODEs to calculate algebraic variables 
for tstep=1:length(t),
    [~,algvars(tstep,:)]=toyModel2ODEfunc(t(tstep),y(tstep,:),params);
end
subplot(1,3,2);
plot(t,algvars(:,1),t,algvars(:,2),t,algvars(:,3),t,algvars(:,4),t,algvars(:,5));
xlabel('Time (sec)'); ylabel('fluxes(t) (\muM/s)');
legend('react1','react2','react3','react4','react5');

%% Run dose response over a range of total ligand concentrations
paramRange = 10.^[-2:.1:2];
for i=1:length(paramRange)
    L = paramRange(i);
    params = {k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L};
    y0 = [Rtot; 0; Etot; 0; Stot; 0; 0];
    tspan = [0 10];
    options = [];
    [t,y] = ode15s(@toyModel2ODEfunc,tspan,y0,options,params);
    Pnucfinal(i) = y(end,end);
end
subplot(1,3,3);
semilogx(paramRange,Pnucfinal);
xlabel('Ligand (\muM)'); ylabel('Steady state nuclear Product (\muM)');