function compartmentscAMP() 
% 2 compartment model of cAMP gradients

% parameters
vmaxAC = 1;     % [uM^-1 s^-1]
vmaxPDE = 1;    % [uM^-1 s^-1]
KmPDE = 2;      % [uM]
kcAMPdiff = 1;  % [s^-1]
volRatio = 0.1;
params = {vmaxAC,vmaxPDE,KmPDE,kcAMPdiff,volRatio};
 
% Run single simulation 
y0 = [0, 0, ]; 
tspan = [0 20]; 
options = []; 
[t,y] = ode23(@ODEfun,tspan,y0,options,params); 

plot(t,y); 
xlabel('Time (sec)'); ylabel('cAMP (\muM)'); 
legend('membrane','cytoplasm');
 
function dydt=ODEfun(t,y,params) 
% Assign names for parameters and variables
[vmaxAC,vmaxPDE,KmPDE,kcAMPdiff,volRatio]=params{:}; 
cAMPmem = y(1); 
cAMPcyt = y(2);  

dydt_cAMPmem = vmaxAC - kcAMPdiff*(cAMPmem-cAMPcyt);
dydt_cAMPcyt = kcAMPdiff*volRatio*(cAMPmem-cAMPcyt) - vmaxPDE*cAMPcyt/(KmPDE+cAMPcyt);
dydt = [dydt_cAMPmem; dydt_cAMPcyt];