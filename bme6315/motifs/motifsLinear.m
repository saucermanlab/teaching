function motifsLinear() 
% linear pathway
% updated 1/10/2020 by JS

% species parameters 
speciesNames = {'A','B'}; 
tau = [0.1, 5, ]; 
ymax = [1, 1, ];
 
% reaction parameters
% r1: input => A
% r2: A => B
w = [1, 1, ]; 
n = [0, 1, ]; 
EC50 = [0, 0.5,]; 
rpar = [w;n;EC50];
 
params = {rpar,tau,ymax,speciesNames};
 
% Run single simulation 
y0 = [0, 0, ]; 
tspan = [0 15]; 
options = []; 
[t,y] = ode23(@ODEfun,tspan,y0,options,params); 
% figure; 
plot(t,y); 
xlabel('Time (sec)'); ylabel('Normalized activity'); legend(speciesNames); 
 
function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
A = 1; 
B = 2;  
% input = rpar(1,1)*(t<5);
input = rpar(1,1)*(mod(t,2)<1);
dydt = zeros(2,1); 
dydt(A) = (input*ymax(A) - y(A))/tau(A); 
dydt(B) = (act(y(A),rpar(:,2))*ymax(B) - y(B))/tau(B); 
   
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    fact = w.*(x.^n)./(EC50.^n + x.^n); 
    
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = rpar(1) - act(x,rpar);
 
function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
 
function z = AND(x,y) 
% AND logic gate, multiplying the reactants together 
    z = x*y;