function motifsNegFeedbackOscillations() 
% negative feedback
% updates 1/10/2020 by JS

% species parameters 
speciesNames = {'A','B','C'}; 
tau = [1, 1, 1]; 
ymax = [1, 1, 1];
 
% reaction parameters
% r1: input
% r2: A => B
% r3: B => C
% r4: input & !C => A
w = [1, 1, 1, 1]; 
n = [0, 5, 5, 5]; 
EC50 = [0, 0.5, 0.5, 0.5]; 
rpar = [w;n;EC50];
 
params = {rpar,tau,ymax,speciesNames};
 
% Run single simulation 
y0 = [0, 0, 0]; 
tspan = [0 50]; 
options = []; 
[t,y] = ode23(@ODEfun,tspan,y0,options,params); 
% figure; 
plot(t,y); 
xlabel('Time (sec)'); ylabel('Normalized activity'); legend(speciesNames);
 
function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax]=params{:}; 
A = 1; 
B = 2;  
C = 3;
input = rpar(1,1);
% input = rpar(1,1)*(mod(t,2)<1);
dydt = zeros(2,1); 
dydt(A) = (AND(input,inhib(y(C),rpar(:,4)))*ymax(A) - y(A))/tau(A); 
dydt(B) = (act(y(A),rpar(:,2))*ymax(B) - y(B))/tau(B); 
dydt(C) = (act(y(B),rpar(:,3))*ymax(C) - y(C))/tau(C);
   
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