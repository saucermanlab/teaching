function motifsPosFeedback() 
% positive feedback
% update 1/10/2020 by JS

% species parameters 
speciesNames = {'A','B'}; 
tau = [1, 1, ]; 
ymax = [0.972, 1, ];
 
% reaction parameters
% r1: input => A
% r2: A => B
% r3: B => A
w = [0.5, 1, 1, ]; 
n = [0, 4, 4, ]; 
EC50 = [0, 0.5, 0.5, ]; 
rpar = [w;n;EC50];
 
params = {rpar,tau,ymax,speciesNames};
 
% Run single simulation 
y0 = [0, 0, ]; 
tspan = [0 50]; 
options = []; 
[t,y] = ode23(@ODEfun,tspan,y0,options,params); 
figure(1); 
plot(t,y(:,1)); 
xlabel('Time (sec)'); ylabel('Normalized activity');
 
% stability analysis based on approximated Jacobian
figure(2);
subplot(1,2,1);
x=[0:0.01:1];
n=5;
plot(x,(x.^n)./(0.5^n+x.^n),(x.^n)./(0.5^n+x.^n),x) % plot nullclines
title('Nullclines'); xlabel('A'); ylabel('B'); legend('dA/dt = 0','dB/dt = 0');

subplot(1,2,2);
y(end,:)
Jac = approxJacobian(@ODEfun,0,y(end,:),params);
% Jac = approxJacobian(@ODEfun,0,[.5,.5],params); % unstable point w/ feedback
% Jac = approxJacobian(@ODEfun,0,[0.9195,0.9198],params); % stable point w/ feedback
imagesc(Jac);
title('Jacobian matrix at final time');
colorbar
disp('Eigenvalues of the Jacobian:'); 
eig(Jac)
% stable if both eigenvalues have negative real parts

function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
A = 1; 
B = 2;  
input = rpar(1,1)*(t<5);

dydt = zeros(2,1); 
dydt(A) = (OR(input,act(y(B),rpar(:,3)))*ymax(A) - y(A))/tau(A); 
% dydt(A) = (input*ymax(A) - y(A))/tau(A); % no feedback
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
    
function Jac = approxJacobian(ODEfunc,t0,y0,p);
% estimate Jacobian matrix by finite difference
    epsilon = 0.001; % perturb params by this relative amount
    func0 = ODEfunc(t0,y0,p);       % evaluate ODEfunc at baseline
    Jac = zeros(numel(y0));           % initializes Jacobian matrix
    for i=1:numel(y0),
        y1 = y0;
        dy = max(epsilon*y0(i),sqrt(eps));  % minimum value of 1.49e-8
        y1(i) = y0(i) + dy;                 % perturb y(i) by dy
        func1 = ODEfunc(t0,y1,p);           % evaluate ODEfunc with perturbed y
        Jac(:,i) = (func1-func0)/dy;        % update Jacobian matrix
    end