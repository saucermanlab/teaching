% LRparamEst2.m
% parameter estimation for 2 fitted parameters

function output = LRparamEst2
%% Parameter assignments
L = 0.9;
kon = 1.1;    % fitting these, put best guess
koff = 0.1; % fitting these, put best guess
% params = [L,kon,koff];

% Simulation parameters
y0 = [1; 0];    % initial conditions
tspan = [0 6];
simOptions = [];

% Load experimental data (synthetic, using L = 1, kon = 1, koff = 0.1)
data = [0:6; 0,0.58,0.79,0.93,0.88,0.94,0.80]';  % 1st column: times; 2nd column: LR measurements

paramsFixed = [koff];     % designate parameters that will not be estimated
lowerbnd = [];              % optional lower bound on parameter estimate
upperbnd = [];              % optional upper bound on parameter estimate
params0 = [L,kon];            % initial guesses of parameter estimate
% params0 = (upperbnd-lowerbnd).*rand(size(params0))+lowerbnd;    % use random initial guesses (overwrites above line)

optimOptions = optimset('Display','iter');
objFcn = @(params0) objectiveFcn(params0,paramsFixed,tspan,y0,data,simOptions);
[paramsEst,resnorm,residual,exitflag,output,lambda,Jacobian]=lsqnonlin(objFcn,params0,lowerbnd,upperbnd,optimOptions);
disp(['Algorithm: ',output.algorithm]);
disp(['Exit flag: ',num2str(exitflag),', resnorm = ',num2str(resnorm)]);
paramsVariance = resnorm*inv(Jacobian'*Jacobian)/length(data(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
paramsStd = sqrt(diag(paramsVariance));   % standard deviations for parameter estimates
paramsCI = 1.96*paramsStd;                % mean +/- 1.96*stdev in 95% confidence interval
output = [paramsEst;paramsCI']

%% Run test of parameter estimation
params = [paramsEst(1),paramsEst(2),koff];
[tSingle,ySingle] = ode23(@LRodeFunc,tspan,y0,simOptions,params);
LRsingle = ySingle(:,2);
figure(1);
subplot(1,3,1);
plot(tSingle,LRsingle,data(:,1),data(:,2),'o');
xlabel('Time (sec)'); ylabel('LR (\muM)');
legend('model','data');
subplot(1,3,2);
barweb([[1.0; 1.0],paramsEst'],[[0; 0], paramsStd],[],{'L','kon'}); % note need 'barweb' from MatlabCentral
axis([0.5 2.5 0 1.5])
legend('Original','Estimated'); ylabel('Estimated parameters');

%% Plot entire SSE surface to verify that paramEst is in fact a minimum
LEst = paramsEst(1);
konEst = paramsEst(2);
SSEest = sum(objFcn(paramsEst).^2);
param1Range = [0.5:.05:2];
param2Range = [0.5:.05:2];
for i = 1:length(param1Range)
    for j = 1:length(param2Range)
        param = [param1Range(i),param2Range(j)];
        SSE(i,j) = sum(objFcn(param).^2);
    end
end
subplot(1,3,3);
imagesc(param1Range,param2Range,SSE); 
xlabel('L'); ylabel('kon'); 
c = colorbar;
c.Label.String = 'SSE';

function error = objectiveFcn(paramsEst,paramsFixed,tspan,y0,data,simOptions)
% objective function for performing least squares minimization
% This error function is weighted by the std dev of the experimental
% measurements, which affects confidence intervals and relative weighting 
% of multi-objective function estimates.
L = paramsEst(1);
kon = paramsEst(2);
koff = paramsFixed(1);
params = [L,kon,koff];

[tSim,ySim] = ode23(@LRodeFunc,tspan,y0,simOptions,params);   % run the model for a given sample of paramsEst
ySimInterp = interp1(tSim,ySim,data(:,1)); % resample simulated data before performing subtraction
stdevExp = .05; % standard deviation of experimental data used in weights (set to 1 by default)
error = (1/stdevExp).*(ySimInterp(:,2) - data(:,end));     % note that LSQNONLIN computes the sum of squares for this automatically
