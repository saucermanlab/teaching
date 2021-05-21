% LRparamEst.m
% parameter estimation for 1 fitted parameter

function LRparamEst
%% Parameter assignments
L = 1;
kon = 1;
koff = 0.1;
params = [L,kon,koff];

% Simulation parameters
y0 = [1; 0];    % initial conditions
tspan = [0 6];
simOptions = [];

% Load experimental data (synthetic, using L = 1, kon = 1, koff = 0.1)
data = [0:6; 0,0.58,0.79,0.93,0.88,0.94,0.80]';  % 1st column: times; 2nd column: LR measurements

paramsFixed = [L,koff];     % designate parameters that will not be estimated
lowerbnd = [];              % optional lower bound on parameter estimate
upperbnd = [];              % optional upper bound on parameter estimate
params0 = [1.5];            % initial guesses of parameter estimate
% params0 = (upperbnd-lowerbnd).*rand(size(params0))+lowerbnd;    % use random initial guesses (overwrites above line)

optimOptions = optimset('Display','iter');
objFcn = @(params0) objectiveFcn(params0,paramsFixed,tspan,y0,data,simOptions);
[paramsEst,resnorm,residual,exitflag,output,lambda,Jacobian]=lsqnonlin(objFcn,params0,lowerbnd,upperbnd,optimOptions);
disp(['Algorithm: ',output.algorithm]);
disp(['Exit flag: ',num2str(exitflag),', resnorm = ',num2str(resnorm)]);
paramsVariance = resnorm*inv(Jacobian'*Jacobian)/length(data(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
paramsStd = sqrt(diag(paramsVariance));   % standard deviations for parameter estimates: 1.96*std is 95% conf interval
output = [paramsEst;paramsStd']

%% Run test of parameter estimation
params = [L,paramsEst(1),koff];
[tSingle,ySingle] = ode23(@LRodeFunc,tspan,y0,simOptions,params);
LRsingle = ySingle(:,2);
figure(1);
subplot(1,3,1);
plot(tSingle,LRsingle,data(:,1),data(:,2),'o');
xlabel('Time (sec)'); ylabel('LR (\muM)');
legend('model','data'); title('Least Squares Optimization');
subplot(1,3,2);
b = barweb([1 1],[0 1.96*paramsStd],[],'kon')
legend('Original','Estimated'); ylabel('Estimated parameters');

%% Plot entire SSE surface to verify that paramEst is in fact a minimum
konEst = paramsEst(1);
SSEest = sum(objFcn(konEst).^2);
konRange = [0.5:.01:2];
for runNumber = 1:length(konRange)
    SSE(runNumber) = sum(objFcn(konRange(runNumber)).^2);
end
subplot(1,3,3);
plot(konRange,SSE,konEst,SSEest,'*');
xlabel('kon'); ylabel('SSE');

%% Brute force plot
konRange = [0.5:.01:2];
for runNumber = 1:length(konRange)
    SSE(runNumber) = sum(objFcn(konRange(runNumber)).^2);
end
konRangeBruteForce = [0.5:0.1:2];
SSEbruteForce = interp1(konRange,SSE,konRangeBruteForce);
figure(2);
subplot(1,3,1);
plot(data(:,1),data(:,2),'o'); xlabel('Time (sec)'); ylabel('LR (\muM)');
subplot(1,3,2);
plot(konRange,SSE,konRangeBruteForce,SSEbruteForce,'o');
xlabel('kon'); ylabel('SSE'); title('Brute Force Estimation');
for i=1:length(konRangeBruteForce),
    text(konRangeBruteForce(i)-0.01,SSEbruteForce(i)+2,num2str(i),...
        'FontSize',8);
end
% Run test of parameter estimation
params = [L,1,koff];
[tSingle,ySingle] = ode23(@LRodeFunc,tspan,y0,simOptions,params);
LRsingle = ySingle(:,2);
subplot(1,3,3);
plot(tSingle,LRsingle,data(:,1),data(:,2),'o');
xlabel('Time (sec)'); ylabel('LR (\muM)');
legend('model','data');


    

function error = objectiveFcn(paramsEst,paramsFixed,tspan,y0,data,simOptions)
% objective function for performing least squares minimization
% This error function is weighted by the std dev of the experimental
% measurements, which affects confidence intervals and relative weighting 
% of multi-objective function estimates.
L = paramsFixed(1);
kon = paramsEst(1);
koff = paramsFixed(2);
params = [L,kon,koff];

[tSim,ySim] = ode23(@LRodeFunc,tspan,y0,simOptions,params);   % run the model for a given sample of paramsEst
ySimInterp = interp1(tSim,ySim,data(:,1)); % resample simulated data before performing subtraction
stdevExp = .1; % standard deviation of experimental data used in weights (set to 1 by default)
error = (1/stdevExp).*(ySimInterp(:,2) - data(:,end));     % note that LSQNONLIN computes the sum of squares for this automatically
