% LRparamEst3.m
% fitting parameters to two datasets

function LRparamEst3
%% Parameter assignments
L1 = 1;
L2 = 0.05;
kon = .5;  % initial guess for parameter value
koff = 0.1;
% params = [L,kon,koff];

% Simulation parameters
y0 = [1; 0];    % initial conditions
tspan = [0 6];
simOptions = [];

% Load experimental data
data = [0:6; ... % 1st column: timepoints
    0,0.58,0.79,0.93,0.88,0.94,0.80;... % 2nd column: data from L = 1
    0,0.12,0.24,0.25,0.33,0.35,0.33;]';  % 3rd column: data from L = 0.1

paramsFixed = [L1,L2,koff];     % designate parameters that will not be estimated
lowerbnd = [];              % optional lower bound on parameter estimate
upperbnd = [];              % optional upper bound on parameter estimate
params0 = [kon];            % initial guesses of parameter estimate

optimOptions = optimset('Display','iter');
objFcn = @(params0) objectiveFcn(params0,paramsFixed,tspan,y0,data,simOptions);
[paramsEst,resnorm,residual,exitflag,output,lambda,Jacobian]=lsqnonlin(objFcn,params0,lowerbnd,upperbnd,optimOptions);
disp(['Algorithm: ',output.algorithm]);
disp(['Exit flag: ',num2str(exitflag),', resnorm = ',num2str(resnorm)]);
paramsVariance = resnorm*inv(Jacobian'*Jacobian)/length(data(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
paramsStd = sqrt(diag(paramsVariance));   % standard deviations for parameter estimates
paramsCI = 1.96*paramsStd;                % mean +/- 1.96*stdev in 95% confidence interval
output = [paramsEst;paramsCI']

%% Run test of paremeter estimation
params = [L1,paramsEst(1),koff];
[tSingle,ySingle] = ode23(@LRodeFunc,tspan,y0,simOptions,params);
LRsingle = ySingle(:,2);
params2 = [L2,paramsEst(1),koff];
[tSingle2,ySingle2] = ode23(@LRodeFunc,tspan,y0,simOptions,params2);
LRsingle2 = ySingle2(:,2);

plot(tSingle,LRsingle,data(:,1),data(:,2),'o',...
    tSingle2,LRsingle2,data(:,1),data(:,3),'o');
xlabel('Time (sec)'); ylabel('LR (\muM)');
legend('model 1','data 1','model 2','data 2');

% %% Plot entire SSE surface
% konEst = paramsEst(1);
% SSEest = sum(objFcn(konEst).^2);
% konRange = [0.5:.01:2];
% for runNumber = 1:length(konRange)
%     SSE(runNumber) = sum(objFcn(konRange(runNumber)).^2);
% end
% subplot(1,2,2);
% plot(konRange,SSE,konEst,SSEest,'*');
% xlabel('kon'); ylabel('SSE');
    

function error = objectiveFcn(paramsEst,paramsFixed,tspan,y0,data,simOptions)
% objective function for performing least squares minimization
% This error function is weighted by the std dev of the experimental
% measurements, which affects confidence intervals and relative weighting 
% of multi-objective function estimates.
L1 = paramsFixed(1);
L2 = paramsFixed(2);
kon = paramsEst(1);
koff = paramsFixed(3);
tData = data(:,1);
LRdata1 = data(:,2);
LRdata2 = data(:,3);

% simulation for L = 1
params = [L1,kon,koff];
[tSim,ySim] = ode23(@LRodeFunc,tspan,y0,simOptions,params);   % run the model for a given sample of paramsEst

ySimInterp = interp1(tSim,ySim,tData); % resample simulated data before performing subtraction
stdevExp = 1; % standard deviation of experimental data used in weights (set to 1 by default)
error1 = (1/stdevExp).*(ySimInterp(:,2) - LRdata1);     % note that LSQNONLIN computes the sum of squares for this automatically

% simulation for L = 0.1
params2 = [L2,kon,koff];
[tSim,ySim] = ode23(@LRodeFunc,tspan,y0,simOptions,params2);   % run the model for a given sample of paramsEst
ySimInterp = interp1(tSim,ySim,tData); % resample simulated data before performing subtraction
stdevExp = 1; % standard deviation of experimental data used in weights (set to 1 by default)
error2 = (1/stdevExp).*(ySimInterp(:,2) - LRdata2);     % note that LSQNONLIN computes the sum of squares for this automatically
error = [error1,error2];
% error = error1;
% error = error2;