function paramEstEnsembleImportanceMCMC
% estimating a parameter ensemble using importance sampling MCMC.
% Which can estimate nonlinear relationships between parameters.
% 5/21/2021 - Jeff Saucerman, based on 'toyModelParamEstEnsemble.m', which
% I wrote in ~2012.
% Note: This code seems to work, but its predictions are worse than the
% covariance approach in paramEstEnsembleCovariance.
%
% Modeling the function y = a*x + b
% simulated data

xdata = [0:0.1:1]';
ydata = [0:0.1:1]' + 0.2*randn(length(xdata),1);
 
% define the model

params0 = [1, 1];   
fun = @(params) params(1)*xdata + params(2); % f(x) = a*x+b
% fun = @(params) params(1)*params(2)*xdata; % f(x) = a*b*x; (singular matrix, undefined variance matrix)

% perform parameter estimation and plot the fit to data

objFun = @(params) fun(params) - ydata;
[paramsEst,resnorm,~,~,~,~,Jacobian]=lsqnonlin(objFun,params0);

subplot(2,2,1);
plot(xdata,fun(paramsEst),'-',xdata,ydata,'o')
title('fit');
legend('model','data');

% calculate confidence in parameter estimates

paramsVariance = resnorm*inv(Jacobian'*Jacobian)/length(ydata(:,1));  % estimates the parameter estimate covariance using the Cramer-Rao inequality
paramsStd = sqrt(diag(paramsVariance));   % standard deviations for parameter estimates: 2*std is 95% conf interval

subplot(2,2,2);
bar([1,2],paramsEst);
title('parameters');
hold on
er = errorbar([1,2],paramsEst,paramsStd);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off

% generate parameter ensemble based on the covariance matrix
n = 10;
Temp = 0.1;
[paramEnsemble,accept] = calcParamEnsemble(objFun,paramsEst,Jacobian,n,Temp); % use importance-sampled MCMC
% [paramEnsemble,accept] = calcIndepParamEnsemble(paramsEst,paramsStd',n)        % directly from paramsStd


disp([num2str(size(paramEnsemble,1)),' param sets generated, acceptance rate: ',num2str(accept*100),'%']);

subplot(2,2,3);
plot(paramEnsemble(:,1),paramEnsemble(:,2),'o');
xlabel('param 1'); ylabel('param 2');
title('parameter ensemble');

subplot(2,2,4); title('ensemble predictions');
hold on;
for i = 1:n
    plot(xdata,fun(paramEnsemble(i,:)));
end
hold off;

function [paramEnsemble,accept] = calcParamEnsemble(objFun,paramsEst,Jacobian,nSamples,Temp)
    % estimating a parameter ensemble using importance-sampled MCMC
    % method for ensemble estimation from: Brown, Sethna Phys Rev E 2004
    % This uses the Jacobian to identify sloppy directions in parameter
    % space.

    SSE = @(params) sum(objFun(params).^2);
    mylogpdf = @(params) -0.5*SSE(params)/Temp;     % use logpdf
    
    Hessian = full(Jacobian'*Jacobian); % Gauss-Newton approximation of the Hessian
    [U,S,V] = svd(0.5*Hessian); % X = U*S*V', 
    Scutoff = 0.001*max(diag(S)); % may need to adjust this cutoff
    D = max(S,Scutoff*eye(length(S)))^(-1); % dimensions length(S)*length(S)
    samplingMatrix = V'*sqrt(D)/sqrt(length(S));
    proprnd = @(x) x+randn(1,length(S))*samplingMatrix;
    
    [paramEnsemble,accept] = mhsample(paramsEst,nSamples,'logpdf',mylogpdf,'proprnd',proprnd,'symmetric',1);   
    
function [paramEnsemble,accept] = calcIndepParamEnsemble(paramsEst,paramsStd,nSamples)
% Generate parameter sets using a normal distribution based on paramsStd.
% This approach assumes parameters are indepenent.
    accept = 1;
    for i=1:nSamples
        paramEnsemble(i,:) = paramsEst + paramsStd.*randn(size(paramsEst));
    end

