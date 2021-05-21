function linearParamEstEnsemble
% estimating a parameter ensemble using the Fisher Information Matrix to 
% estimate covariance between parameters.
% 5/21/2021 - Jeff Saucerman
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
paramEnsemble = repmat(paramsEst,n,1) + randn(n,length(paramsEst))*chol(paramsVariance); % see randn documentation

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