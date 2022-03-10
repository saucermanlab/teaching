function [Y,S] = toyModelSensitivity
% various sensitivity analyses using the toy model
% 09/13/2011 - Jeff Saucerman
% updated 20181009 by Jeff Saucerman
% Modeling the following reactions:
% 1) L+R<->LR, 2) LR+E<->LRE, 3) LRE catalyzes S->P, 4) P->S
% runs sensitivity analyses
%% define parameters
k1f = 1;    % [uM^-1 s^-1] react1 forward rate constant
k1r = 1;    % [s^-1] react1 reverse rate constant
k2f = 1;    % [uM^-1 s^-1] react2 forward rate constant
k2r = 1;    % [s^-1] react2 reverse rate constant
kcat = 1;   % [s^-1] catalytic rate constant for enzyme
Km = 1;     % [uM] Michaelis constant for enzyme
k4 = 1;     % [s^-1] react4 rate constant
k5 = 1;     % [s^-1] react5 rate constant
Vratio = 0.1; % dimensionless
L = 1;   % [uM] concentration of ligand
Rtot = 1;   % [uM] total concentration of receptor
Etot = 1;   % [uM] total concentration of enzyme
Stot = 1;   % [uM] total concentration of substrate
params = {k1f,k1r,k2f,k2r,kcat,Km,k4,k5,Vratio,L,Rtot,Etot,Stot};

%% Sensitivity analysis coefficients and response curves
% Calculates log sensitivity coefficients for output Y, with parameter increments
% specified in deltaParams for parameters listed in paramList
y0 = [0; 0; 0; 0];   % or you could load saved data like: y0 = load('yfinal.dat');
tspan = [0 100];
options = [];
[t,y] = ode23(@toyModelODEfunc2,tspan,y0,options,params); % Run model with default params
Ydefault = y(end,end);  % the output value that you care about

deltaParams = [0.1,1.1,10];
% deltaParams = 10.^[-1:.1:2];    % more refined range for response curves
for deltaNum=1:length(deltaParams)
    for paramNum = 1:length(params)    
        paramsNew = params;
        paramsNew{paramNum} = deltaParams(deltaNum)*params{paramNum};
        try
            [t,y] = ode23(@toyModelODEfunc2,tspan,y0,options,paramsNew);
            Y(paramNum,deltaNum) = y(end,end);
            S(paramNum,deltaNum) = (Y(paramNum,deltaNum)-Ydefault)./(paramsNew{paramNum}-params{paramNum})*params{paramNum}/Ydefault;            
        catch       % in case the tested parameters generate numerical problems
            Y(paramNum,deltaNum) = NaN;
            S(paramNum,deltaNum) = NaN;
            disp(['Error at paramNum = ',num2str(paramNum),', delta = ',num2str(deltaParams(deltaNum))]);
        end
    end
end

% Plotting sensitivity coefficients
figure(1);
subplot(1,2,1); plot(t,y(:,3)); xlabel('Time (sec)'); ylabel('P (\muM)'); 
subplot(1,2,2); bar(S); xlabel('Parameter number'); ylabel('S=dY/dp*p/Y for Product'); legend('0.1x','1.1x','10x');

% Plot response curves
figure(2);
size = ceil(sqrt(length(params)));
paramNames = {'k1f','k1r','k2f','k2r','kcat','Km','k4','k5','Vratio','Ltot','Rtot','Etot','Stot'};
for i=1:length(params)
    subplot(size,size,i);
    semilogx(deltaParams,Y(i,:));
    axis([min(deltaParams) max(deltaParams) 0 max(max(Y))]);
    xlabel(paramNames{i});
end

%% Uncertainty analysis
% Example of how to generate new randomized parameter values:
% >>x = [1,2,3,4,5];
% >>cv = 0.1*ones(size(x));
% >>x+x.*cv.*randn(1,length(x))
cv = 0.3*ones(1,length(params));  % coefficient of variation for each parameter: cv = stdev/mean
 y0 = [0; 0; 0; 0];   % or you could load saved data like: y0 = load('yfinal.dat');
tspan = [0 100];
tinterp = [0:1:100];
options = [];
figure(3); subplot(2,1,1); 
for run=1:10
    try
        paramsNew = num2cell([params{:}] + cv.*[params{:}].*randn(1,length(params)));               % randomize new params with cv
        [t,y] = ode23(@toyModelODEfunc2,tspan,y0,options,paramsNew);
        P = y(:,end);   % entire timecourse of the last state variable, P
        Pinterp(:,run) = interp1(t,P,tinterp);  % holds matrix of all P timecourses
        plot(tinterp,Pinterp(:,run)); hold on
    catch       % in case the tested parameters generate numerical problems
        disp(['Error with params = ',num2str(paramsNew)]);
    end
end
xlabel('time'); ylabel('P');
Pmean = mean(Pinterp,2);
Pstd = std(Pinterp,0,2);
subplot(2,1,2);
% length(tinterp)
% length(Pmean)
% length(Pstd)
plot(tinterp,Pmean,'k-',tinterp,Pmean+2.*Pstd,'b--',tinterp,Pmean-2.*Pstd,'b--');
xlabel('time'); ylabel('P');

%% 2D sensitivity matrix
y0 = [0; 0; 0; 0];   % or you could load saved data like: y0 = load('yfinal.dat');
tspan = [0 100];
options = [];
[t,y] = ode23(@toyModelODEfunc2,tspan,y0,options,params); % run default simulation
y_fin_orig = y(end,:)';
for i=1:length(params)
    paramsNew = params; % Initialize new parameters to old ones
    y_initial = y0;
    paramsNew{i} = params{i}*1.1; % Increase parameter by 10%
    [~,y_new] = ode15s(@toyModelODEfunc2,tspan,y_initial,options,paramsNew);
    y_final(:,i) = y_new(end,:)'; % Make column of y_final equal to steady state values of y_new
    
    % Generate normalized sensitivity coefficients
    S2d(:,i) = (y_final(:,i) - y_fin_orig(:,1))./(paramsNew{i}-params{i}).*params{i}./y_fin_orig(:,1);
end
figure(4);

imagesc(S2d)
colorbar;
maxVal = max(max(S2d));
minVal = -maxVal;
caxis([-2, 2])
xlabel('Perturbed Parameter'); ylabel('Sensitivity of State Variable');
yNames = {'LR','LRE','P','Pnuc'};
paramNames = {'k1f','k1r','k2f','k2r','kcat','Km','k4','k5','Vratio','Ltot','Rtot','Etot','Stot'};
set(gca,'YTick',1:length(yNames));
set(gca,'YTickLabel',yNames);
set(gca,'XTick',1:length(paramNames));
set(gca,'XTickLabel',paramNames);