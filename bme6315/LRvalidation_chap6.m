%% Validation of algebraic ligand/receptor model
% model predictions
L = 10.^[-2:0.1:2];
Kd = 1.8;
fracLR = L./(Kd+L);
semilogx(L,fracLR);
xlabel('[Ligand] (\muM)');
ylabel('Fraction of Receptors Bound');

% experimental data
Lexp = 10.^[-2:1:2];
fracLRexp = [0.0327,0.067,0.478,0.892,0.931];
stdError = 0.10*[0.02 0.03 0.08 0.07 0.05];
hold on;
errorbar(Lexp,fracLRexp,stdError,'o')
axis([10^-2.5 10^2.5 0 1]);
hold off;

% calculate sum of squares error (SSE) and standard error of estimate (SEE)
n = 4; % repeats at each data point
k = 0; % number of fitted parameters
stdDev = stdError*sqrt(n); % 
fracLRinterp = interp1(log(L),fracLR,log(Lexp)); % interpolate simulation results where you have data
SSE = sum( 1./stdDev.^2 .* (fracLRexp - fracLRinterp).^2 )
SEE = sqrt(SSE/(n-k))