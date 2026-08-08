% This function estimates the volumetric flow rate V in cubic cm per sec
%Inputs: T0 (temperature in celsius), p0 (pressure in kPa)
%Inputs: dp (diameter of pipe in cm), Lp (length of pipe in cm)
function [V] = pvcInputs(T0,p0,dp,Lp)
T = T0 + 273.15;         % T = absolute temperature in K
p = p0*1000;             % p = outlet pressure in Pa
a = (dp/2)/100;          % a = radius of the pipe in m
L = Lp/100;              % L = pipe length in m
n = 0.00875*T - 3.51;    % n = power law index
Q = 146000;              % Q = activation energy in J/mol
R = 8.314;               % R = gas constant in J/mol*K
eta0 = (2.5*10^-13)*exp(Q/(R*T));  % eta0 = pre-exponential factor for power law viscosity function
V = ((pi*n)/(3*n+1))*((p/(2*eta0*L))^(1/n))*(a^((3*n+1)/n));  % solves for V in cubic m per sec
V = V*1000000;           % solves for V in cubic cm per sec
end

% Call a test case using this format: V = pvcInputs(T0,p0,dp,Lp)
% If another variable missing, guess & check until V is correct
V = pvcInputs()


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":22.9}
%---
