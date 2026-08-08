% CHILLER SYSTEM DESIGN CODE
% Usage: Call chillerInputs() with the design variables listed below. See function definition for input/output details.
% Load water properties from the imported CSV table
% Columns: T(K), rho(kg/m^3), eta(Pa·s), kf(W/m/K), Pr(non-dimensional)

waterprops1 = table2array(waterprops); % Convert table to numeric array
T = waterprops1(:,1);                  % Water temperature [K] 
rho = waterprops1(:,2);                % Water density [kg/m^3]
eta = waterprops1(:,3);                % Water viscosity [Pa·s]
kf = waterprops1(:,4);                 % Water thermal conductivity [W/m/K]
Pr = waterprops1(:,5);                 % Prandtl number [non-dimensional]


% The function below computes the time to cool PVC to a target temperature... 
% inside a cylindrical mold with forced water convection.

% INPUTS:
%   T_mPVC  — Initial PVC temperature [°C]
%   T_fPVC  — Final (target) PVC temperature [°C]
%   T_w     — Cooling water temperature [°C]
%   v_w     — Water velocity over mold [cm/s]
%   d_m     — Mold diameter [cm]
%   T       — Water property table: temperature [K]
%   rho     — Water property table: density [kg/m^3]
%   eta     — Water property table: dynamic viscosity [Pa·s]
%   kf      — Water property table: thermal conductivity [W/m/K]
%   Pr      — Water property table: Prandtl number [non-dimensional]

% OUTPUT:
%   t_f     — Time to cool PVC from T_mPVC to T_fPVC [s]

function [t_f] = chillerInputs(T_mPVC,T_fPVC,T_w,v_w,d_m,T,rho,eta,kf,Pr)

% Unit conversions — convert all inputs to SI units

T_mPVC_K = T_mPVC + 273.15;   % Initial PVC temperature: °C -> K
T_w_K = T_w + 273.15;         % Water temperature: °C -> K
d_m_SI = d_m / 100;           % Mold diameter: cm -> m
v_w_SI = v_w / 100;           % Water velocity: cm/s -> m/s


T_f = (T_mPVC_K + T_w_K)/2;   % Film temperature [K]

% Interpolate water properties at the film temperature

rho_f = interp1(T, rho, T_f, 'linear');    
eta_f = interp1(T, eta, T_f, 'linear');
kf_f  = interp1(T, kf,  T_f, 'linear');
Pr_f  = interp1(T, Pr,  T_f, 'linear');

Re_d = (rho_f*v_w_SI*d_m_SI)/eta_f;   % Reynolds number

% Get constants C & n from Reynolds number based on table in the assignment

if Re_d >= 0.4 && Re_d < 4    
    C = 0.989;
    n = 0.330;
elseif Re_d < 40
    C = 0.911;
    n = 0.385;
elseif Re_d < 4000
    C = 0.683;
    n = 0.466;
elseif Re_d < 40000
    C = 0.193;
    n = 0.618;
elseif Re_d <= 400000
    C = 0.027;
    n = 0.805;
end

Nu_d = C*((Re_d)^n)*((Pr_f)^(1/3));   % Average Nusselt number

h = (Nu_d*kf_f)/d_m_SI;       % Convective heat transfer coefficient
                              % Derived from Nu = h*d/kf -> h = Nu*kf/d

alpha = 0.11e-6; % Thermal diffusivity of PVC [m^2/s] (converted from 0.11 mm^2/s)
k = 0.2;         % Thermal conductivity of PVC [W/m/K]

% Cooling time from lumped capacitance model
% Solving for t when T(t) = T_fPVC:
% t_f = -(k*dm)/(4*h*alpha) * ln((T_fPVC - T_w)/(T_mPVC - T_w))

t_f = -((k*d_m_SI)/(4*h*alpha))*log((T_fPVC - T_w)/(T_mPVC - T_w));

end

% Call a test case using this format: t_f = chillerInputs(T_mPVC,T_fPVC,T_w,v_w,d_m,T,rho,eta,kf,Pr)
% If another variable missing, guess & check until t_f is correct
t_f = chillerInputs(T_mPVC,T_fPVC,T_w,v_w,d_m,T,rho,eta,kf,Pr)