% This code models the thermal degradation of PVC

% Load temperature-time history data from CSV files
temp1 = table2array(PVC_time_temp1_csv);
temp2 = table2array(PVC_time_temp2_csv);

% Extract time (s) and temperature (°C) columns for each dataset
t1 = temp1(:,1);
T1 = temp1(:,2);
t2 = temp2(:,1);
T2 = temp2(:,2);

% Arrhenius/reaction parameters
A = 6.6e15;
E_A = 190000;
R = 8.314;
n = 1.6;

% Total cooling duration: 5 minutes = 300 seconds
t_end = 300;

% Time-step sizes (in secs) to demonstrate solution convergence.
dt_values = [10,5,1,0.1];

line_styles = {'-', '--', ':', '-.'};

% CASE STUDY 1
figure(1);
hold on;
title('Case Study 1: Mass Fraction \alpha vs Time');
xlabel('Time (s)');
ylabel('Mass Fraction \alpha = m(t)/m_0');
grid on;

%Set up for loop to iteratively compute delta alpha and alpha
for i = 1:length(dt_values)
    dt = dt_values(i);         % Current time-step size (s)
    t = 0:dt:t_end;            % Time array from 0 to 300 s
    N = length(t);             % Number of time steps

    alpha = zeros(1,N);        % Pre-allocate mass fraction array
    alpha(1) = 1;              % Initial condition: α(t=0) = 1

    % Inner loop: iterates over each time step to compute α at the next time
    for j = 1:N-1

        % Interpolate temperature at current time t(j) from the CSV data
        T_current = interp1(t1, T1, t(j), 'linear');   

        % Convert temperature from Celsius to Kelvin
        T_K = T_current + 273.15;

        % Compute temperature-dependent rate constant k using Arrhenius equation
        k = A * exp(-E_A / (R * T_K));
        
        % Compute the reaction rate dα/dt = -k(T) * α^n
        d_alpha_dt = -k * alpha(j)^n;
        
        % Numerical time integration; assume that dα/dt is constant over the small interval Δt
        alpha(j+1) = alpha(j) + d_alpha_dt * dt;

    end

% Plot results for this dt value
plot(t, alpha,line_styles{i}, 'LineWidth', 1.5,'DisplayName', ['\Deltat = ' num2str(dt) ' s']);

end

legend('Location', 'southwest');

% CASE STUDY 2 (rinse and repeat Case 1 using PVC_time_temp2_csv instead)

figure(2);
hold on;
title('Case Study 2: Mass Fraction \alpha vs Time');
xlabel('Time (s)');
ylabel('Mass Fraction \alpha = m(t)/m_0');
grid on;

for i = 1:length(dt_values)
    dt = dt_values(i);
    t = 0:dt:t_end;
    N = length(t);

    alpha = zeros(1,N);
    alpha(1) = 1;

    for j = 1:N-1
        T_current = interp1(t2, T2, t(j), 'linear');
        T_K = T_current + 273.15;
        k = A * exp(-E_A / (R * T_K));

        d_alpha_dt = -k * alpha(j)^n;

        alpha(j+1) = alpha(j) + d_alpha_dt * dt;

        
    end

plot(t, alpha,line_styles{i}, 'LineWidth', 1.5,'DisplayName', ['\Deltat = ' num2str(dt) ' s']);
end

legend('Location', 'southwest');