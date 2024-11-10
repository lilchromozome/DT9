% Parameters for Continuous Peritoneal Dialysis System (Enhanced Diffusion Model)
blood_vol = 80;       % mL/kg (Blood volume per kg of body weight)
dialysate_vol = 15;    % mL/kg (Dialysate volume per kg of body weight)
infant_weight = 5;    % kg (Infant weight)

V_b_max = blood_vol * infant_weight;   % Blood volume in mL
V_d_max = dialysate_vol * infant_weight; % Maximum dialysate volume in mL (constant)
P = 0.02;                               % Permeability coefficient (mL/min/cm²)
A = 1000;                               % Surface area (cm²)
C_b0 = 100;                             % Initial blood toxin concentration (arbitrary units)
drain_rate = 5;                         % Drain rate (mL/min)

% Time settings
dt = 0.01;                               % Time step (minutes)
total_time = 400;                       % Total simulation time (minutes)
time = 0:dt:total_time;                 % Total time vector

% Initial concentrations
C_b_current = C_b0;                     % Initial blood toxin concentration
C_d_current = 0;                        % Initial dialysis fluid toxin concentration
V_d_current = V_d_max;                  % Initial dialysate volume

% Initialize arrays to store concentrations and volumes
C_b = [];                               % Blood concentration over time
C_d = [];                               % Dialysis fluid concentration over time
V_d = [];                               % Dialysate volume over time

% Simulation loop (continuous fill and drain)
for t = time
    % Calculate the diffusion rate based on Fick's law
    diffusion_rate = (P * A) * (C_b_current - C_d_current);
    
    % Update concentrations due to diffusion
    dC_b = -diffusion_rate * dt / V_b_max;      % Change in blood concentration
    dC_d = diffusion_rate * dt / V_d_current;   % Change in dialysate concentration
    
    % Update concentrations
    C_b_current = C_b_current + dC_b;           % Update blood concentration
    C_d_current = (C_d_current + dC_d) * ((V_d_current - drain_rate * dt) / V_d_current); % Adjust for drainage
    
    % Store concentrations and volume at each time step
    C_b = [C_b, C_b_current];
    C_d = [C_d, C_d_current];
    V_d = [V_d, V_d_current];                   % Dialysate volume remains constant for this model
end

% Generate a time vector that matches the size of the concentration arrays
time_b = linspace(0, total_time, length(C_b)); % Time vector for blood concentration
time_d = linspace(0, total_time, length(C_d)); % Time vector for dialysate concentration

% Plotting the results
figure;

% Blood and Dialysate concentration plot (on the same plot)
subplot(3,1,1);
plot(time_b, C_b, 'b', 'LineWidth', 1.5); hold on;
plot(time_d, C_d, 'r', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Toxin Concentration (mg/mL)');
legend('Blood (C_b)', 'Dialysis Fluid (C_d)');
title('Blood and Dialysis Fluid Toxin Concentrations with Enhanced Diffusion Model');
grid on;

% Dialysate volume plot (constant volume)
subplot(3,1,2);
plot(time_d, V_d, 'g', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Dialysate Volume (mL)');
title('Constant Dialysate Volume (V_d)');
grid on;

% Initialize array to store the ratio of C_d to C_b
C_ratio = C_d ./ C_b;

% Generate a time vector matching the size of the C_ratio array
time_ratio = linspace(0, total_time, length(C_ratio));

% Plotting the ratio of C_d to C_b
subplot(3,1,3);
plot(time_ratio, C_ratio, 'm', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Ratio of C_d to C_b');
title('C_d / C_b');
grid on;
