% Realistic Parameters for Infant
blood_vol = 80; % ml/kg
dialysate_vol = 20; % ml/kg
infant_weight = 5; % kg

V_b_max = blood_vol * infant_weight;             % Blood volume in mL
V_d_max = dialysate_vol * infant_weight;         % Dialysis fluid volume in mL
k = 0.01;                                        % Diffusion rate constant (mL/min)
C_b0 = 1;                                      % Initial blood toxin concentration (arbitrary units)
num_cycles = 8;                                 % Number of dialysis cycles
fill_duration = 5;                               % Fill phase duration (minutes)
dwell_duration = 45;                             % Dwell phase duration (minutes)
drain_duration = 10;                             % Drain phase duration (minutes)

% Time settings
dt = 0.01;                                        % Time step (minutes)
total_time = (fill_duration + dwell_duration + drain_duration) * num_cycles;
time = 0:dt:total_time;                          % Total time vector

% Initialize concentrations and dialysate volume
C_b = [];                                        % Blood concentration over time
C_d = [];                                        % Dialysis fluid concentration over time
V_d = [];                                        % Dialysate volume over time

% Initial conditions
C_b_current = C_b0;                              % Initial blood toxin concentration
C_d_current = 0;                                 % Initial dialysis fluid toxin concentration
V_d_current = 0;                                 % Initial dialysate volume (at 0)

% Concentration update function
update_concentration = @(C_b_current, C_d_current, V_d_current, V_d_max) ...
    deal( ...
        C_b_current - k * (V_d_current / V_d_max) * (C_b_current - C_d_current) * dt, ...
        C_d_current + (V_b_max / V_d_max) * k * (C_b_current - C_d_current) * (V_d_current / V_d_max) * dt ...
    );

% Simulation loop for each cycle
for cycle = 1:num_cycles
    % Fill phase
    for t = 0:dt:fill_duration
        [C_b_current, C_d_current] = update_concentration(C_b_current, C_d_current, V_d_current, V_d_max);
        C_b = [C_b, C_b_current];
        C_d = [C_d, C_d_current];
        V_d_current = min(V_d_max, V_d_current + (V_d_max / fill_duration) * dt);  % Increment until max volume
        V_d = [V_d, V_d_current];  % Store dialysate volume
    end
    
    % Dwell phase (Dialysate volume stays constant)
    for t = 0:dt:dwell_duration
        [C_b_current, C_d_current] = update_concentration(C_b_current, C_d_current, V_d_current, V_d_max);
        C_b = [C_b, C_b_current];
        C_d = [C_d, C_d_current];
        V_d = [V_d, V_d_current];  % Dialysate volume constant
    end
    
    % Drain phase (Dialysate volume decreases)
    for t = 0:dt:drain_duration
        [C_b_current, C_d_current] = update_concentration(C_b_current, C_d_current, V_d_current, V_d_max);
        C_b = [C_b, C_b_current];
        C_d = [C_d, C_d_current];
        V_d_current = max(0.01, V_d_current - (V_d_max / drain_duration) * dt);  % Decrease dialysate volume
        V_d = [V_d, V_d_current];  % Store dialysate volume
    end
    
    % Reset dialysis fluid concentration for next cycle
    C_d_current = 0;
end

% Generate a new time vector matching the size of C_b, C_d, and V_d_current
time_b = linspace(0, total_time, length(C_b)); % Ensure time matches C_b's length
time_d = linspace(0, total_time, length(V_d)); % Ensure time matches V_d's length

% Plotting the results
figure;

% Blood concentration plot
subplot(2,1,1);
plot(time_b, C_b, 'b', 'LineWidth', 1.5); hold on;
plot(time_b, C_d, 'r--', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Toxin Concentration (arbitrary units)');
legend('Blood (C_b)', 'Dialysis Fluid (C_d)');
title('Peritoneal Dialysis Toxin Clearance in an Infant');

% Dialysate volume plot
subplot(2,1,2);
plot(time_d, V_d, 'g', 'LineWidth', 1.5);
xlabel('Time (minutes)');
ylabel('Dialysate Volume (mL)');
title('Dialysate Volume (V_d)');
grid on;
