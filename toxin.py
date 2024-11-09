import numpy as np
import matplotlib.pyplot as plt

def toxin_cycling(
    Vb_per_kg, Vd_per_kg, infant_weight, SA_per_kg, Vd_full_per_BSA, P, 
    Cb_0, num_cycles, fill_duration, dwell_duration, drain_duration, dt, plot=False
):
    # Calculated parameters
    SA = SA_per_kg * infant_weight  # (cm^2)
    Vb_max = Vb_per_kg * infant_weight  # Blood volume in mL
    Vd_max = Vd_per_kg * infant_weight  # Dialysis fluid volume in mL
    Vd_full = Vd_full_per_BSA * SA  # Maximum dialysate volume for full diffusion (mL)
    k = P * SA  # Diffusion rate constant (mL/min)
    fill_rate = Vd_max / fill_duration
    drain_rate = Vd_max / drain_duration

    # Total simulation time
    total_time = (fill_duration + dwell_duration + drain_duration) * num_cycles
    time = np.arange(0, total_time, dt)  # Time vector

    # Initialize
    Cb = [] 
    Cd = [] 
    Vd = [] 

    Cb_t = Cb_0 
    Cd_t = 0 
    Vd_t = 0 

    # Concentration update function
    def update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt):
        if Vd_t > 0:
            SA_coverage = Vd_t / Vd_full # percent of surface area covered by dialysate
            delta = k * SA_coverage * (Cb_t - Cd_t) * dt
        else:
            delta = 0
        Cb_new = Cb_t - (delta / Vb_max) # Assume blood volume is constant
        Cd_new = Cd_t + (delta / Vd_t)
        return Cb_new, Cd_new

    # Simulation loop for each cycle
    for cycle in range(num_cycles):
        # Fill phase
        for t in np.arange(0, fill_duration, dt):
            dVd = fill_rate * dt
            Vd_t = min(Vd_max, Vd_t + dVd)
            Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
            Cd_t = Cd_t * (Vd_t / (Vd_t + dVd)) # Adjust Cd for dilution

            Cb.append(Cb_t)
            Cd.append(Cd_t)
            Vd.append(Vd_t)

        # Dwell phase (Dialysate volume stays constant)
        for t in np.arange(0, dwell_duration, dt):
            Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
            Cb.append(Cb_t)
            Cd.append(Cd_t)
            Vd.append(Vd_t)

        # Drain phase (Dialysate volume decreases)
        for t in np.arange(0, drain_duration, dt):
            dVd = -drain_rate * dt
            Vd_t = max(0.01, Vd_t + dVd)
            Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
            Cb.append(Cb_t)
            Cd.append(Cd_t)
            Vd.append(Vd_t)

    # Generate a new time vector matching the size of Cb, Cd, and Vd
    time_b = np.linspace(0, total_time, len(Cb))
    time_d = np.linspace(0, total_time, len(Vd))

    print(f"Cycle = {Cb[-1]}")
    
    if plot:
        DP_ratio = []
        for i in range(len(Cb)):
            DP_ratio.append(Cd[i] / Cb[i])
        # Plotting the results
        plt.figure(figsize=(12, 8))

        # Blood concentration plot
        plt.subplot(2, 1, 1)
        plt.plot(time_b, Cb, 'b', linewidth=1.5, label='Blood (C_b)')
        plt.plot(time_b, Cd, 'r--', linewidth=1.5, label='Dialysis Fluid (C_d)')
        plt.plot(time_b, DP_ratio, 'g--', linewidth=1.5, label='D/P Ratio')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Toxin Concentration (arbitrary units)')
        plt.legend()
        plt.title(f'Peritoneal Dialysis Toxin Clearance in an Infant: Final = {Cb[-1]}')

        # Dialysate volume plot
        plt.subplot(2, 1, 2)
        plt.plot(time_d, Vd, 'g', linewidth=1.5)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Dialysate Volume (mL)')
        plt.title('Dialysate Volume (V_d)')
        plt.grid(True)

        plt.tight_layout()
        plt.show()


def toxin_continuous(
    Vb_per_kg, Vd_per_kg, infant_weight, SA_per_kg, Vd_full_per_BSA, P, 
    Cb_0, fill_duration, drain_rate, dt, total_time, plot=False
):
    # Calculated parameters
    SA = SA_per_kg * infant_weight  # (cm^2)
    Vb_max = Vb_per_kg * infant_weight  # Blood volume in mL
    Vd_max = Vd_per_kg * infant_weight  # Dialysis fluid volume in mL
    Vd_full = Vd_full_per_BSA * SA  # Maximum dialysate volume for full diffusion (mL/)
    k = P * SA  # Diffusion rate constant (mL/min)
    fill_rate = Vd_max / fill_duration
    drain_duration = Vd_max / drain_rate
    dwell_duration = total_time - fill_duration - drain_duration

    # Total simulation time
    time = np.arange(0, total_time, dt)  # Time vector

    # Initialize
    Cb = [] 
    Cd = [] 
    Vd = [] 

    Cb_t = Cb_0 
    Cd_t = 0 
    Vd_t = 0 

    # Concentration update function
    def update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt):
        if Vd_t > 0:
            SA_coverage = Vd_t / Vd_full # percent of surface area covered by dialysate
            delta = k * SA_coverage * (Cb_t - Cd_t) * dt
        else:
            delta = 0
        Cb_new = Cb_t - (delta / Vb_max) # Assume blood volume is constant
        Cd_new = Cd_t + (delta / Vd_t)
        return Cb_new, Cd_new


    for t in np.arange(0, fill_duration, dt):
        dVd = fill_rate * dt
        Vd_t = min(Vd_max, Vd_t + dVd)
        Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
        Cd_t = Cd_t * (Vd_t - dVd) / Vd_t # Adjust Cd for dilution

        Cb.append(Cb_t)
        Cd.append(Cd_t)
        Vd.append(Vd_t)

    # Continuous cycling phase
    for t in np.arange(0, dwell_duration, dt):
        Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
        Cd_t = Cd_t * (Vd_t - drain_rate * dt) / Vd_t # Adjust Cd for dilution

        Cb.append(Cb_t)
        Cd.append(Cd_t)
        Vd.append(Vd_t)

    for t in np.arange(0, drain_duration, dt):
        dVd = -drain_rate * dt
        Vd_t = max(0.01, Vd_t + dVd)
        Cb_t, Cd_t = update_concentration(Cb_t, Cd_t, Vd_t, Vd_max, k, dt)
        Cb.append(Cb_t)
        Cd.append(Cd_t)
        Vd.append(Vd_t)



    # Generate a new time vector matching the size of Cb, Cd, and Vd
    time_b = np.linspace(0, total_time, len(Cb))
    time_d = np.linspace(0, total_time, len(Vd))



    print(f"Cont. = {Cb[-1]}")

    if plot:
        DP_ratio = []
        for i in range(len(Cb)):
            DP_ratio.append(Cd[i] / Cb[i])
        # Plotting the results
        plt.figure(figsize=(12, 8))

        # Blood concentration plot
        plt.subplot(2, 1, 1)
        plt.plot(time_b, Cb, 'b', linewidth=1.5, label='Blood (C_b)')
        plt.plot(time_b, Cd, 'r--', linewidth=1.5, label='Dialysis Fluid (C_d)')
        plt.plot(time_b, DP_ratio, 'g--', linewidth=1.5, label='D/P Ratio')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Toxin Concentration (arbitrary units)')
        plt.legend()
        plt.title(f'Peritoneal Dialysis Toxin Clearance in an Infant: Final = {Cb[-1]}')

        # Dialysate volume plot
        plt.subplot(2, 1, 2)
        plt.plot(time_d, Vd, 'g', linewidth=1.5)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Dialysate Volume (mL)')
        plt.title('Dialysate Volume (V_d)')
        plt.grid(True)

        plt.tight_layout()
        plt.show()


if __name__ == '__main__':

    plot = True

    # Parameters for Cycling

    Vb_per_kg = 80  # Blood volume per kg of infant (mL/kg)
    infant_weight_kg = 5  # Infant weight in kg
    SA_constant = 533  # Surface area constant (cm^2)
    Vd_full_per_BSA = 800 / 10000  # Maximum dialysate volume (mL/cm^2)
    Cb_0 = 1  # Initial blood toxin concentration (mg/mL)
    P = 0.0005  # Permeability constant (mL/min/cm^2) #calculated from 60 min D/P ratio for creatinine from Pediatric Dialysis textbook page 200
    dt = 0.1  # Time step for simulation (minutes)

    Vd_per_kg = 30  # Dialysate volume per kg of infant (mL/kg)
    num_cycles = 10  # Number of dialysis cycles
    fill_duration = 5  # Fill phase duration (minutes)
    dwell_duration = 45  # Dwell phase duration (minutes)
    drain_duration = 10  # Drain phase duration (minutes)

    toxin_cycling(Vb_per_kg, Vd_per_kg, infant_weight_kg, SA_constant, Vd_full_per_BSA, P, Cb_0, num_cycles, fill_duration, dwell_duration, drain_duration, dt, plot=plot)

    # Parameters for Continuous 

    Vd_per_kg = 23  # Dialysate volume per kg of infant (mL/kg)
    fill_duration = 5  # Fill phase duration (minutes)
    drain_rate = 20  # Drain rate (mL/min)
    total_time = num_cycles * 60  # Total simulation time (minutes)

    toxin_continuous(Vb_per_kg, Vd_per_kg, infant_weight_kg, SA_constant, Vd_full_per_BSA, P, Cb_0, fill_duration, drain_rate, dt, total_time, plot=plot)


