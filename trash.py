import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define parameters
MTAC = 20
V1 = 800
epsilon = 1e-6

# Initial conditions
C1_0 = 10
C2_0 = 0
y0 = [C1_0, C2_0]

def V2(t):
    t_mod = t % 60  # For repeating every 60 minutes
    if t_mod < 10:
        return 30 * t_mod
    elif 10 <= t_mod < 40:
        return 300
    elif 40 <= t_mod < 60:
        return 300 - 15 * (t_mod - 40)

def dV2(t):
    t_mod = t % 60  # For repeating every 60 minutes
    if t_mod < 10:
        return 30
    elif 10 <= t_mod < 40:
        return 0
    elif 40 <= t_mod < 60:
        return -15

# Define the system of differential equations
def equations(t, y):
    C1, C2 = y
    current_V2 = V2(t)
    
    # Reset C2 to 0 at the start of each 60-minute cycle
    if abs(t % 60) < 1e-3:  # Use a small tolerance to approximate multiples of 60
        C2 = 0
    
    dC1_dt = -MTAC / V1 * (C1 - C2)
    dC2_dt = MTAC / (current_V2 + epsilon) * (C1 - C2) - dV2(t) / (current_V2 + epsilon) * C2
    
    print(f"t={t}, C1={C1}, C2={C2}, V2={current_V2}, dC1_dt={dC1_dt}, dC2_dt={dC2_dt}")
    return [dC1_dt, dC2_dt]

# Time range for the solution
t_span = (0, 120)  # Solve over a range of 120 minutes to see multiple cycles
t_eval = np.linspace(0, 120, 200)

# Solve the differential equations
sol = solve_ivp(equations, t_span, y0, t_eval=t_eval)

# Plot the results
plt.figure(figsize=(10, 5))
plt.plot(sol.t, sol.y[0], label="C1 (Blood)")
plt.plot(sol.t, sol.y[1], label="C2 (Dialysate)")
plt.xlabel("Time (minutes)")
plt.ylabel("Concentration")
plt.legend()
plt.title("Concentration vs. Time")
plt.show()
