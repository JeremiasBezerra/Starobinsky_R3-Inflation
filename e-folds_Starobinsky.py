# -*- coding: utf-8 -*-
"""
@author: J. Bezerra-Sobrinho
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# ----------------------------
# Constants
# ----------------------------
PI = np.pi
a = np.sqrt(16 * PI / 3)  # Model-dependent constant in the Starobinsky potential
V0 = 1.0                   # Overall normalization

# ----------------------------
# Initial conditions
# phi0 : initial field value
# pi0  : initial field velocity (pi = dphi/dt)
# N0   : initial number of e-folds
# ----------------------------
phi0, pi0, N0 = 1.54, 0.0, 0.0

# ----------------------------
# Potential V(phi) and its derivative dV/dphi
# ----------------------------
def V(phi):
    return V0 * (1 - np.exp(-a * phi))**2

def dV(phi):
    return 2 * a * V0 * np.exp(-a * phi) * (1 - np.exp(-a * phi))

# ----------------------------
# System of ODEs:
# y = (phi, pi, N)
#
# dphi/dt = pi
# dpi/dt  = -3H*pi - dV/dphi
# dN/dt   = H
#
# H is computed from the Friedmann equation:
#    H^2 = (pi^2/2 + V)/ (3/(8π))
# ----------------------------
def eom(t, y):
    phi, pi, N = y
    Hsq = (pi**2 / 2 + V(phi)) * 8 * PI / 3
    H = np.sqrt(max(Hsq, 0.0))  # Numerical safety
    return [pi, -3.0 * H * pi - dV(phi), H]

# ----------------------------
# Event function: inflation ends when epsilon = 1
# epsilon = 4π pi^2 / H^2
# ----------------------------
def end_of_inflation(t, y):
    phi, pi, _ = y
    Hsq = (pi**2 / 2 + V(phi)) * 8 * PI / 3
    H = np.sqrt(max(Hsq, 0.0))
    eps = 4 * PI * pi**2 / (H**2) if H > 0 else np.inf
    return eps - 1.0  # Inflation ends when this crosses zero

end_of_inflation.terminal = True   # Stop integration at epsilon = 1
end_of_inflation.direction = 0     # Detect crossing in any direction

# ----------------------------
# Integrate the system until inflation ends
# ----------------------------
sol = solve_ivp(
    eom,
    [0.0, 1e5],                   # Large upper bound; event will terminate earlier
    [phi0, pi0, N0],
    events=end_of_inflation,
    dense_output=True,
    max_step=0.05,
    rtol=1e-9,
    atol=1e-12
)

# Retrieve the event time (end of inflation)
if sol.t_events and sol.t_events[0].size > 0:
    t_end = float(sol.t_events[0][0])
else:
    t_end = float(sol.t[-1])

phi_end, pi_end, N_end = sol.sol(t_end)

# ----------------------------
# Sample the numerical solution on a uniform grid in time
# ----------------------------
n_samples = 5000
t_samples = np.linspace(0.0, t_end, n_samples)
phi_samples, N_samples = sol.sol(t_samples)[0], sol.sol(t_samples)[2]

# ----------------------------
# Sort by N and remove duplicates (needed for interpolation)
# ----------------------------
order = np.argsort(N_samples)
N_sorted, phi_sorted = N_samples[order], phi_samples[order]

tol = 1e-12
mask = np.concatenate(([True], np.abs(np.diff(N_sorted)) > tol))

N_u, phi_u = N_sorted[mask], phi_sorted[mask]

# ----------------------------
# Build cubic interpolation phi(N)
# ----------------------------
phi_of_N = interp1d(N_u, phi_u, kind='cubic', bounds_error=True)

# ----------------------------
# Function phi_k(N): field value when N e-folds remain
# We evaluate phi at (N_end - Nk).
# ----------------------------
def phik(Nk):
    Nk_arr = np.asarray(Nk, dtype=float)
    if np.any(Nk_arr < 0) or np.any(Nk_arr > N_end):
        raise ValueError(f"N outside allowed range [0, {N_end:.6f}]")
    return phi_of_N(N_end - Nk_arr)

# ----------------------------
# Build a grid for N in the range [N_min, N_max]
# ----------------------------
N_min = 20
N_max_user = 200
num_points = 1000000

# Cannot exceed the numerical inflation length N_end
N_max = min(N_max_user, N_end)

if N_min > N_max:
    raise ValueError(f"N_min={N_min} > N_end={N_end:.6f}")

N_grid = np.linspace(N_min, N_max, num_points)
phi_grid = phik(N_grid)

# ----------------------------
# Save table to disk: two columns [N, phi(N)]
# ----------------------------
tabela = np.vstack([N_grid, phi_grid]).T
np.savetxt("Python_Starobinsky.dat", tabela, fmt=["%.6f","%.6f"])

print("File Python_Starobinsky.dat saved.")


