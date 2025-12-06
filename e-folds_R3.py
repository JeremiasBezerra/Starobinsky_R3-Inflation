# -*- coding: utf-8 -*-
"""
@author: J. Bezerra-Sobrinho
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d


# -------------------------------------------------------
# Constants
# -------------------------------------------------------
PI = np.pi
a = np.sqrt(16 * PI / 3)   # Model-dependent constant
V0 = 1.0                   # Overall normalization


# -------------------------------------------------------
# R^3 potential V(phi, α) and its derivative V'(phi)
#
# The potential contains a square root term:
#   sqrt( 1 - 4α (1 - e^{aφ}) )
#
# For numerical stability, we check that the inside of the
# square root is non-negative; otherwise return NaN.
# -------------------------------------------------------
def V(phi, alpha):
    exp_term = np.exp(a * phi)
    inside_sqrt = 1 - 4 * alpha * (1 - exp_term)
    
    if np.any(inside_sqrt < 0):
        return np.nan  # Prevent invalid sqrt
    
    sqrt_term = np.sqrt(inside_sqrt)
    
    return (
        (1/12) * V0 / alpha**2 * np.exp(-2*a*phi)
        * (1 - sqrt_term)
        * (-1 + 8*alpha*(1 - exp_term) + sqrt_term)
    )


def Vprime(phi, alpha):
    exp_term = np.exp(a * phi)
    inside_sqrt = 1 - 4 * alpha * (1 - exp_term)
    
    if np.any(inside_sqrt < 0):
        return np.nan
    
    sqrt_term = np.sqrt(inside_sqrt)
    
    return (
        (1/6) * V0 / alpha**2 * a * np.exp(-2*a*phi)
        * (1 - sqrt_term)
        * (1 - 2*alpha*(4 - exp_term) - sqrt_term)
    )


# -------------------------------------------------------
# Solve the inflationary dynamics for the R^3 model.
#
# The system of equations is:
#   dφ/dt = π
#   dπ/dt = -3Hπ - V'(φ)
#   dN/dt = H
#
# The Hubble rate is given by:
#   H^2 = (8π/3) [ 1/2 π^2 + V(φ) ]
#
# Inflation ends when ε_H = 1, where:
#   ε_H = 4π π^2 / H^2
# -------------------------------------------------------
def solve_R3(alpha, phi0=1.54, pi0=0.0, N0=0.0):

    # Avoid α = 0, which makes the potential ill-defined
    if alpha == 0.0:
        alpha = 1e-10

    # System of ODEs
    def derivs(t, y):
        phi, pi, Nefold = y
        H = np.sqrt((8 * PI / 3) * (0.5*pi**2 + V(phi, alpha)))
        dphi = pi
        dpi  = -3*H*pi - Vprime(phi, alpha)
        dN   = H
        return [dphi, dpi, dN]

    # End of inflation: ε_H = 1
    def stop_event(t, y):
        phi, pi, Nefold = y
        H = np.sqrt((8 * PI / 3) * (0.5*pi**2 + V(phi, alpha)))
        epsilon_H = 4 * PI * pi**2 / H**2
        return 1 - epsilon_H  # Crosses zero at inflation end

    stop_event.terminal = True    # Stop integration
    stop_event.direction = -1     # Detect decreasing to zero

    # Solve system
    y0 = [phi0, pi0, N0]
    sol = solve_ivp(
        derivs, [0, 1e6], y0,
        events=stop_event,
        max_step=0.8,
        rtol=1e-9, atol=1e-12
    )

    phi_t = sol.y[0]
    N_t   = sol.y[2]
    Nend  = N_t[-1]  # Total number of e-folds generated

    # ---------------------------------------------------
    # Build cubic interpolation:
    #   φ(N) with N = Nend - N(t)
    #
    # This returns the field value when N e-folds remain
    # before the end of inflation.
    # ---------------------------------------------------
    phi_of_N_interp = interp1d(
        Nend - N_t,
        phi_t,
        kind='cubic',
        bounds_error=False,
        fill_value="extrapolate"
    )

    return Nend, phi_of_N_interp


# -------------------------------------------------------
# α scan (prior range)
# -------------------------------------------------------
alpha_min = -1e-9
alpha_max = -1e-10
num_alpha = 100000

alpha_grid = np.linspace(alpha_min, alpha_max, num_alpha)

# -------------------------------------------------------
# Prior for N (remaining e-folds)
# -------------------------------------------------------
Nmin = 20
Nmax_user = 200
num_points = 100000

tabela_completa = []


# -------------------------------------------------------
# Compute φ(N, α) for each α in the grid
# -------------------------------------------------------
for alpha in alpha_grid:
    print(f"Calculating for alpha = {alpha:.2e}")
    
    Nend, phi_of_N = solve_R3(alpha)
    Nmax = min(Nmax_user, Nend)
    
    Ngrid = np.linspace(Nmin, Nmax, num_points)
    
    for N in Ngrid:
        phi_val = phi_of_N(N)
        tabela_completa.append([N, alpha, phi_val])


# -------------------------------------------------------
# Export final table: columns = [N, α, φ(N, α)]
# -------------------------------------------------------
tabela_completa = np.array(tabela_completa)
np.savetxt("Python_R3.dat", tabela_completa,
           fmt="%.6f %.6e %.6f",
           comments="")

print("Table successfully exported!")
