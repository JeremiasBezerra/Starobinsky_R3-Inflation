# Starobinsky_R3-Inflation
This repository provides the implementation and numerical tools used to evaluate the inflationary dynamics of the Starobinsky and 𝑅3 models within the CLASS Boltzmann code framework.

The main objective is to supply the inflaton potential in CLASS using a local Taylor expansion around the field value 𝜙_* corresponding to the horizon-crossing of a given mode k. Since CLASS requires the potential and its derivatives evaluated at a single field point, the code computes the coefficients V0, V1, V2, V3 and V4, directly from the exact analytical expressions of the potentials.

A key ingredient of both inflationary models is the determination of 𝜙_*, defined by the condition N_k = Integrate H * dt.

which represents the number of e-folds between horizon crossing and the end of inflation. Although in the literature this mapping is commonly obtained through slow-roll approximations, here we avoid such assumptions and compute 𝑁(𝜙) exactly by numerically solving the full background equations.

Two Python scripts included in the repository perform this calculation:

Starobinsky model:
The script solves the system (𝜙, d𝜙/dt, N) using solve_ivp, identifies the end of inflation using the Hubble slow-roll parameter ϵ=1, and then constructs an interpolating function 𝜙(N). A table containing pairs (N, ϕ) is saved, and CLASS reads this table to determine 𝜙_* for any chosen value of N_k.

R3 model:
The script follows the same procedure but evaluates the exact R3 potential and its derivative, which contain a square-root structure depending on the parameter α. A table containing pairs (N, α, ϕ) is saved, and CLASS reads this table to determine 𝜙_* for any chosen value of N_k.

Inside CLASS, the edited primordial.c reads the precomputed 𝜙_* value from the tables and substitutes it into the Taylor-expanded potential. This modular approach cleanly separates the numerical background integration (handled externally in Python) from the perturbation and primordial spectrum computation performed internally by CLASS.
