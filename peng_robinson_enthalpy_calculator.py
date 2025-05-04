
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Constants
R = 8.314  # J/mol.K

def Cp_ideal_gas(T, A, B, C, D):
    return A + B*T + C*T**2 + D*T**3  # J/mol·K

def calculate_enthalpy(T, P, Tc, Pc, omega, Cp_coeffs, Tref=298.15):
    # Step 1: Calculate Peng-Robinson constants a, b
    a = 0.45724 * R**2 * Tc**2 / Pc
    b = 0.07780 * R * Tc / Pc

    # Step 2: Calculate alpha and kappa
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    alpha = (1 + kappa * (1 - np.sqrt(T / Tc)))**2

    # Step 3: Calculate A and B
    A = a * alpha * P / (R**2 * T**2)
    B = b * P / (R * T)

    # Step 4: Solve cubic EOS for Z
    coeffs = [1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)]
    Z_roots = np.roots(coeffs)
    Z_real = Z_roots[np.isreal(Z_roots)].real
    Z = max(Z_real)  # Assume vapor phase

    # Step 5: Departure enthalpy
    term1 = R * T * (Z - 1)
    term2 = (a * alpha) / (2 * np.sqrt(2) * b)
    term3 = np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
    H_dep = term1 + term2 * term3  # J/mol

    # Step 6: Ideal gas enthalpy
    T_range = np.linspace(Tref, T, 100)
    Cp_vals = Cp_ideal_gas(T_range, *Cp_coeffs)
    H_ideal = np.trapz(Cp_vals, T_range)

    # Step 7: Total enthalpy
    H_total = H_ideal + H_dep

    return {
        "Z": Z,
        "H_dep (J/mol)": H_dep,
        "H_ideal (J/mol)": H_ideal,
        "H_total (J/mol)": H_total
    }

# Example for methane
Tc = 190.6     # K
Pc = 4.599e6   # Pa
omega = 0.011
T = 300        # K
P = 5e6        # Pa
Cp_coeffs = [19.89, 5.024e-2, 1.269e-5, -11.01e-9]  # A, B, C, D

results = calculate_enthalpy(T, P, Tc, Pc, omega, Cp_coeffs)
for k, v in results.items():
    print(f"{k}: {v:.4f}")
