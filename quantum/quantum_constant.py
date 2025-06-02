# -*- coding: utf-8 -*-
"""
Created for quantum correction implementation

@author:Xiaoyan Xu
        (XIDIAN UNIVERSITY)
Quantum correction constants for Density Gradient model implementation.
Contains all physical constants and parameters needed for quantum potential calculation.
"""

import numpy as np
import constants as const

delta_n = 1e-15
delta_p = 1e-15

# Fundamental quantum constants
hbar = 1.054571817e-34  # Reduced Planck's constant, J·s
hbar_eV = 6.582119569e-16  # Reduced Planck's constant, eV·s

# Default effective masses (in units of free electron mass m0)
# These should be adjusted based on the specific semiconductor material
m0 = 9.1093837015e-31  # Free electron mass, kg
m_n_star = 0.26 * m0   # Electron effective mass (typical for GaAs), kg
m_p_star = 0.38 * m0   # Hole effective mass (typical for GaAs), kg

# Quantum correction scaling factors
# These parameters r_n and r_p are material-dependent scaling factors
# For most semiconductors, they are close to unity
r_n = 1.0  # Electron quantum correction scaling factor
r_p = 1.0  # Hole quantum correction scaling factor

# Calculate quantum correction coefficients b_n and b_p
# b_n = (ℏ²r_n)/(12*q*m_n*)
# b_p = (ℏ²r_p)/(12*q*m_p*)

def calculate_b_n():
    """
    Calculate electron quantum correction coefficient b_n
    
    Formula: b_n = (ℏ²r_n)/(12*q*m_n*)
    
    Returns:
    --------
    float : b_n in units of J·m²
    """
    b_n = (hbar**2 * r_n) / (12.0 * const.q * m_n_star)
    return b_n

def calculate_b_p():
    """
    Calculate hole quantum correction coefficient b_p
    
    Formula: b_p = (ℏ²r_p)/(12*q*m_p*)
    
    Returns:
    --------
    float : b_p in units of J·m²
    """
    b_p = (hbar**2 * r_p) / (12.0 * const.q * m_p_star)
    return b_p

# Pre-calculate the coefficients for efficiency
b_n = calculate_b_n()
b_p = calculate_b_p()

# Convert to more convenient units for numerical computation
# Since we work with thermal voltage Vt = kT/q, it's convenient to express
# b_n and b_p in terms of Vt and length units
b_n_normalized = b_n / (const.Vt * const.q)  # in units of m²
b_p_normalized = b_p / (const.Vt * const.q)  # in units of m²

print(f"Quantum correction coefficients calculated:")
print(f"b_n = {b_n:.6e} J·m²")
print(f"b_p = {b_p:.6e} J·m²")
print(f"b_n_normalized = {b_n_normalized:.6e} m²")
print(f"b_p_normalized = {b_p_normalized:.6e} m²")