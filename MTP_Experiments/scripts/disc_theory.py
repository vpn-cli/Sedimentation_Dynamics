"""
MTP_Experiments/scripts/disc_theory.py

Analytical solutions for circular discs and oblate spheroids in low-Reynolds-number Stokes flow.
References:
- Happel, J. and Brenner, H., "Low Reynolds number hydrodynamics", Martinus Nijhoff (1983).
- Kim, S. and Karrila, S. J., "Microhydrodynamics: Principles and Selected Applications", Butterworth-Heinemann (1991).
- Perrin, F., J. Phys. Radium 5, 497-511 (1934).
"""

import numpy as np

def thin_disc_resistance_perp(R, eta=1.0):
    """
    Translational resistance for an infinitely thin circular disc moving perpendicular
    to its face (broadside-on / face-on) in unbounded Stokes flow.
    F_perp = 16 * eta * R * U_perp
    """
    return 16.0 * eta * R

def thin_disc_resistance_parallel(R, eta=1.0):
    """
    Translational resistance for an infinitely thin circular disc moving parallel
    to its face (edge-on) in unbounded Stokes flow.
    F_parallel = (32 / 3) * eta * R * U_parallel
    """
    return (32.0 / 3.0) * eta * R

def thin_disc_velocity_perp(F_mag, R, eta=1.0):
    """Face-on settling velocity: U = F / (16 * eta * R)"""
    return F_mag / thin_disc_resistance_perp(R, eta)

def thin_disc_velocity_parallel(F_mag, R, eta=1.0):
    """Edge-on settling velocity: U = 3 * F / (32 * eta * R)"""
    return F_mag / thin_disc_resistance_parallel(R, eta)

def thin_disc_anisotropy_ratio():
    """
    Theoretical ratio of edge-on to face-on settling velocities under equal force:
    U_parallel / U_perp = R_perp / R_parallel = 16 / (32/3) = 1.5
    """
    return 1.5

def thin_disc_rotational_resistance_normal(R, eta=1.0):
    """Rotational resistance around symmetry axis (spinning in-plane): T = (32/3) * eta * R^3 * Omega"""
    return (32.0 / 3.0) * eta * (R**3)

def thin_disc_rotational_resistance_diametral(R, eta=1.0):
    """Rotational resistance around diametral in-plane axis (tumbling): T = (32/3) * eta * R^3 * Omega"""
    return (32.0 / 3.0) * eta * (R**3)

def oblate_spheroid_resistance(a, c, eta=1.0):
    """
    Translational resistance for an oblate spheroid with semi-major axis a and semi-minor axis c (c < a).
    Returns (R_perp, R_parallel).
    """
    if c >= a:
        raise ValueError("For oblate spheroid, semi-minor axis c must be strictly less than a.")
    e = np.sqrt(1.0 - (c / a)**2)
    asin_e = np.arcsin(e)
    sqrt_1_e2 = np.sqrt(1.0 - e**2)
    
    # Perpendicular (along axis of symmetry c)
    denom_perp = (2.0 * e**2 - 1.0) * asin_e + e * sqrt_1_e2
    R_perp = 16.0 * np.pi * eta * a * (e**3) / denom_perp
    
    # Parallel (perpendicular to axis of symmetry, along a)
    denom_parallel = (2.0 * e**2 + 1.0) * asin_e - e * sqrt_1_e2
    R_parallel = 32.0 * np.pi * eta * a * (e**3) / denom_parallel
    
    return R_perp, R_parallel

def tilted_disc_unbounded_velocity(force_vec, normal_vec, R, eta=1.0):
    """
    Calculates the exact theoretical velocity vector for a disc at arbitrary orientation
    in an unbounded Stokes fluid under applied force vector force_vec.
    
    Decomposes force into normal and in-plane tangential components:
    F = F_perp * n + F_parallel
    U = (F_perp / R_perp) * n + (F_parallel / R_parallel)
    """
    F = np.asarray(force_vec, dtype=np.float64)
    n = np.asarray(normal_vec, dtype=np.float64)
    n = n / np.linalg.norm(n)
    
    F_perp_mag = np.dot(F, n)
    F_perp_vec = F_perp_mag * n
    F_parallel_vec = F - F_perp_vec
    
    R_perp = thin_disc_resistance_perp(R, eta)
    R_parallel = thin_disc_resistance_parallel(R, eta)
    
    U_perp_vec = F_perp_vec / R_perp
    U_parallel_vec = F_parallel_vec / R_parallel
    U_total = U_perp_vec + U_parallel_vec
    
    return U_total, U_perp_vec, U_parallel_vec
