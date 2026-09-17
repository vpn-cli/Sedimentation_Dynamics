"""
MTP_Experiments/scripts/cylinder_theory.py

Analytical and theoretical hydrodynamic benchmarks for straight circular cylinders
in low-Reynolds-number Stokes flow.

CRITICAL THEORETICAL FOUNDATION:
Unlike a sphere (Stokes 1851: F = 6*pi*eta*R*U) or an infinitely thin circular disc
(Oberbeck 1876, Happel & Brenner 1983: F_perp = 16*eta*R*U_perp, F_par = 32/3*eta*R*U_par),
there is NO EXACT CLOSED-FORM ANALYTICAL SOLUTION for the Stokes resistance of a finite
circular cylinder with flat end caps.

This module provides established theoretical benchmarks:
1. Short-Cylinder Bead-Shell Theory (Tirado & Garcia de la Torre 1979, 1984)
2. Slender-Body Theory (Broersma 1960) [and documentation of its breakdown for L/D = 1]
3. Equivalent Surface-Area Sphere (Boundary Shear / Multiblob Analogy)
4. Equivalent Volume Sphere (Volume Displacement Analogy)
5. Hubbard-Douglas Electrostatic Capacitance Analogy (Douglas, Mansfield & Garboczi 1996)

Orientation Conventions:
- Axial motion (parallel to symmetry axis z): Force F_z || e_axis, velocity U_parallel
- Transverse motion (perpendicular to symmetry axis z): Force F_x perp e_axis, velocity U_perp
"""

import numpy as np

def cylinder_geometry_properties(R=1.0, L=2.0):
    """Computes geometric properties of a circular cylinder of radius R and length L."""
    D = 2.0 * R
    aspect_ratio_p = L / D  # p = L / (2R)
    slenderness_lambda = L / R
    
    # Surface areas
    A_lat = 2.0 * np.pi * R * L
    A_caps = 2.0 * np.pi * (R**2)
    A_total = A_lat + A_caps
    
    # Volume
    Volume = np.pi * (R**2) * L
    
    # Projected areas
    Area_proj_axial = np.pi * (R**2)
    Area_proj_transverse = 2.0 * R * L
    
    # Theoretical surface radius of gyration
    # <r^2> = (A_lat * (R^2 + L^2/12) + A_caps * (R^2/2 + L^2/4)) / A_total
    mean_r2 = (A_lat * (R**2 + (L**2)/12.0) + A_caps * ((R**2)/2.0 + (L**2)/4.0)) / A_total
    Rg_surface = np.sqrt(mean_r2)
    
    return {
        'R': R,
        'L': L,
        'D': D,
        'aspect_ratio_p': aspect_ratio_p,
        'slenderness_lambda': slenderness_lambda,
        'A_lat': A_lat,
        'A_caps': A_caps,
        'A_total': A_total,
        'Volume': Volume,
        'Area_proj_axial': Area_proj_axial,
        'Area_proj_transverse': Area_proj_transverse,
        'Area_proj_ratio': Area_proj_transverse / Area_proj_axial,
        'Rg_surface': Rg_surface
    }

def equivalent_surface_sphere_benchmark(R=1.0, L=2.0, eta=1.0):
    """
    Theoretical benchmark based on a sphere having the exact same surface area as the cylinder.
    A_total = 2*pi*R*L + 2*pi*R^2 = 6*pi (for R=1, L=2).
    A_sphere = 4*pi*R_eff^2 = 6*pi => R_eff = sqrt(1.5) * R approx 1.22474 * R.
    Resistance: R_h = 6 * pi * eta * R_eff.
    Settling speed: U = F / R_h.
    """
    A_total = 2.0 * np.pi * R * L + 2.0 * np.pi * (R**2)
    R_eff = np.sqrt(A_total / (4.0 * np.pi))
    Rh = 6.0 * np.pi * eta * R_eff
    mobility = 1.0 / Rh
    return {
        'model': 'equivalent_surface_sphere',
        'R_eff': float(R_eff),
        'hydrodynamic_radius': float(R_eff),
        'resistance': float(Rh),
        'mobility': float(mobility),
        'velocity_unit_force': float(mobility)
    }

def equivalent_volume_sphere_benchmark(R=1.0, L=2.0, eta=1.0):
    """
    Theoretical benchmark based on a sphere having the exact same volume as the cylinder.
    V_cyl = pi * R^2 * L = 2*pi (for R=1, L=2).
    V_sphere = (4/3)*pi*R_eff^3 = 2*pi => R_eff = (3/2)^(1/3) * R approx 1.144714 * R.
    """
    V_cyl = np.pi * (R**2) * L
    R_eff = (3.0 * V_cyl / (4.0 * np.pi))**(1.0 / 3.0)
    Rh = 6.0 * np.pi * eta * R_eff
    mobility = 1.0 / Rh
    return {
        'model': 'equivalent_volume_sphere',
        'R_eff': float(R_eff),
        'hydrodynamic_radius': float(R_eff),
        'resistance': float(Rh),
        'mobility': float(mobility),
        'velocity_unit_force': float(mobility)
    }

def hubbard_douglas_capacitance_benchmark(R=1.0, L=2.0, eta=1.0):
    """
    Hubbard-Douglas capacitance benchmark (Mansfield, Douglas & Garboczi 1996; Hubbard & Douglas 1993).
    For a closed cylinder of aspect ratio L/D = 1, electrostatic capacitance calculation yields:
    C / R_vol = 1.0121(6).
    Hydrodynamic resistance: R_h = 6 * pi * eta * C = 6 * pi * eta * (1.0121 * R_vol).
    """
    vol_bench = equivalent_volume_sphere_benchmark(R, L, eta)
    R_vol = vol_bench['R_eff']
    C_cap = 1.0121 * R_vol
    Rh = 6.0 * np.pi * eta * C_cap
    mobility = 1.0 / Rh
    return {
        'model': 'hubbard_douglas_capacitance',
        'C_capacitance': float(C_cap),
        'hydrodynamic_radius': float(C_cap),
        'resistance': float(Rh),
        'mobility': float(mobility),
        'velocity_unit_force': float(mobility)
    }

def tirado_de_la_torre_benchmark(R=1.0, L=2.0, eta=1.0):
    """
    Tirado & Garcia de la Torre (1979, 1984) bead-shell extrapolation model for short cylinders.
    Aspect ratio: p = L / (2R) = 1.0.
    End-effect correction: nu = 0.312 + 0.565/p - 0.100/p^2 = 0.7770 for p=1.0.
    Mean translational friction:
    xi_mean = 3 * pi * eta * L / (ln(p) + nu).
    For p = 1.0, ln(p) = 0 => xi_mean = 3 * pi * eta * L / nu.
    """
    p = L / (2.0 * R)
    nu = 0.312 + 0.565 / p - 0.100 / (p**2)
    xi_mean = (3.0 * np.pi * eta * L) / (np.log(p) + nu)
    Rh = xi_mean
    mobility = 1.0 / Rh
    return {
        'model': 'tirado_garcia_de_la_torre_1979_1984',
        'p': float(p),
        'nu_correction': float(nu),
        'resistance': float(Rh),
        'mobility': float(mobility),
        'velocity_unit_force': float(mobility),
        'hydrodynamic_radius': float(Rh / (6.0 * np.pi * eta))
    }

def broersma_benchmark(R=1.0, L=2.0, eta=1.0):
    """
    Broersma (1960) slender-body formulas.
    sigma = ln(2a/b) = ln(L/R) = ln(2) approx 0.69315.
    End corrections:
    gamma_par = 1.30 - 8 * (1/sigma - 0.30)^2
    gamma_perp = 0.35 - 4 * (1/sigma - 0.43)^2
    IMPORTANT NOTE: Broersma explicitly notes validity requires sigma > 2 (i.e. L/R > 7.4).
    For L/R = 2, 1/sigma approx 1.44, outside validity range.
    """
    sigma = np.log(L / R)
    gamma_par = 1.30 - 8.0 * (1.0 / sigma - 0.30)**2
    gamma_perp = 0.35 - 4.0 * (1.0 / sigma - 0.43)**2
    
    # Denominators
    denom_par = sigma - gamma_par
    denom_perp = sigma - gamma_perp
    
    F_par = (2.0 * np.pi * eta * L) / denom_par if denom_par > 0 else np.nan
    F_perp = (4.0 * np.pi * eta * L) / denom_perp if denom_perp > 0 else np.nan
    
    return {
        'model': 'broersma_1960',
        'sigma': float(sigma),
        'valid_range': 'sigma > 2 (L/R > 7.4)',
        'is_valid_for_current_geometry': False,
        'gamma_par': float(gamma_par),
        'gamma_perp': float(gamma_perp),
        'resistance_par': float(F_par),
        'resistance_perp': float(F_perp)
    }

def cylinder_theoretical_summary(R=1.0, L=2.0, eta=1.0):
    """Consolidated summary of all theoretical approximations."""
    geom = cylinder_geometry_properties(R, L)
    surf_bench = equivalent_surface_sphere_benchmark(R, L, eta)
    vol_bench = equivalent_volume_sphere_benchmark(R, L, eta)
    hd_bench = hubbard_douglas_capacitance_benchmark(R, L, eta)
    tirado_bench = tirado_de_la_torre_benchmark(R, L, eta)
    broersma_bench = broersma_benchmark(R, L, eta)
    
    return {
        'geometry': geom,
        'equivalent_surface_sphere': surf_bench,
        'equivalent_volume_sphere': vol_bench,
        'hubbard_douglas': hd_bench,
        'tirado_de_la_torre': tirado_bench,
        'broersma': broersma_bench,
        # Recommended primary benchmark for isotropic average: Equivalent Surface Sphere (Rh = 23.0859)
        # Because multiblob surface discretization directly models boundary skin friction
        'primary_benchmark': surf_bench
    }

if __name__ == '__main__':
    summary = cylinder_theoretical_summary(1.0, 2.0, 1.0)
    print("=" * 80)
    print("THEORETICAL BENCHMARKS FOR RIGID CIRCULAR CYLINDER (R=1.0, L=2.0, eta=1.0)")
    print("=" * 80)
    print("1. Equivalent Surface-Area Sphere (Boundary Shell Analogy):")
    print(f"   R_eff = {summary['equivalent_surface_sphere']['R_eff']:.5f} | Rh = {summary['equivalent_surface_sphere']['resistance']:.4f} | U = {summary['equivalent_surface_sphere']['velocity_unit_force']:.6f}")
    print("2. Hubbard-Douglas Capacitance / BEM Analogy (Mansfield & Douglas 1996):")
    print(f"   C_cap = {summary['hubbard_douglas']['C_capacitance']:.5f} | Rh = {summary['hubbard_douglas']['resistance']:.4f} | U = {summary['hubbard_douglas']['velocity_unit_force']:.6f}")
    print("3. Tirado & Garcia de la Torre (1979/1984 Short Cylinder Bead-Shell):")
    print(f"   nu = {summary['tirado_de_la_torre']['nu_correction']:.4f} | Rh = {summary['tirado_de_la_torre']['resistance']:.4f} | U = {summary['tirado_de_la_torre']['velocity_unit_force']:.6f}")
    print("4. Equivalent Volume Sphere (Displacement Analogy):")
    print(f"   R_vol = {summary['equivalent_volume_sphere']['R_eff']:.5f} | Rh = {summary['equivalent_volume_sphere']['resistance']:.4f} | U = {summary['equivalent_volume_sphere']['velocity_unit_force']:.6f}")
    print("5. Broersma (1960 Slender-Body Formulation):")
    print(f"   sigma = {summary['broersma']['sigma']:.4f} | Status: OUT OF VALIDITY RANGE ({summary['broersma']['valid_range']})")
    print("=" * 80)
