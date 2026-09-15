"""
MTP_Experiments/scripts/phase1a_baseline.py

Phase 1A: Single Sphere Baseline Validation
Evaluates N=42 and N=162 under unit force F=1.0, eta=1.0 in unbounded Stokes fluid.
Verifies Omega approx 0, U || F, and saves comprehensive metadata to MTP_Experiments/.
"""

import sys
import os
import csv
import json
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts.common import (
    load_sphere_structure, solve_sphere_unbounded, stokes_velocity,
    angle_between_vectors, PHASE1_RAW_DATA_DIR, PHASE1_BASELINE_DIR,
    PHASE1_LOGS_DIR
)

def run_phase1a(N_list=[42, 162], sphere_type='Rh_calibrated', force_direction=(0.0, 0.0, 1.0)):
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 70)
    log(f"PHASE 1A — SINGLE SPHERE BASELINE VALIDATION ({sphere_type})")
    log("=" * 70)
    
    eta = 1.0
    F_mag = 1.0
    force_vec = F_mag * np.array(force_direction, dtype=np.float64)
    force_vec = force_vec / np.linalg.norm(force_vec) * F_mag
    orientation_quaternion = [1.0, 0.0, 0.0, 0.0]
    
    results = []
    
    for N in N_list:
        log(f"\n>>> Running Baseline for N = {N} Blobs <<<")
        r_conf, meta = load_sphere_structure(sphere_type, N)
        Rg = meta['Rg']
        Rh = meta['Rh']
        a_blob = meta['blob_radius']
        
        U, Omega, runtime = solve_sphere_unbounded(r_conf, a_blob, eta, force_vec, method='cholesky')
        
        U_mag = np.linalg.norm(U)
        Omega_mag = np.linalg.norm(Omega)
        angle_deg = angle_between_vectors(U, force_vec)
        U_theory = stokes_velocity(F_mag, eta, Rh)
        abs_error = abs(U_mag - U_theory)
        rel_error_pct = abs_error / U_theory * 100.0
        
        # Tolerances
        TOL_OMEGA = 1e-12
        TOL_ANGLE_DEG = 1e-10
        TOL_REL_ERR = 0.5
        
        pass_omega = Omega_mag < TOL_OMEGA
        pass_align = angle_deg < TOL_ANGLE_DEG
        pass_speed = rel_error_pct < TOL_REL_ERR
        overall_pass = pass_omega and pass_align and pass_speed
        
        log(f"Shape                 : sphere")
        log(f"N_blobs               : {N}")
        log(f"Blob Radius           : {a_blob:.6f}")
        log(f"R_g (Geometric)       : {Rg:.4f}")
        log(f"R_h (Hydrodynamic)    : {Rh:.4f}")
        log(f"Applied Force         : {force_vec.tolist()}")
        log(f"Fluid Viscosity (eta) : {eta:.4f}")
        log(f"Orientation (q)       : {orientation_quaternion}")
        log(f"Velocity Vector U     : {U.tolist()}")
        log(f"Angular Velocity Omega: {Omega.tolist()}")
        log(f"Simulated Speed |U|   : {U_mag:.8e}")
        log(f"Theoretical Speed U_th: {U_theory:.8e}")
        log(f"Absolute Error        : {abs_error:.4e}")
        log(f"Relative Error [%]    : {rel_error_pct:.6f}%")
        log(f"Rotation Mag |Omega|  : {Omega_mag:.4e}")
        log(f"Colinearity Angle     : {angle_deg:.4e} deg")
        log(f"Runtime               : {runtime:.4f} s")
        log(f"Check |Omega| < 1e-12 : {'PASS' if pass_omega else 'FAIL'}")
        log(f"Check U || F          : {'PASS' if pass_align else 'FAIL'}")
        log(f"Check Speed Error     : {'PASS' if pass_speed else 'FAIL'}")
        log(f"Single Test Verdict   : {'PASSED' if overall_pass else 'FAILED'}")
        
        row = {
            'shape': 'sphere',
            'N_blobs': int(N),
            'sphere_type': str(sphere_type),
            'blob_radius': float(a_blob),
            'R_g': float(Rg),
            'R_h': float(Rh),
            'force_x': float(force_vec[0]),
            'force_y': float(force_vec[1]),
            'force_z': float(force_vec[2]),
            'force_magnitude': float(F_mag),
            'viscosity': float(eta),
            'orientation': str(orientation_quaternion),
            'velocity_x': float(U[0]),
            'velocity_y': float(U[1]),
            'velocity_z': float(U[2]),
            'simulation_velocity': float(U_mag),
            'omega_x': float(Omega[0]),
            'omega_y': float(Omega[1]),
            'omega_z': float(Omega[2]),
            'angular_velocity': float(Omega_mag),
            'angle_U_F_deg': float(angle_deg),
            'theoretical_velocity': float(U_theory),
            'absolute_error': float(abs_error),
            'relative_error_percent': float(rel_error_pct),
            'runtime_seconds': float(runtime),
            'pass_omega': bool(pass_omega),
            'pass_align': bool(pass_align),
            'pass_speed': bool(pass_speed),
            'verdict': 'PASSED' if overall_pass else 'FAILED'
        }
        results.append(row)
        
    # Save CSV
    csv_file = os.path.join(PHASE1_RAW_DATA_DIR, f'phase1a_baseline_{sphere_type}.csv')
    fieldnames = list(results[0].keys())
    with open(csv_file, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    log(f"\n[Raw Data Saved] -> {csv_file}")
    
    # Save JSON summary
    json_file = os.path.join(PHASE1_BASELINE_DIR, f'phase1a_baseline_summary_{sphere_type}.json')
    with open(json_file, 'w') as fp:
        json.dump(results, fp, indent=2)
    log(f"[Summary Saved]  -> {json_file}")
    
    # Save Log
    log_file = os.path.join(PHASE1_LOGS_DIR, f'phase1a_baseline_{sphere_type}.log')
    with open(log_file, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    log(f"[Log Saved]      -> {log_file}")
    log("=" * 70)
    
    return results

if __name__ == '__main__':
    # Run for calibrated spheres
    run_phase1a(N_list=[42, 162], sphere_type='Rh_calibrated')
    # Run for geometric spheres
    run_phase1a(N_list=[42, 162], sphere_type='Rg_geometric')
