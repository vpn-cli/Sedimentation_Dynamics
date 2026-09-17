"""
MTP_Experiments/scripts/phase2c_single_cylinder.py

Phase 2C: Single Rigid Circular Cylinder Mobility, Orientation & Conditioning Audit.
Evaluates Axial (parallel to axis) and Transverse (perpendicular to axis) motion
across the established MTP refinement sequence N in {12, 42, 162, 642, 2562}.
Performs rigorous numerical conditioning audit on the RPY and Schur complement systems.
"""

import sys
import os
import time
import json
import numpy as np
import scipy.linalg

# Ensure repository and MTP scripts are on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts.common import REPO_PATH, EXPERIMENTS_DIR
from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from read_input.read_vertex_file import read_vertex_file

from scripts.cylinder_theory import (
    cylinder_theoretical_summary,
    equivalent_surface_sphere_benchmark,
    tirado_de_la_torre_benchmark,
    hubbard_douglas_capacitance_benchmark
)

# Directory structure
PHASE2_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_shapes')
STRUCTURES_DIR = os.path.join(PHASE2_DIR, 'structures')
RAW_DATA_DIR = os.path.join(PHASE2_DIR, 'raw_data')
PROCESSED_DATA_DIR = os.path.join(PHASE2_DIR, 'processed_data')
PLOTS_DIR = os.path.join(PHASE2_DIR, 'plots')
REPORTS_DIR = os.path.join(PHASE2_DIR, 'reports')
LOGS_DIR = os.path.join(PHASE2_DIR, 'logs')

for d in [RAW_DATA_DIR, PROCESSED_DATA_DIR, PLOTS_DIR, REPORTS_DIR, LOGS_DIR]:
    os.makedirs(d, exist_ok=True)

RESOLUTIONS = [12, 42, 162, 642, 2562]
R_GEOM = 1.0
L_GEOM = 2.0
ETA = 1.0
RHO = 1.0

def load_cylinder_geometry(N):
    """Loads cylinder coordinates and blob radius for given N."""
    filepath = os.path.join(STRUCTURES_DIR, f"cylinder_N_{N}_R_{R_GEOM:.1f}_L_{L_GEOM:.1f}.vertex")
    if not os.path.exists(filepath):
        # Generate on the fly if missing
        from scripts.generate_cylinder import generate_mtp_cylinder_hierarchy
        generate_mtp_cylinder_hierarchy(STRUCTURES_DIR, R_GEOM, L_GEOM)
        
    r_conf = read_vertex_file(filepath)
    with open(filepath, 'r') as fp:
        for line in fp:
            if not line.startswith('#'):
                a_blob = float(line.split()[1]) / 2.0
                break
    return r_conf, a_blob, filepath

def solve_cylinder_mobility_unbounded(r_conf, a_blob, orientation_q, force_vec, torque_vec=None, eta=ETA):
    """
    Computes instantaneous rigid-body velocity and angular velocity for a multiblob cylinder
    in unbounded Stokes flow, along with full conditioning metrics.
    """
    t0 = time.time()
    N_blobs = len(r_conf)
    b = Body(np.zeros(3), orientation_q, r_conf, a_blob)
    r_lab = b.get_r_vectors()
    K = b.calc_K_matrix()
    
    # Evaluate 3N x 3N RPY mobility matrix
    M = mb.rotne_prager_tensor(r_lab, eta, a_blob)
    
    # Cholesky factorization
    L_chol, lower = scipy.linalg.cho_factor(M)
    cholesky_success = True
    diag_L = np.diag(L_chol)
    min_diag_L = float(np.min(diag_L))
    max_diag_L = float(np.max(diag_L))
    
    # Schur complement: 6x6 resistance matrix Kt * M^-1 * K
    Minv_K = scipy.linalg.cho_solve((L_chol, lower), K, check_finite=False)
    Kt_Minv_K = np.dot(K.T, Minv_K)
    
    # Conditioning of 6x6 Schur complement
    eigvals = np.linalg.eigvalsh(Kt_Minv_K)
    min_eig = float(np.min(eigvals))
    max_eig = float(np.max(eigvals))
    cond_num = float(max_eig / min_eig)
    
    # 6x6 body mobility matrix
    N_body = np.linalg.pinv(Kt_Minv_K)
    
    wrench = np.zeros(6, dtype=np.float64)
    wrench[0:3] = force_vec
    if torque_vec is not None:
        wrench[3:6] = torque_vec
        
    U_Omega = np.dot(N_body, wrench)
    U = U_Omega[0:3]
    Omega = U_Omega[3:6]
    runtime = time.time() - t0
    
    diag_info = {
        'N_blobs': int(N_blobs),
        'matrix_dim': f"{3*N_blobs}x{3*N_blobs}",
        'cholesky_success': bool(cholesky_success),
        'min_diag_L': min_diag_L,
        'max_diag_L': max_diag_L,
        'schur_cond_num': cond_num,
        'min_eigenvalue': min_eig,
        'max_eigenvalue': max_eig,
        'runtime_seconds': float(runtime)
    }
    
    return U, Omega, diag_info, N_body

def run_single_cylinder_audit():
    """Runs single cylinder orientation and conditioning audit across all resolutions."""
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 95)
    log("PHASE 2C: SINGLE RIGID CYLINDER MOBILITY, ORIENTATION & CONDITIONING AUDIT")
    log("=" * 95)
    
    theory = cylinder_theoretical_summary(R_GEOM, L_GEOM, ETA)
    U_th_surf = theory['equivalent_surface_sphere']['velocity_unit_force']
    Rh_th_surf = theory['equivalent_surface_sphere']['resistance']
    U_th_tirado = theory['tirado_de_la_torre']['velocity_unit_force']
    U_th_hd = theory['hubbard_douglas']['velocity_unit_force']
    
    log(f"Cylinder Dimensions: Radius R={R_GEOM:.2f}, Length L={L_GEOM:.2f}, Aspect Ratio L/(2R)={theory['geometry']['aspect_ratio_p']:.2f}")
    log(f"Fluid Properties: Viscosity eta={ETA:.2f}, Density rho={RHO:.2f}")
    log(f"Theoretical Benchmarks (Unit Force F=1.0):")
    log(f"  Equivalent Surface Sphere (Primary) : U = {U_th_surf:.6f} | Rh = {Rh_th_surf:.4f}")
    log(f"  Hubbard-Douglas Capacitance         : U = {U_th_hd:.6f} | Rh = {theory['hubbard_douglas']['resistance']:.4f}")
    log(f"  Tirado & Garcia de la Torre (Mean)  : U = {U_th_tirado:.6f} | Rh = {theory['tirado_de_la_torre']['resistance']:.4f}\n")
    
    audit_results = []
    
    log(f"{'N':>6} | {'a_blob':>8} | {'U_axial':>10} | {'U_trans':>10} | {'Ratio U_tr/U_ax':>14} | {'Err Surf%':>10} | {'|Omega|':>10} | {'Cond(K^tM^-1K)':>14} | {'Time(s)':>8}")
    print("-" * 105)
    
    q_identity = Quaternion([1.0, 0.0, 0.0, 0.0])
    F_mag = 1.0
    F_axial = np.array([0.0, 0.0, F_mag])  # Axial: along z-axis
    F_trans = np.array([F_mag, 0.0, 0.0])  # Transverse: along x-axis
    
    for N in RESOLUTIONS:
        r_conf, a_blob, vfile = load_cylinder_geometry(N)
        
        # Case A: Axial
        U_ax, Om_ax, diag_ax, N_body_ax = solve_cylinder_mobility_unbounded(
            r_conf, a_blob, q_identity, F_axial, eta=ETA
        )
        u_ax = float(np.linalg.norm(U_ax))
        om_ax = float(np.linalg.norm(Om_ax))
        err_ax_surf = abs(u_ax - U_th_surf) / U_th_surf * 100.0
        
        # Case B: Transverse
        U_tr, Om_tr, diag_tr, N_body_tr = solve_cylinder_mobility_unbounded(
            r_conf, a_blob, q_identity, F_trans, eta=ETA
        )
        u_tr = float(np.linalg.norm(U_tr))
        om_tr = float(np.linalg.norm(Om_tr))
        err_tr_surf = abs(u_tr - U_th_surf) / U_th_surf * 100.0
        
        ratio = u_tr / u_ax
        mean_u = (u_ax + 2.0 * u_tr) / 3.0
        mean_err_surf = abs(mean_u - U_th_surf) / U_th_surf * 100.0
        
        log(f"{N:6d} | {a_blob:8.4f} | {u_ax:10.5f} | {u_tr:10.5f} | {ratio:14.4f} | {mean_err_surf:9.2f}% | {max(om_ax, om_tr):10.2e} | {diag_ax['schur_cond_num']:14.2f} | {diag_ax['runtime_seconds']:8.2f}")
        
        audit_results.append({
            'N': int(N),
            'a_blob': float(a_blob),
            'U_axial': float(u_ax),
            'U_transverse': float(u_tr),
            'U_isotropic_mean': float(mean_u),
            'anisotropy_ratio': float(ratio),
            'rel_error_axial_surf_pct': float(err_ax_surf),
            'rel_error_trans_surf_pct': float(err_tr_surf),
            'rel_error_mean_surf_pct': float(mean_err_surf),
            'rel_error_mean_tirado_pct': float(abs(mean_u - U_th_tirado) / U_th_tirado * 100.0),
            'rel_error_mean_hd_pct': float(abs(mean_u - U_th_hd) / U_th_hd * 100.0),
            'omega_axial_norm': float(om_ax),
            'omega_trans_norm': float(om_tr),
            'solver_diag': diag_ax
        })
        
    log("=" * 95)
    
    # Save log
    log_path = os.path.join(LOGS_DIR, 'phase2c_single_cylinder_audit.log')
    with open(log_path, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    print(f"Saved audit log: {log_path}")
    
    return audit_results

if __name__ == '__main__':
    run_single_cylinder_audit()
