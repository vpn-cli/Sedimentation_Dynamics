"""
MTP_Experiments/scripts/phase2a_single_disc.py

Phase 2A: Single Rigid Disc Mobility & Resolution Validation
Evaluates Face-on, Edge-on, and Tilted orientations across resolutions N in {20, 39, 95, 177, 347, 755}.
Compares numerical terminal velocities and anisotropy ratio against analytical Stokes theory.
Saves comprehensive metadata, CSVs, summary JSON, and publication plots.
"""

import sys
import os
import glob
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt

# Setup path and Python 3.12+ imp shim via common.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts.common import REPO_PATH, EXPERIMENTS_DIR
from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from read_input.read_vertex_file import read_vertex_file
import scipy.linalg

from scripts.disc_theory import (
    thin_disc_resistance_perp, thin_disc_resistance_parallel,
    thin_disc_velocity_perp, thin_disc_velocity_parallel,
    thin_disc_anisotropy_ratio, tilted_disc_unbounded_velocity
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

def solve_disc_mobility_unbounded(r_conf, a_blob, orientation_q, force_vec, torque_vec=None, eta=1.0):
    """
    Computes instantaneous rigid-body velocity and angular velocity for a multiblob disc
    with given orientation quaternion and applied wrench in unbounded Stokes flow.
    """
    t0 = time.time()
    b = Body(np.zeros(3), orientation_q, r_conf, a_blob)
    r_lab = b.get_r_vectors()
    K = b.calc_K_matrix()
    
    # Evaluate Rotne-Prager-Yamakawa mobility matrix
    M = mb.rotne_prager_tensor(r_lab, eta, a_blob)
    L, lower = scipy.linalg.cho_factor(M)
    
    # 6x6 body mobility matrix
    Kt_Minv_K = np.dot(K.T, scipy.linalg.cho_solve((L, lower), K, check_finite=False))
    N_body = np.linalg.pinv(Kt_Minv_K)
    
    wrench = np.zeros(6, dtype=np.float64)
    wrench[0:3] = force_vec
    if torque_vec is not None:
        wrench[3:6] = torque_vec
        
    U_Omega = np.dot(N_body, wrench)
    U = U_Omega[0:3]
    Omega = U_Omega[3:6]
    runtime = time.time() - t0
    
    return U, Omega, runtime, N_body

def run_phase2a_resolution_study():
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 80)
    log("PHASE 2A: SINGLE RIGID DISC MOBILITY & RESOLUTION CONVERGENCE")
    log("=" * 80)
    
    R_geom = 1.0
    eta = 1.0
    F_mag = 1.0
    gravity_force = np.array([0.0, 0.0, -F_mag])
    
    # Theory
    U_th_perp = thin_disc_velocity_perp(F_mag, R_geom, eta)        # 1/16 = 0.0625
    U_th_parallel = thin_disc_velocity_parallel(F_mag, R_geom, eta) # 3/32 = 0.09375
    anisotropy_th = thin_disc_anisotropy_ratio()                    # 1.5000
    
    log(f"Physical Parameters:")
    log(f"  Disc Radius R       : {R_geom}")
    log(f"  Fluid Viscosity eta : {eta}")
    log(f"  Force Magnitude F   : {F_mag}")
    log(f"  Theoretical U_perp  : {U_th_perp:.6f}")
    log(f"  Theoretical U_par   : {U_th_parallel:.6f}")
    log(f"  Theoretical Ratio   : {anisotropy_th:.4f}\n")
    
    # Find all generated disc vertex files
    vertex_files = sorted(
        glob.glob(os.path.join(STRUCTURES_DIR, "disc_N_*.vertex")),
        key=lambda p: int(os.path.basename(p).split('_N_')[1].split('_')[0])
    )
    
    if not vertex_files:
        raise FileNotFoundError("No disc vertex files found in structures directory.")
        
    results_records = []
    summary_by_resolution = []
    
    for vfile in vertex_files:
        filename = os.path.basename(vfile)
        r_conf = read_vertex_file(vfile)
        N = len(r_conf)
        
        # Read a_blob from header
        with open(vfile, 'r') as fp:
            for line in fp:
                if not line.startswith('#'):
                    a_blob = float(line.split()[1]) / 2.0
                    break
                    
        log("-" * 80)
        log(f">>> Discretization: N = {N:4d} Blobs | a_blob = {a_blob:.5f} | File: {filename} <<<")
        
        # -------------------------------------------------------------
        # Case A: Face-on (normal along z: q = [1, 0, 0, 0])
        # -------------------------------------------------------------
        q_face = Quaternion([1.0, 0.0, 0.0, 0.0])
        U_face, Om_face, rt_face, _ = solve_disc_mobility_unbounded(
            r_conf, a_blob, q_face, gravity_force, eta=eta
        )
        Uz_face = abs(U_face[2])
        err_face = abs(Uz_face - U_th_perp) / U_th_perp * 100.0
        
        # -------------------------------------------------------------
        # Case B: Edge-on (rotated 90 deg about y: q = [cos(pi/4), 0, sin(pi/4), 0])
        # -------------------------------------------------------------
        q_edge = Quaternion([np.cos(np.pi/4), 0.0, np.sin(np.pi/4), 0.0])
        U_edge, Om_edge, rt_edge, _ = solve_disc_mobility_unbounded(
            r_conf, a_blob, q_edge, gravity_force, eta=eta
        )
        Uz_edge = abs(U_edge[2])
        err_edge = abs(Uz_edge - U_th_parallel) / U_th_parallel * 100.0
        
        # Anisotropy ratio
        ratio_sim = Uz_edge / Uz_face
        err_ratio = abs(ratio_sim - anisotropy_th) / anisotropy_th * 100.0
        
        # Reynolds numbers (Re = rho * U * 2R / eta with rho=1)
        Re_face = 1.0 * Uz_face * (2.0 * R_geom) / eta
        Re_edge = 1.0 * Uz_edge * (2.0 * R_geom) / eta
        
        log(f"  [Face-on]  Uz = {Uz_face:.6f} | Theory = {U_th_perp:.6f} | Rel Error = {err_face:6.2f}% | |Omega| = {np.linalg.norm(Om_face):.2e} | Re = {Re_face:.4e}")
        log(f"  [Edge-on]  Uz = {Uz_edge:.6f} | Theory = {U_th_parallel:.6f} | Rel Error = {err_edge:6.2f}% | |Omega| = {np.linalg.norm(Om_edge):.2e} | Re = {Re_edge:.4e}")
        log(f"  [Anisotropy Ratio] U_edge / U_face = {ratio_sim:.4f} | Theory = {anisotropy_th:.4f} | Rel Error = {err_ratio:6.2f}%")
        
        # -------------------------------------------------------------
        # Case C: Tilted (theta = 45 degrees about y)
        # -------------------------------------------------------------
        theta_tilt_deg = 45.0
        theta_rad = np.radians(theta_tilt_deg)
        q_tilt = Quaternion([np.cos(theta_rad / 2.0), 0.0, np.sin(theta_rad / 2.0), 0.0])
        U_tilt, Om_tilt, rt_tilt, _ = solve_disc_mobility_unbounded(
            r_conf, a_blob, q_tilt, gravity_force, eta=eta
        )
        drift_angle_deg = np.degrees(np.arctan2(abs(U_tilt[0]), abs(U_tilt[2])))
        
        # Theoretical drift angle for this ratio:
        # tan(phi) = (ratio - 1)*sin(theta)*cos(theta) / (sin^2(theta) + ratio*cos^2(theta))
        tan_phi_th_ideal = (anisotropy_th - 1.0) * np.sin(theta_rad) * np.cos(theta_rad) / (
            np.sin(theta_rad)**2 + anisotropy_th * np.cos(theta_rad)**2
        )
        phi_th_ideal_deg = np.degrees(np.arctan(tan_phi_th_ideal))
        
        log(f"  [Tilted 45°] Ux = {U_tilt[0]:+.6f} | Uz = {U_tilt[2]:.6f} | Drift Angle = {drift_angle_deg:.2f}° (Ideal Th: {phi_th_ideal_deg:.2f}°)")
        
        res_summary = {
            'N': int(N),
            'a_blob': float(a_blob),
            'a_blob_over_R': float(a_blob / R_geom),
            'Uz_face': float(Uz_face),
            'Uz_edge': float(Uz_edge),
            'ratio_edge_to_face': float(ratio_sim),
            'rel_error_face_pct': float(err_face),
            'rel_error_edge_pct': float(err_edge),
            'rel_error_ratio_pct': float(err_ratio),
            'drift_angle_45deg': float(drift_angle_deg),
            'Re_face': float(Re_face),
            'Re_edge': float(Re_edge),
            'runtime_face_s': float(rt_face),
            'runtime_edge_s': float(rt_edge)
        }
        summary_by_resolution.append(res_summary)
        
        # Save detailed case rows
        for case_name, U_vec, Om_vec, q_val in [
            ('face_on', U_face, Om_face, q_face),
            ('edge_on', U_edge, Om_edge, q_edge),
            ('tilted_45', U_tilt, Om_tilt, q_tilt)
        ]:
            results_records.append({
                'shape': 'disc',
                'N_blobs': int(N),
                'blob_radius': float(a_blob),
                'case': case_name,
                'force_mag': float(F_mag),
                'viscosity': float(eta),
                'quaternion': str(list(q_val.entries)),
                'Ux': float(U_vec[0]),
                'Uy': float(U_vec[1]),
                'Uz': float(U_vec[2]),
                'U_mag': float(np.linalg.norm(U_vec)),
                'Omega_x': float(Om_vec[0]),
                'Omega_y': float(Om_vec[1]),
                'Omega_z': float(Om_vec[2]),
                'Omega_mag': float(np.linalg.norm(Om_vec)),
                'Reynolds_number': float(1.0 * np.linalg.norm(U_vec) * 2.0 * R_geom / eta)
            })

    # Save CSV and JSON
    csv_path = os.path.join(RAW_DATA_DIR, 'phase2a_disc_cases.csv')
    with open(csv_path, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(results_records[0].keys()))
        writer.writeheader()
        writer.writerows(results_records)
    log(f"\n[Raw Data Saved]   -> {csv_path}")
    
    json_path = os.path.join(PROCESSED_DATA_DIR, 'phase2a_disc_resolution_summary.json')
    with open(json_path, 'w') as fp:
        json.dump(summary_by_resolution, fp, indent=2)
    log(f"[Summary Saved]    -> {json_path}")
    
    # -------------------------------------------------------------
    # Tilt Angle Sweep for finest resolution (N = 755)
    # -------------------------------------------------------------
    log("\n" + "=" * 80)
    log("TILT ANGLE SWEEP (0° to 90°) FOR RESOLUTION N=755")
    log("=" * 80)
    vfile_fine = vertex_files[-1]
    r_conf_fine = read_vertex_file(vfile_fine)
    with open(vfile_fine, 'r') as fp:
        for line in fp:
            if not line.startswith('#'):
                a_blob_fine = float(line.split()[1]) / 2.0
                break
                
    tilt_angles_deg = np.linspace(0.0, 90.0, 19)
    tilt_sweep_data = []
    
    for ang_deg in tilt_angles_deg:
        ang_rad = np.radians(ang_deg)
        q_t = Quaternion([np.cos(ang_rad / 2.0), 0.0, np.sin(ang_rad / 2.0), 0.0])
        U_t, Om_t, _, _ = solve_disc_mobility_unbounded(r_conf_fine, a_blob_fine, q_t, gravity_force, eta=eta)
        drift_ang = np.degrees(np.arctan2(abs(U_t[0]), abs(U_t[2])))
        
        # Analytical prediction
        n_vec = np.array([np.sin(ang_rad), 0.0, np.cos(ang_rad)])
        U_th_vec, _, _ = tilted_disc_unbounded_velocity(gravity_force, n_vec, R_geom, eta=eta)
        drift_ang_th = np.degrees(np.arctan2(abs(U_th_vec[0]), abs(U_th_vec[2])))
        
        tilt_sweep_data.append({
            'tilt_angle_deg': float(ang_deg),
            'Ux_sim': float(U_t[0]),
            'Uz_sim': float(U_t[2]),
            'drift_angle_sim': float(drift_ang),
            'Ux_theory': float(U_th_vec[0]),
            'Uz_theory': float(U_th_vec[2]),
            'drift_angle_theory': float(drift_ang_th),
            'Omega_mag': float(np.linalg.norm(Om_t))
        })
        log(f"Tilt {ang_deg:4.1f}° | Ux = {U_t[0]:+.5f} | Uz = {U_t[2]:.5f} | Drift Angle = {drift_ang:5.2f}° | Theory Drift = {drift_ang_th:5.2f}° | |Omega| = {np.linalg.norm(Om_t):.2e}")
        
    tilt_csv_path = os.path.join(RAW_DATA_DIR, 'phase2a_disc_tilt_sweep.csv')
    with open(tilt_csv_path, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(tilt_sweep_data[0].keys()))
        writer.writeheader()
        writer.writerows(tilt_sweep_data)
    log(f"[Tilt Sweep Saved] -> {tilt_csv_path}")
    
    # -------------------------------------------------------------
    # Plot Generation
    # -------------------------------------------------------------
    generate_publication_plots(summary_by_resolution, tilt_sweep_data, PLOTS_DIR)
    log(f"[Plots Generated]  -> {PLOTS_DIR}")
    
    # Save log
    log_path = os.path.join(LOGS_DIR, 'phase2a_single_disc.log')
    with open(log_path, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    log(f"[Log Saved]        -> {log_path}")
    log("=" * 80)
    
    return summary_by_resolution, tilt_sweep_data

def generate_publication_plots(summary, tilt_data, output_dir):
    """Generates high-resolution publication plots."""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    N_vals = [s['N'] for s in summary]
    a_over_R = [s['a_blob_over_R'] for s in summary]
    Uz_face = [s['Uz_face'] for s in summary]
    Uz_edge = [s['Uz_edge'] for s in summary]
    ratios = [s['ratio_edge_to_face'] for s in summary]
    err_face = [s['rel_error_face_pct'] for s in summary]
    err_edge = [s['rel_error_edge_pct'] for s in summary]
    err_ratio = [s['rel_error_ratio_pct'] for s in summary]
    
    # 1. Velocities vs Resolution
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(N_vals, Uz_face, 'o-', color='#1f77b4', lw=2, ms=6, label=r'Simulated Face-on $U_\perp$')
    ax.axhline(1.0/16.0, color='#1f77b4', ls='--', alpha=0.7, label=r'Theory Face-on $U_\perp = 1/16$ (0.0625)')
    ax.plot(N_vals, Uz_edge, 's-', color='#d62728', lw=2, ms=6, label=r'Simulated Edge-on $U_\parallel$')
    ax.axhline(3.0/32.0, color='#d62728', ls='--', alpha=0.7, label=r'Theory Edge-on $U_\parallel = 3/32$ (0.09375)')
    ax.set_xlabel('Number of Blobs $N$', fontsize=12)
    ax.set_ylabel('Settling Speed $|U_z|$', fontsize=12)
    ax.set_title(r'Settling Speed Convergence Across Blob Resolutions ($R=1.0, F=1.0, \eta=1.0$)', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'disc_velocity_vs_resolution.png'))
    plt.close(fig)
    
    # 2. Anisotropy Ratio Convergence
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(N_vals, ratios, 'o-', color='#2ca02c', lw=2.2, ms=7, label=r'Simulated Ratio $U_\parallel / U_\perp$')
    ax.axhline(1.5, color='black', ls='--', lw=1.8, label=r'Theoretical Thin Disc Ratio = 1.5000')
    ax.set_xlabel('Number of Blobs $N$', fontsize=12)
    ax.set_ylabel(r'Velocity Ratio $U_\parallel / U_\perp$', fontsize=12)
    ax.set_title(r'Hydrodynamic Anisotropy Ratio Convergence ($U_\parallel / U_\perp \to 1.5$)', fontsize=12, fontweight='bold')
    ax.set_ylim(1.2, 1.55)
    ax.legend(frameon=True, fontsize=11)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'disc_anisotropy_ratio_convergence.png'))
    plt.close(fig)
    
    # 3. Relative Errors vs Blob Size (a_blob / R)
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(a_over_R, err_face, 'o-', color='#1f77b4', lw=2, ms=6, label=r'Face-on Error (%)')
    ax.plot(a_over_R, err_edge, 's-', color='#d62728', lw=2, ms=6, label=r'Edge-on Error (%)')
    ax.plot(a_over_R, err_ratio, '^-', color='#2ca02c', lw=2, ms=6, label=r'Ratio Error (%)')
    ax.set_xlabel(r'Normalized Blob Radius $a_{\rm blob} / R$', fontsize=12)
    ax.set_ylabel('Relative Error (%)', fontsize=12)
    ax.set_title(r'Convergence with Discretization Parameter $a_{\rm blob}/R \to 0$', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'disc_error_vs_blob_size.png'))
    plt.close(fig)
    
    # 4. Tilt Angle Sweep (Drift Angle vs Tilt Angle)
    tilts = [t['tilt_angle_deg'] for t in tilt_data]
    drift_sim = [t['drift_angle_sim'] for t in tilt_data]
    drift_th = [t['drift_angle_theory'] for t in tilt_data]
    Ux_sim = [t['Ux_sim'] for t in tilt_data]
    Uz_sim = [t['Uz_sim'] for t in tilt_data]
    
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(tilts, drift_sim, 'o-', color='#9467bd', lw=2, ms=5, label=r'Multiblob Simulation ($N=755$)')
    ax.plot(tilts, drift_th, '--', color='#ff7f0e', lw=2, label=r'Analytical Stokes Theory (Thin Disc)')
    ax.set_xlabel(r'Disc Normal Tilt Angle $\theta$ (degrees)', fontsize=12)
    ax.set_ylabel(r'Sedimentation Drift Angle $\phi$ (degrees)', fontsize=12)
    ax.set_title(r'Horizontal Drift Angle vs Initial Disc Tilt $\theta$', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=11)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'disc_drift_angle_vs_tilt.png'))
    plt.close(fig)

if __name__ == '__main__':
    run_phase2a_resolution_study()
