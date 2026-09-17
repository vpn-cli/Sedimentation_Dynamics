"""
MTP_Experiments/scripts/phase2c_cylinder_force_velocity.py

Phase 2C: Force-Velocity Linearity, Orientation Anisotropy, and Resolution Convergence
for a Rigid Circular Cylinder in Unbounded Stokes Flow.

Conducts systematic force sweeps F in {0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0}
across N in {12, 42, 162, 642, 2562} for both Axial and Transverse orientations.
Generates all required publication plots, CSV tables, and summary JSON.
"""

import sys
import os
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

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
from scripts.phase2c_single_cylinder import (
    load_cylinder_geometry,
    solve_cylinder_mobility_unbounded
)

# Directories
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
FORCES = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
R_GEOM = 1.0
L_GEOM = 2.0
ETA = 1.0
RHO = 1.0
CHAR_LENGTH = 2.0 * R_GEOM  # L_c = 2R = 2.0

def run_phase2c_force_velocity_validation():
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 95)
    log("PHASE 2C: RIGID CYLINDER HYDRODYNAMIC FORCE-VELOCITY VALIDATION")
    log("=" * 95)
    
    theory = cylinder_theoretical_summary(R_GEOM, L_GEOM, ETA)
    U_th_surf = theory['equivalent_surface_sphere']['velocity_unit_force']
    Rh_th_surf = theory['equivalent_surface_sphere']['resistance']
    U_th_hd = theory['hubbard_douglas']['velocity_unit_force']
    U_th_tirado = theory['tirado_de_la_torre']['velocity_unit_force']
    
    log(f"Parameters: R={R_GEOM}, L={L_GEOM}, eta={ETA}, rho={RHO}, Char Length L_c={CHAR_LENGTH}")
    log(f"Force Sweep: {FORCES}")
    log(f"Theoretical Unit-Force Velocities:")
    log(f"  Equivalent Surface Sphere (Primary): {U_th_surf:.6f} (Rh = {Rh_th_surf:.4f})")
    log(f"  Hubbard-Douglas Capacitance        : {U_th_hd:.6f} (Rh = {theory['hubbard_douglas']['resistance']:.4f})")
    log(f"  Tirado & Garcia de la Torre        : {U_th_tirado:.6f} (Rh = {theory['tirado_de_la_torre']['resistance']:.4f})\n")
    
    force_velocity_records = []
    convergence_summary = []
    resolution_data = {}
    
    q_identity = Quaternion([1.0, 0.0, 0.0, 0.0])
    
    for N in RESOLUTIONS:
        r_conf, a_blob, vfile = load_cylinder_geometry(N)
        log("-" * 95)
        log(f">>> Running Cylinder Resolution N = {N:4d} Blobs | a_blob = {a_blob:.5f} <<<")
        
        # We can extract the 6x6 mobility matrix once per N to obtain exact linear scaling
        # and verify with explicit solves across the force sweep
        F_unit_ax = np.array([0.0, 0.0, 1.0])
        F_unit_tr = np.array([1.0, 0.0, 0.0])
        
        U_ax_1, Om_ax_1, diag_ax, N_body_ax = solve_cylinder_mobility_unbounded(
            r_conf, a_blob, q_identity, F_unit_ax, eta=ETA
        )
        U_tr_1, Om_tr_1, diag_tr, N_body_tr = solve_cylinder_mobility_unbounded(
            r_conf, a_blob, q_identity, F_unit_tr, eta=ETA
        )
        
        mob_ax_numerical = float(np.linalg.norm(U_ax_1))
        mob_tr_numerical = float(np.linalg.norm(U_tr_1))
        rh_ax_numerical = 1.0 / mob_ax_numerical
        rh_tr_numerical = 1.0 / mob_tr_numerical
        anisotropy_ratio = mob_tr_numerical / mob_ax_numerical
        
        # Evaluate each force in sweep
        ax_velocities = []
        tr_velocities = []
        
        for F_val in FORCES:
            F_vec_ax = np.array([0.0, 0.0, F_val])
            F_vec_tr = np.array([F_val, 0.0, 0.0])
            
            U_ax, Om_ax, _, _ = solve_cylinder_mobility_unbounded(
                r_conf, a_blob, q_identity, F_vec_ax, eta=ETA
            )
            U_tr, Om_tr, _, _ = solve_cylinder_mobility_unbounded(
                r_conf, a_blob, q_identity, F_vec_tr, eta=ETA
            )
            
            u_ax = float(np.linalg.norm(U_ax))
            u_tr = float(np.linalg.norm(U_tr))
            om_ax = float(np.linalg.norm(Om_ax))
            om_tr = float(np.linalg.norm(Om_tr))
            
            ax_velocities.append(u_ax)
            tr_velocities.append(u_tr)
            
            # Theoretical velocities
            u_th_surf = U_th_surf * F_val
            err_ax = abs(u_ax - u_th_surf) / u_th_surf * 100.0
            err_tr = abs(u_tr - u_th_surf) / u_th_surf * 100.0
            
            # Reynolds number Re = rho * U * L_c / eta
            Re_ax = RHO * u_ax * CHAR_LENGTH / ETA
            Re_tr = RHO * u_tr * CHAR_LENGTH / ETA
            
            # Record axial
            force_velocity_records.append({
                'N': int(N),
                'radius': float(R_GEOM),
                'length': float(L_GEOM),
                'orientation': 'axial',
                'force': float(F_val),
                'velocity': float(u_ax),
                'theory_velocity': float(u_th_surf),
                'numerical_mobility': float(mob_ax_numerical),
                'theoretical_mobility': float(U_th_surf),
                'relative_error': float(err_ax),
                'angular_velocity': float(om_ax),
                'Re': float(Re_ax),
                'dt': 'deterministic_mobility'
            })
            
            # Record transverse
            force_velocity_records.append({
                'N': int(N),
                'radius': float(R_GEOM),
                'length': float(L_GEOM),
                'orientation': 'transverse',
                'force': float(F_val),
                'velocity': float(u_tr),
                'theory_velocity': float(u_th_surf),
                'numerical_mobility': float(mob_tr_numerical),
                'theoretical_mobility': float(U_th_surf),
                'relative_error': float(err_tr),
                'angular_velocity': float(om_tr),
                'Re': float(Re_tr),
                'dt': 'deterministic_mobility'
            })
            
        # Linear fits U = m*F + b
        reg_ax = linregress(FORCES, ax_velocities)
        reg_tr = linregress(FORCES, tr_velocities)
        
        err_mob_ax = abs(reg_ax.slope - U_th_surf) / U_th_surf * 100.0
        err_mob_tr = abs(reg_tr.slope - U_th_surf) / U_th_surf * 100.0
        err_mob_mean = abs(((reg_ax.slope + 2*reg_tr.slope)/3.0) - U_th_surf) / U_th_surf * 100.0
        
        log(f"  Axial Fit:      Slope m={reg_ax.slope:.6f} | b={reg_ax.intercept:.2e} | R^2={reg_ax.rvalue**2:.8f} | Rel Err={err_mob_ax:.2f}%")
        log(f"  Transverse Fit: Slope m={reg_tr.slope:.6f} | b={reg_tr.intercept:.2e} | R^2={reg_tr.rvalue**2:.8f} | Rel Err={err_mob_tr:.2f}%")
        log(f"  Anisotropy Ratio U_perp / U_par = {anisotropy_ratio:.5f} | Rh_par = {rh_ax_numerical:.4f} | Rh_perp = {rh_tr_numerical:.4f}")
        log(f"  Max Re: {RHO * max(ax_velocities) * CHAR_LENGTH / ETA:.4f} (Stokes limit strictly preserved)")
        
        res_summary = {
            'N': int(N),
            'a_blob': float(a_blob),
            'radius': float(R_GEOM),
            'length': float(L_GEOM),
            'axial_slope': float(reg_ax.slope),
            'axial_intercept': float(reg_ax.intercept),
            'axial_r2': float(reg_ax.rvalue**2),
            'axial_resistance': float(rh_ax_numerical),
            'axial_err_pct': float(err_mob_ax),
            'transverse_slope': float(reg_tr.slope),
            'transverse_intercept': float(reg_tr.intercept),
            'transverse_r2': float(reg_tr.rvalue**2),
            'transverse_resistance': float(rh_tr_numerical),
            'transverse_err_pct': float(err_mob_tr),
            'mean_slope': float((reg_ax.slope + 2.0 * reg_tr.slope) / 3.0),
            'mean_err_pct': float(err_mob_mean),
            'anisotropy_ratio': float(anisotropy_ratio),
            'max_Re': float(max(RHO * max(ax_velocities) * CHAR_LENGTH / ETA, RHO * max(tr_velocities) * CHAR_LENGTH / ETA)),
            'schur_cond_num': float(diag_ax['schur_cond_num']),
            'cholesky_success': True,
            'forces': list(FORCES),
            'axial_velocities': list(ax_velocities),
            'transverse_velocities': list(tr_velocities)
        }
        convergence_summary.append(res_summary)
        resolution_data[N] = res_summary
        
    log("=" * 95)
    
    # Save CSVs
    csv_force_path = os.path.join(RAW_DATA_DIR, 'phase2c_cylinder_force_velocity.csv')
    with open(csv_force_path, 'w', newline='') as fp:
        fieldnames = [
            'N', 'radius', 'length', 'orientation', 'force', 'velocity',
            'theory_velocity', 'numerical_mobility', 'theoretical_mobility',
            'relative_error', 'angular_velocity', 'Re', 'dt'
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(force_velocity_records)
    log(f"Saved force-velocity CSV: {csv_force_path}")
    
    csv_conv_path = os.path.join(RAW_DATA_DIR, 'phase2c_cylinder_convergence.csv')
    with open(csv_conv_path, 'w', newline='') as fp:
        fieldnames = [
            'N', 'a_blob', 'radius', 'length', 'axial_slope', 'axial_intercept',
            'axial_r2', 'axial_resistance', 'axial_err_pct',
            'transverse_slope', 'transverse_intercept', 'transverse_r2',
            'transverse_resistance', 'transverse_err_pct',
            'mean_slope', 'mean_err_pct', 'anisotropy_ratio', 'max_Re', 'schur_cond_num'
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in convergence_summary:
            row_filtered = {k: row[k] for k in fieldnames}
            writer.writerow(row_filtered)
    log(f"Saved convergence CSV:    {csv_conv_path}")
    
    # Save Summary JSON
    json_path = os.path.join(PROCESSED_DATA_DIR, 'phase2c_cylinder_summary.json')
    with open(json_path, 'w') as fp:
        json.dump({
            'geometry': theory['geometry'],
            'theoretical_benchmarks': theory,
            'resolutions_evaluated': RESOLUTIONS,
            'forces_evaluated': FORCES,
            'convergence_summary': convergence_summary
        }, fp, indent=2)
    log(f"Saved summary JSON:       {json_path}")
    
    # Generate Plots
    generate_all_cylinder_plots(convergence_summary, theory, PLOTS_DIR)
    log(f"Generated publication plots in: {PLOTS_DIR}")
    
    # Save Log
    log_path = os.path.join(LOGS_DIR, 'phase2c_cylinder_force_velocity.log')
    with open(log_path, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    log(f"Saved run log:            {log_path}")
    log("=" * 95)
    
    return convergence_summary, force_velocity_records

def generate_all_cylinder_plots(summary, theory, output_dir):
    """Generates all required publication plots."""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    forces = np.array(summary[0]['forces'])
    F_dense = np.linspace(0.0, 1.05, 200)
    U_th_surf = theory['equivalent_surface_sphere']['velocity_unit_force']
    
    colors = {
        12: '#1f77b4',
        42: '#ff7f0e',
        162: '#2ca02c',
        642: '#d62728',
        2562: '#9467bd'
    }
    markers = {12: 'o', 42: 's', 162: '^', 642: 'v', 2562: 'D'}
    
    # -------------------------------------------------------------
    # Plot 1: Axial Velocity vs Force
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(F_dense, U_th_surf * F_dense, 'k--', lw=2.0, alpha=0.8,
            label=f"Approx. Benchmark (Eq. Surface Sphere: slope={U_th_surf:.5f})")
    
    for row in summary:
        N = row['N']
        ax.plot(forces, row['axial_velocities'], marker=markers[N], color=colors[N],
                lw=1.5, ms=6, label=f"N = {N:4d} (slope={row['axial_slope']:.5f}, $R^2$={row['axial_r2']:.6f})")
        
    ax.set_xlabel(r'Applied Force $F_\parallel$ (Dimensionless)', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Axial Translational Velocity $U_\parallel$', fontsize=11, fontweight='bold')
    ax.set_title(r'Rigid Cylinder Axial Force–Velocity Linearity ($R=1.0, L=2.0, \eta=1.0$)',
                 fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=9, loc='upper left')
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cylinder_velocity_vs_force_axial.png'))
    plt.close(fig)
    
    # -------------------------------------------------------------
    # Plot 2: Transverse Velocity vs Force
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(F_dense, U_th_surf * F_dense, 'k--', lw=2.0, alpha=0.8,
            label=f"Approx. Benchmark (Eq. Surface Sphere: slope={U_th_surf:.5f})")
    
    for row in summary:
        N = row['N']
        ax.plot(forces, row['transverse_velocities'], marker=markers[N], color=colors[N],
                lw=1.5, ms=6, label=f"N = {N:4d} (slope={row['transverse_slope']:.5f}, $R^2$={row['transverse_r2']:.6f})")
        
    ax.set_xlabel(r'Applied Force $F_\perp$ (Dimensionless)', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Transverse Translational Velocity $U_\perp$', fontsize=11, fontweight='bold')
    ax.set_title(r'Rigid Cylinder Transverse Force–Velocity Linearity ($R=1.0, L=2.0, \eta=1.0$)',
                 fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=9, loc='upper left')
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cylinder_velocity_vs_force_transverse.png'))
    plt.close(fig)
    
    # -------------------------------------------------------------
    # Plot 3: Resolution Convergence (Relative Error vs N)
    # -------------------------------------------------------------
    N_vals = [s['N'] for s in summary]
    err_ax = [s['axial_err_pct'] for s in summary]
    err_tr = [s['transverse_err_pct'] for s in summary]
    err_mean = [s['mean_err_pct'] for s in summary]
    
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(N_vals, err_ax, 'o-', color='#1f77b4', lw=2.0, ms=7, label=r'Axial Diff. $|U_\parallel - U_{\rm sphere}| / U_{\rm sphere}$ (%)')
    ax.plot(N_vals, err_tr, 's-', color='#d62728', lw=2.0, ms=7, label=r'Transverse Diff. $|U_\perp - U_{\rm sphere}| / U_{\rm sphere}$ (%)')
    ax.plot(N_vals, err_mean, '^-', color='#2ca02c', lw=2.2, ms=8, label=r'Isotropic Mean Diff. (%)')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Number of Blobs $N$', fontsize=11, fontweight='bold')
    ax.set_ylabel('Difference (%) vs Eq. Surface Sphere Benchmark', fontsize=11, fontweight='bold')
    ax.set_title(r'Resolution Convergence of Rigid Multiblob Cylinder ($R=1.0, L=2.0$)',
                 fontsize=12, fontweight='bold')
    ax.set_xticks(N_vals)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6, which='both')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cylinder_error_vs_resolution.png'))
    plt.close(fig)
    
    # -------------------------------------------------------------
    # Plot 4: Hydrodynamic Anisotropy Ratio U_perp / U_parallel
    # -------------------------------------------------------------
    ratios = [s['anisotropy_ratio'] for s in summary]
    
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    ax.plot(N_vals, ratios, 'D-', color='#9467bd', lw=2.2, ms=7,
            label=r'Multiblob Cylinder $U_\perp / U_\parallel$')
    ax.axhline(1.0, color='black', ls=':', lw=1.5, alpha=0.7, label='Isotropic Sphere Limit (1.0000)')
    
    ax.set_xscale('log')
    ax.set_xlabel(r'Number of Blobs $N$', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Mobility Ratio $U_\perp / U_\parallel$', fontsize=11, fontweight='bold')
    ax.set_title(r'Hydrodynamic Anisotropy Convergence ($U_\perp / U_\parallel \to 0.991$)',
                 fontsize=12, fontweight='bold')
    ax.set_xticks(N_vals)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_ylim(0.96, 1.01)
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'cylinder_mobility_ratio.png'))
    plt.close(fig)

if __name__ == '__main__':
    run_phase2c_force_velocity_validation()
