"""
MTP_Experiments/scripts/phase2b_two_disc_sweep.py

Phase 2B: Two-Disc Hydrodynamic Force-Velocity Validation across N in {12, 42, 162, 642, 2562}.
1. Simulates two identical rigid discs falling under gravity in unbounded Stokes fluid.
2. Computes the full two-body mobility matrix including self and mutual interactions.
3. Conducts force sweeps F in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0} across all five resolutions.
4. Verifies symmetry (U1 approx U2), linearity (U proportional to F), and Re << 1.
5. Performs a separation sanity check across S/R in {5, 10, 20, 40, 80, 160}.
6. Generates the REQUIRED consolidated Velocity vs Force figure and clean CSV output.
"""

import sys
import os
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts.common import REPO_PATH, EXPERIMENTS_DIR
from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from read_input.read_vertex_file import read_vertex_file
import scipy.linalg

from scripts.disc_theory import thin_disc_velocity_perp

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
FORCES = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
PRIMARY_SEPARATION = 40.0  # S / R = 40.0
R_GEOM = 1.0
ETA = 1.0
RHO = 1.0

def load_disc_geometry(N):
    """Loads disc coordinates and blob radius for given N."""
    filepath = os.path.join(STRUCTURES_DIR, f"disc_N_{N}_R_{R_GEOM:.1f}.vertex")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Missing disc vertex file: {filepath}")
    r_conf = read_vertex_file(filepath)
    with open(filepath, 'r') as fp:
        for line in fp:
            if not line.startswith('#'):
                a_blob = float(line.split()[1]) / 2.0
                break
    return r_conf, a_blob

def solve_two_disc_system(r_conf, a_blob, separation=PRIMARY_SEPARATION, eta=ETA):
    """
    Assembles the 12x12 body mobility matrix for two identical discs separated by S along x-axis.
    Returns (N_body_12x12, runtime).
    """
    t0 = time.time()
    N_blobs = len(r_conf)
    
    pos1 = np.array([-separation / 2.0, 0.0, 0.0], dtype=np.float64)
    pos2 = np.array([+separation / 2.0, 0.0, 0.0], dtype=np.float64)
    q = Quaternion([1.0, 0.0, 0.0, 0.0])  # Face-on
    
    b1 = Body(pos1, q, r_conf, a_blob)
    b2 = Body(pos2, q, r_conf, a_blob)
    
    r_lab1 = b1.get_r_vectors()
    r_lab2 = b2.get_r_vectors()
    r_lab_total = np.vstack([r_lab1, r_lab2])  # Shape (2N, 3)
    
    K1 = b1.calc_K_matrix()  # Shape (3N, 6)
    K2 = b2.calc_K_matrix()  # Shape (3N, 6)
    
    # Block-diagonal geometric matrix K: (6N, 12)
    K_total = np.zeros((6 * N_blobs, 12), dtype=np.float64)
    K_total[0:3*N_blobs, 0:6] = K1
    K_total[3*N_blobs:6*N_blobs, 6:12] = K2
    
    # 6N x 6N RPY fluid mobility matrix
    M_total = mb.rotne_prager_tensor(r_lab_total, eta, a_blob)
    
    # Cholesky solve
    L, lower = scipy.linalg.cho_factor(M_total)
    Kt_Minv_K = np.dot(K_total.T, scipy.linalg.cho_solve((L, lower), K_total, check_finite=False))
    N_body_12x12 = np.linalg.pinv(Kt_Minv_K)
    
    runtime = time.time() - t0
    return N_body_12x12, runtime

def run_phase2b_validation():
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 85)
    log("PHASE 2B: TWO-DISC HYDRODYNAMIC FORCE-VELOCITY VALIDATION")
    log("=" * 85)
    log(f"Refinement Hierarchy: N in {RESOLUTIONS}")
    log(f"Force Sweep         : F in {FORCES}")
    log(f"Primary Separation  : S = {PRIMARY_SEPARATION} R (S/R = {PRIMARY_SEPARATION})")
    log(f"Orientation         : Face-on (broadside-on, normal along z)")
    log(f"Theoretical Law     : U_th = F / (16 * eta * R) = F / 16\n")
    
    results_records = []
    summary_by_resolution = {}
    
    for N in RESOLUTIONS:
        log("-" * 85)
        log(f">>> Processing Resolution N = {N:4d} Blobs <<<")
        r_conf, a_blob = load_disc_geometry(N)
        log(f"  Loaded geometry: a_blob = {a_blob:.6f}, N_total_blobs = {2*N}")
        
        # Assemble two-disc mobility matrix
        N_body_12x12, t_setup = solve_two_disc_system(r_conf, a_blob, separation=PRIMARY_SEPARATION, eta=ETA)
        log(f"  Mobility system assembled and inverted in {t_setup:.3f} s")
        
        res_data = []
        for F_val in FORCES:
            # Applied gravitational force on both discs along -z
            wrench_12 = np.zeros(12, dtype=np.float64)
            wrench_12[2] = -F_val  # Disc 1 Fz
            wrench_12[8] = -F_val  # Disc 2 Fz
            
            U_Omega = np.dot(N_body_12x12, wrench_12)
            U1_vec = U_Omega[0:3]
            Om1_vec = U_Omega[3:6]
            U2_vec = U_Omega[6:9]
            Om2_vec = U_Omega[9:12]
            
            U1 = abs(U1_vec[2])
            U2 = abs(U2_vec[2])
            U_avg = (U1 + U2) / 2.0
            diff_U1_U2 = abs(U1 - U2)
            
            U_th = thin_disc_velocity_perp(F_val, R_GEOM, eta=ETA)
            rel_error = abs(U_avg - U_th) / U_th * 100.0
            
            Re = RHO * U_avg * (2.0 * R_GEOM) / ETA
            
            record = {
                'N': int(N),
                'force': float(F_val),
                'velocity_disc_1': float(U1),
                'velocity_disc_2': float(U2),
                'velocity_average': float(U_avg),
                'theory_velocity': float(U_th),
                'relative_error': float(rel_error),
                'diff_U1_U2': float(diff_U1_U2),
                'separation': float(PRIMARY_SEPARATION),
                'orientation': 'face_on',
                'Re': float(Re),
                'Omega1_mag': float(np.linalg.norm(Om1_vec)),
                'Omega2_mag': float(np.linalg.norm(Om2_vec))
            }
            results_records.append(record)
            res_data.append(record)
            
            log(f"  F = {F_val:4.1f} | U1 = {U1:.6f} | U2 = {U2:.6f} | diff = {diff_U1_U2:.2e} | U_th = {U_th:.6f} | err = {rel_error:5.2f}% | Re = {Re:.4f}")
            
        summary_by_resolution[N] = {
            'a_blob': float(a_blob),
            'setup_time_s': float(t_setup),
            'mean_error_pct': float(np.mean([r['relative_error'] for r in res_data])),
            'max_diff_U1_U2': float(np.max([r['diff_U1_U2'] for r in res_data])),
            'fitted_mobility_slope': float(np.polyfit(FORCES, [r['velocity_average'] for r in res_data], 1)[0])
        }

    # Save CSV
    csv_file = os.path.join(RAW_DATA_DIR, 'phase2b_two_disc_force_velocity.csv')
    with open(csv_file, 'w', newline='') as fp:
        fieldnames = [
            'N', 'force', 'velocity_disc_1', 'velocity_disc_2', 'velocity_average',
            'theory_velocity', 'relative_error', 'separation', 'orientation', 'Re'
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(results_records)
    log(f"\n[Raw CSV Saved] -> {csv_file}")
    
    # Save JSON Summary
    json_file = os.path.join(PROCESSED_DATA_DIR, 'phase2b_two_disc_summary.json')
    with open(json_file, 'w') as fp:
        json.dump(summary_by_resolution, fp, indent=2)
    log(f"[Summary Saved] -> {json_file}")
    
    # -------------------------------------------------------------
    # Separation Sanity Check
    # -------------------------------------------------------------
    log("\n" + "=" * 85)
    log("SEPARATION SANITY CHECK (S / R sweep for N = 162, F = 1.0)")
    log("=" * 85)
    separations = [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]
    r_conf_162, a_blob_162 = load_disc_geometry(162)
    
    # Isolated single disc velocity at N=162
    b_single = Body(np.zeros(3), Quaternion([1,0,0,0]), r_conf_162, a_blob_162)
    M_s = mb.rotne_prager_tensor(b_single.get_r_vectors(), ETA, a_blob_162)
    K_s = b_single.calc_K_matrix()
    Ls, lows = scipy.linalg.cho_factor(M_s)
    N_s = np.linalg.pinv(np.dot(K_s.T, scipy.linalg.cho_solve((Ls, lows), K_s)))
    U_single_162 = abs(np.dot(N_s, [0, 0, -1.0, 0, 0, 0])[2])
    log(f"Isolated single disc benchmark at N=162: U_single = {U_single_162:.6f}")
    
    sep_records = []
    for S_val in separations:
        N_body_sep, _ = solve_two_disc_system(r_conf_162, a_blob_162, separation=S_val, eta=ETA)
        w_sep = np.zeros(12)
        w_sep[2] = -1.0
        w_sep[8] = -1.0
        U_sep = abs(np.dot(N_body_sep, w_sep)[2])
        diff_from_single = abs(U_sep - U_single_162)
        rel_diff_pct = diff_from_single / U_single_162 * 100.0
        
        sep_records.append({
            'separation': float(S_val),
            'U_two_disc': float(U_sep),
            'U_single_disc': float(U_single_162),
            'diff_from_single': float(diff_from_single),
            'rel_diff_pct': float(rel_diff_pct)
        })
        log(f"  S = {S_val:5.1f} R | U = {U_sep:.6f} | diff vs single = {diff_from_single:+.2e} ({rel_diff_pct:5.2f}%)")
        
    sep_csv = os.path.join(RAW_DATA_DIR, 'phase2b_separation_check.csv')
    with open(sep_csv, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(sep_records[0].keys()))
        writer.writeheader()
        writer.writerows(sep_records)
    log(f"[Separation CSV Saved] -> {sep_csv}")
    
    # -------------------------------------------------------------
    # Plotting
    # -------------------------------------------------------------
    generate_phase2b_plots(results_records, summary_by_resolution, sep_records, PLOTS_DIR)
    log(f"[Plots Generated]      -> {PLOTS_DIR}")
    
    # Save Log
    log_file = os.path.join(LOGS_DIR, 'phase2b_two_disc_force_velocity.log')
    with open(log_file, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    log(f"[Log Saved]            -> {log_file}")
    log("=" * 85)
    
    return results_records, summary_by_resolution, sep_records

def generate_phase2b_plots(results, summary, sep_data, output_dir):
    """Generates the required ONE consolidated plot and diagnostic figures."""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # =========================================================================
    # REQUIRED MAIN DELIVERABLE: ONE Consolidated Velocity vs Force Plot
    # =========================================================================
    fig, ax = plt.subplots(figsize=(8.5, 6.0), dpi=300)
    
    # Theoretical line: U = F / 16 (slope = 0.0625)
    F_dense = np.linspace(0.0, 5.2, 100)
    U_th_dense = F_dense / 16.0
    ax.plot(F_dense, U_th_dense, 'k--', lw=2.2, label=r'Analytical Stokes Theory: $U = F / (16 \eta R)$ (slope $= 0.0625$)')
    
    # Colors and markers for the 5 refinement levels
    styles = {
        12:   {'color': '#9467bd', 'marker': '^', 'label': r'$N = 12$ ($a_{\rm blob}/R = 0.312$)'},
        42:   {'color': '#1f77b4', 'marker': 's', 'label': r'$N = 42$ ($a_{\rm blob}/R = 0.167$)'},
        162:  {'color': '#2ca02c', 'marker': 'o', 'label': r'$N = 162$ ($a_{\rm blob}/R = 0.085$)'},
        642:  {'color': '#ff7f0e', 'marker': 'D', 'label': r'$N = 642$ ($a_{\rm blob}/R = 0.043$)'},
        2562: {'color': '#d62728', 'marker': 'v', 'label': r'$N = 2562$ ($a_{\rm blob}/R = 0.021$)'}
    }
    
    for N_val in RESOLUTIONS:
        sub = [r for r in results if r['N'] == N_val]
        F_arr = [r['force'] for r in sub]
        U_arr = [r['velocity_average'] for r in sub]
        st = styles[N_val]
        ax.plot(F_arr, U_arr, marker=st['marker'], color=st['color'], lw=1.8, ms=6, label=st['label'])
        
    ax.set_xlabel('Applied Force $F$ (code units)', fontsize=13, labelpad=8)
    ax.set_ylabel('Settling Velocity $U$ (code units)', fontsize=13, labelpad=8)
    ax.set_title(r'Consolidated Velocity vs Applied Force across Resolutions ($N \in \{12, 42, 162, 642, 2562\}$)',
                 fontsize=13, fontweight='bold', pad=12)
    ax.set_xlim(-0.1, 5.2)
    ax.set_ylim(-0.01, 0.34)
    ax.legend(frameon=True, fontsize=10, loc='upper left', framealpha=0.95)
    ax.grid(True, ls=':', alpha=0.6)
    
    # Annotate linearity and convergence
    ax.text(3.5, 0.08, 'Creeping-Flow Regime ($\\mathrm{Re} \\ll 1$)\nStrict Linear Response: $U \\propto F$',
            fontsize=10, bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0', edgecolor='#cccccc', alpha=0.9))
    
    plt.tight_layout()
    main_plot_path = os.path.join(output_dir, 'two_disc_velocity_vs_force_consolidated.png')
    fig.savefig(main_plot_path)
    plt.close(fig)
    print(f"  Saved main consolidated plot: {main_plot_path}")
    
    # =========================================================================
    # Diagnostic Plot 1: Relative Error vs Resolution N
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    
    N_list = RESOLUTIONS
    mean_errors = [summary[n]['mean_error_pct'] for n in N_list]
    slopes = [summary[n]['fitted_mobility_slope'] for n in N_list]
    th_slope = 1.0 / 16.0
    
    ax1.plot(N_list, mean_errors, 'o-', color='#d62728', lw=2, ms=6)
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of Blobs $N$ (log scale)', fontsize=12)
    ax1.set_ylabel('Mean Relative Error vs Theory (%)', fontsize=12)
    ax1.set_title('Resolution Convergence: Error vs $N$', fontsize=12, fontweight='bold')
    ax1.grid(True, which='both', ls=':', alpha=0.6)
    
    ax2.plot(N_list, slopes, 's-', color='#1f77b4', lw=2, ms=6, label='Fitted Slope $m_N$')
    ax2.axhline(th_slope, color='k', ls='--', label=r'Theory Slope $1/16 = 0.0625$')
    ax2.set_xscale('log')
    ax2.set_xlabel('Number of Blobs $N$ (log scale)', fontsize=12)
    ax2.set_ylabel(r'Mobility Slope $dU/dF$', fontsize=12)
    ax2.set_title(r'Convergence of Mobility Slope $dU/dF \to 1/16$', fontsize=12, fontweight='bold')
    ax2.legend(frameon=True, fontsize=10)
    ax2.grid(True, which='both', ls=':', alpha=0.6)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'two_disc_error_vs_resolution.png'))
    plt.close(fig)
    
    # =========================================================================
    # Diagnostic Plot 2: Separation Sanity Check
    # =========================================================================
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    sep_arr = [r['separation'] for r in sep_data]
    u_two = [r['U_two_disc'] for r in sep_data]
    u_sing = sep_data[0]['U_single_disc']
    
    ax.plot(sep_arr, u_two, 'o-', color='#2ca02c', lw=2, ms=6, label=r'Two-Disc Velocity $U(S)$ ($N=162, F=1.0$)')
    ax.axhline(u_sing, color='k', ls='--', lw=1.8, label=f'Isolated Single-Disc Limit ($U = {u_sing:.6f}$)')
    ax.set_xlabel('Inter-Disc Separation $S / R$', fontsize=12)
    ax.set_ylabel('Settling Velocity $U$', fontsize=12)
    ax.set_title('Separation Sanity Check: Recovery of Isolated Limit ($S/R \\to \\infty$)', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'two_disc_separation_check.png'))
    plt.close(fig)

if __name__ == '__main__':
    run_phase2b_validation()
