"""
MTP_Experiments/scripts/phase2_ellipsoid_force_resolution.py

Phase 2: Sedimentation Velocity vs Applied Force Across Multiple Blob Resolutions
for a Prolate Spheroid (Ellipsoid) in Low-Reynolds-Number Stokes Flow.

Investigates:
1. Multi-Resolution Force Sweep: N in {12, 42, 162, 642, 2562} across F_z in {0.2, 0.5, 1.0, 2.0, 5.0, 10.0}.
2. Comparison against exact continuum Perrin (1934) theoretical slopes for both:
   - Vertical sedimentation speed |U_z|
   - Lateral oblique drift velocity U_x (at tilt angle theta = 45 deg)
3. Slope convergence analysis: demonstrates that fitted mobility slopes m_N converge toward m_Perrin as N increases.
4. Preserves all canonical MTP physics definitions and scientific language rules.
"""

import sys
import os
import time
import json
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Dynamic relative paths
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
EXPERIMENTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
WORKSPACE_DIR = os.path.abspath(os.path.join(EXPERIMENTS_DIR, '..'))
REPO_PATH = os.path.abspath(os.path.join(WORKSPACE_DIR, 'RigidMultiblobsWall'))

if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)
if REPO_PATH not in sys.path:
    sys.path.insert(0, REPO_PATH)

import phase2_common as p2c
from quaternion_integrator.quaternion import Quaternion

# Configure output directories
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'ellipsoid')
RAW_DIR = os.path.join(OUT_DIR, 'raw_data')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')
REP_DIR = os.path.join(OUT_DIR, 'reports')

for d in [RAW_DIR, PROC_DIR, PLOT_DIR, REP_DIR]:
    os.makedirs(d, exist_ok=True)

# Styling for publication plots (300 DPI)
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9.5,
    'figure.titlesize': 14,
    'lines.linewidth': 2.0,
    'lines.markersize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

RESOLUTIONS = [12, 42, 162, 642, 2562]
FORCES = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
ASPECT_RATIO = 2.0
SEMI_MAJOR = 2.0
SEMI_MINOR = 1.0
ETA = 1.0
THETA_DEG = 45.0


def run_force_resolution_sweep():
    print("\n" + "="*75)
    print(f"PHASE 2: ELLIPSOID FORCE VS SEDIMENTATION VELOCITY ACROSS RESOLUTIONS")
    print(f"Aspect Ratio lambda = {ASPECT_RATIO:.1f} (a={SEMI_MAJOR}, b={SEMI_MINOR}), Tilt theta = {THETA_DEG} deg")
    print(f"Resolutions N in {RESOLUTIONS}, Forces F_z in {FORCES}")
    print("="*75 + "\n")

    # 1. Compute exact analytical Perrin continuum benchmark
    perrin = p2c.perrin_analytical_prolate(SEMI_MAJOR, SEMI_MINOR, eta=ETA)
    mu_par_th = perrin['mu_parallel']
    mu_perp_th = perrin['mu_perp']

    theta_rad = np.radians(THETA_DEG)
    # Exact continuum slopes at theta = 45 deg:
    # U_z = -F_z * [mu_perp + (mu_par - mu_perp) * sin^2(theta)]
    # U_x = +0.5 * F_z * (mu_par - mu_perp) * sin(2*theta)
    slope_z_th = mu_perp_th + (mu_par_th - mu_perp_th) * (np.sin(theta_rad)**2)
    slope_x_th = 0.5 * (mu_par_th - mu_perp_th) * np.sin(2.0 * theta_rad)

    print(f"--- Perrin Continuum Theory Reference ---")
    print(f"mu_parallel (Perrin): {mu_par_th:.8f}")
    print(f"mu_perp     (Perrin): {mu_perp_th:.8f}")
    print(f"Anisotropy Ratio:     {mu_par_th / mu_perp_th:.4f}")
    print(f"Theoretical Slope |U_z| / F_z: {slope_z_th:.8f}")
    print(f"Theoretical Slope  U_x  / F_z: {slope_x_th:.8f}\n")

    # Rotation quaternion for tilt angle theta
    q = Quaternion.from_rotation(np.array([0.0, theta_rad, 0.0]))

    all_data = []
    resolution_summaries = {}

    for N in RESOLUTIONS:
        t0 = time.time()
        r_conf, a_blob, meta = p2c.generate_ellipsoid_mesh(N, SEMI_MAJOR, SEMI_MINOR)

        # Precompute body mobility tensor N_body at this orientation
        # This allows instantaneous exact evaluation across all forces
        b = p2c.Body(np.zeros(3), q, r_conf, a_blob)
        K = b.calc_K_matrix()
        r_lab = b.get_r_vectors()
        M = p2c.mb.rotne_prager_tensor(r_lab, ETA, a_blob)
        L, lower = p2c.scipy.linalg.cho_factor(M)
        Kt_Minv_K = np.dot(K.T, p2c.scipy.linalg.cho_solve((L, lower), K, check_finite=False))
        N_body = np.linalg.pinv(Kt_Minv_K)
        t_factor = time.time() - t0

        Uz_list = []
        Ux_list = []
        F_list = []

        for Fz in FORCES:
            force_vec = np.array([0.0, 0.0, -Fz])
            wrench = np.concatenate([force_vec, np.zeros(3)])
            U_Omega = np.dot(N_body, wrench)
            U = U_Omega[0:3]
            Omega = U_Omega[3:6]

            row = {
                'N': N,
                'a_blob': a_blob,
                'Fz': Fz,
                'Ux': U[0],
                'Uy': U[1],
                'Uz': U[2],
                'abs_Uz': abs(U[2]),
                'U_mag': np.linalg.norm(U),
                'Omega_mag': np.linalg.norm(Omega),
                'Ux_th': slope_x_th * Fz,
                'abs_Uz_th': slope_z_th * Fz
            }
            all_data.append(row)
            Uz_list.append(abs(U[2]))
            Ux_list.append(U[0])
            F_list.append(Fz)

        # Linear regressions
        reg_z = stats.linregress(F_list, Uz_list)
        reg_x = stats.linregress(F_list, Ux_list)

        slope_z_num = reg_z.slope
        slope_x_num = reg_x.slope
        err_slope_z = 100.0 * abs(slope_z_num - slope_z_th) / slope_z_th
        err_slope_x = 100.0 * abs(slope_x_num - slope_x_th) / slope_x_th

        res_summary = {
            'N': N,
            'blob_radius': a_blob,
            't_factor_s': t_factor,
            'slope_z_num': slope_z_num,
            'slope_z_th': slope_z_th,
            'err_slope_z_pct': err_slope_z,
            'R2_z': reg_z.rvalue**2,
            'slope_x_num': slope_x_num,
            'slope_x_th': slope_x_th,
            'err_slope_x_pct': err_slope_x,
            'R2_x': reg_x.rvalue**2
        }
        resolution_summaries[N] = res_summary

        print(f"N = {N:3d} (a_blob={a_blob:.4f}, solve={t_factor:6.4f}s):")
        print(f"  |U_z| Slope = {slope_z_num:.6f} vs Theory = {slope_z_th:.6f} | Error = {err_slope_z:5.2f}% | R^2 = {reg_z.rvalue**2:.8f}")
        print(f"   U_x  Slope = {slope_x_num:.6f} vs Theory = {slope_x_th:.6f} | Error = {err_slope_x:5.2f}% | R^2 = {reg_x.rvalue**2:.8f}")

    # Save detailed CSV
    csv_file = os.path.join(PROC_DIR, "ellipsoid_force_resolution_data.csv")
    headers = ["N", "blob_radius", "Fz", "Ux", "Uy", "Uz", "abs_Uz", "U_mag", "Omega_mag", "Ux_th", "abs_Uz_th"]
    with open(csv_file, 'w') as f:
        f.write(",".join(headers) + "\n")
        for d in all_data:
            f.write(f"{d['N']},{d['a_blob']:.6e},{d['Fz']:.2f},{d['Ux']:.8e},{d['Uy']:.8e},{d['Uz']:.8e},"
                    f"{d['abs_Uz']:.8e},{d['U_mag']:.8e},{d['Omega_mag']:.8e},{d['Ux_th']:.8e},{d['abs_Uz_th']:.8e}\n")
    print(f"\nSaved detailed dataset to: {csv_file}")

    # Save summary JSON
    json_file = os.path.join(PROC_DIR, "ellipsoid_force_resolution_summary.json")
    with open(json_file, 'w') as f:
        json.dump(resolution_summaries, f, indent=2)
    print(f"Saved summary JSON to: {json_file}")

    # =========================================================================
    # PUBLICATION PLOT: Multi-Panel Force vs Sedimentation Velocity
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    colors = {12: '#2ca02c', 42: '#ff7f0e', 162: '#1f77b4', 642: '#9467bd', 2562: '#d62728'}
    markers = {12: 'o', 42: 's', 162: '^', 642: 'D', 2562: 'v'}

    F_dense = np.linspace(0.0, 10.0, 100)
    F_plot = np.array(FORCES)

    # Subplot 1: Vertical Sedimentation Speed |U_z| vs Applied Force F_z
    for N in RESOLUTIONS:
        pts = [d['abs_Uz'] for d in all_data if d['N'] == N]
        s = resolution_summaries[N]
        label_txt = rf'$N = {N}$ ($m = {s["slope_z_num"]:.4f}$, err: {s["err_slope_z_pct"]:.1f}%)'
        ax1.plot(F_plot, pts, marker=markers[N], color=colors[N], linestyle='--', label=label_txt)

    # Theory curve
    ax1.plot(F_dense, slope_z_th * F_dense, 'k-', linewidth=2.5,
             label=rf'Perrin Theory ($m = {slope_z_th:.4f}$)')

    ax1.set_xlabel(r'Applied Sedimentation Load $F_z$')
    ax1.set_ylabel(r'Downward Sedimentation Speed $|U_z|$')
    ax1.set_title(r'(a) Vertical Sedimentation Speed vs. Load ($\theta = 45^\circ$)')
    ax1.set_xlim(0, 10.5)
    ax1.set_ylim(0, 0.45)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper left', framealpha=0.92)

    # Subplot 2: Oblique Lateral Drift Velocity U_x vs Applied Force F_z
    for N in RESOLUTIONS:
        pts = [d['Ux'] for d in all_data if d['N'] == N]
        s = resolution_summaries[N]
        label_txt = rf'$N = {N}$ ($m = {s["slope_x_num"]:.5f}$, err: {s["err_slope_x_pct"]:.1f}%)'
        ax2.plot(F_plot, pts, marker=markers[N], color=colors[N], linestyle='--', label=label_txt)

    # Theory curve
    ax2.plot(F_dense, slope_x_th * F_dense, 'k-', linewidth=2.5,
             label=rf'Perrin Theory ($m = {slope_x_th:.5f}$)')

    ax2.set_xlabel(r'Applied Sedimentation Load $F_z$')
    ax2.set_ylabel(r'Lateral Drift Velocity $U_x$')
    ax2.set_title(r'(b) Lateral Oblique Drift Velocity vs. Load ($\theta = 45^\circ$)')
    ax2.set_xlim(0, 10.5)
    ax2.set_ylim(0, 0.030)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper left', framealpha=0.92)

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, "ellipsoid_force_vs_sedimentation_by_resolution.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Generated publication figure: {plot_file}\n")

    return resolution_summaries, all_data


if __name__ == '__main__':
    t_start = time.time()
    run_force_resolution_sweep()
    print(f"Complete in {time.time() - t_start:.2f} seconds!")
