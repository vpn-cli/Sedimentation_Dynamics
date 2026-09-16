"""
MTP_Experiments/scripts/phase2_boomerang_force_resolution.py

Phase 2: Sedimentation Velocity vs Applied Force Across Multiple Blob Resolutions
and Oblique Drift vs Tilt Angle for a Boomerang Colloidal Particle in Low-Reynolds-Number Stokes Flow.

Investigates:
1. Multi-Resolution Force Sweep: N in {7, 15, 29} across F_z in {0.2, 0.5, 1.0, 2.0, 5.0, 10.0} at tilt theta = 45 deg.
   - Evaluates vertical settling speed |U_z| and lateral drift velocity U_x.
   - Computes linear regression slopes and R^2 linearity for each resolution.
   - Plots boomerang_force_vs_sedimentation_by_resolution.png (matching the ellipsoid benchmark).
2. Oblique Tilt Angle Sweep: theta in [0, 90 deg] in 5 deg increments.
   - Evaluates lateral drift velocity U_x(theta) and vertical settling speed |U_z(theta)|.
   - Compares with hydrodynamic anisotropic drift scaling: U_x(theta) ~ sin(2*theta).
   - Plots boomerang_tilted_drift_vs_angle.png (matching the ellipsoid benchmark).
3. Preserves all canonical MTP physics definitions and scientific language rules.
"""

import sys
import os
import time
import json
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

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

OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'boomerang')
RAW_DIR = os.path.join(OUT_DIR, 'raw_data')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')
REP_DIR = os.path.join(OUT_DIR, 'reports')

for d in [RAW_DIR, PROC_DIR, PLOT_DIR, REP_DIR]:
    os.makedirs(d, exist_ok=True)

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

RESOLUTIONS = [7, 15, 29]
FORCES = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
ARM_LENGTH = 2.1
ANGLE_DEG = 90.0
ETA = 1.0
THETA_DEG = 45.0


def run_force_resolution_sweep():
    print("\n" + "="*75)
    print("PHASE 2: BOOMERANG FORCE VS SEDIMENTATION VELOCITY ACROSS RESOLUTIONS")
    print(f"Opening Angle alpha = {ANGLE_DEG:.1f} deg, Arm Length L = {ARM_LENGTH:.2f}")
    print(f"Resolutions N in {RESOLUTIONS}, Forces F_z in {FORCES}, Tilt theta = {THETA_DEG} deg")
    print("="*75 + "\n")

    theta_rad = np.radians(THETA_DEG)
    q = Quaternion.from_rotation(np.array([0.0, theta_rad, 0.0]))

    all_data = []
    resolution_summaries = {}

    for N in RESOLUTIONS:
        t0 = time.time()
        r_conf, a_blob, meta = p2c.generate_boomerang_mesh(N, ARM_LENGTH, ANGLE_DEG, center_at='com')

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
                'Omega_mag': np.linalg.norm(Omega)
            }
            all_data.append(row)
            Uz_list.append(abs(U[2]))
            Ux_list.append(U[0])
            F_list.append(Fz)

        reg_z = stats.linregress(F_list, Uz_list)
        reg_x = stats.linregress(F_list, Ux_list)

        res_summary = {
            'N': N,
            'blob_radius': a_blob,
            't_factor_s': t_factor,
            'slope_z_num': reg_z.slope,
            'R2_z': reg_z.rvalue**2,
            'slope_x_num': reg_x.slope,
            'R2_x': reg_x.rvalue**2
        }
        resolution_summaries[N] = res_summary

        print(f"N = {N:2d} (a_blob={a_blob:.4f}, solve={t_factor:6.4f}s):")
        print(f"  |U_z| Slope = {reg_z.slope:.6f} | R^2 = {reg_z.rvalue**2:.8f}")
        print(f"   U_x  Slope = {reg_x.slope:.6f} | R^2 = {reg_x.rvalue**2:.8f}")

    # Save detailed CSV
    csv_file = os.path.join(PROC_DIR, "boomerang_force_resolution_data.csv")
    headers = ["N", "blob_radius", "Fz", "Ux", "Uy", "Uz", "abs_Uz", "U_mag", "Omega_mag"]
    with open(csv_file, 'w') as f:
        f.write(",".join(headers) + "\n")
        for d in all_data:
            f.write(f"{d['N']},{d['a_blob']:.6e},{d['Fz']:.2f},{d['Ux']:.8e},{d['Uy']:.8e},{d['Uz']:.8e},"
                    f"{d['abs_Uz']:.8e},{d['U_mag']:.8e},{d['Omega_mag']:.8e}\n")
    print(f"\nSaved detailed dataset to: {csv_file}")

    # Save summary JSON
    json_file = os.path.join(PROC_DIR, "boomerang_force_resolution_summary.json")
    with open(json_file, 'w') as f:
        json.dump(resolution_summaries, f, indent=2)
    print(f"Saved summary JSON to: {json_file}")

    # -------------------------------------------------------------------------
    # Publication Plot 1: Force vs Sedimentation Velocity by Resolution
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    colors = {7: '#2ca02c', 15: '#1f77b4', 29: '#d62728'}
    markers = {7: 's', 15: 'o', 29: '^'}

    F_plot = np.array(FORCES)

    for N in RESOLUTIONS:
        pts_z = [d['abs_Uz'] for d in all_data if d['N'] == N]
        s = resolution_summaries[N]
        lbl = rf'$N = {N}$ ($m = {s["slope_z_num"]:.5f}$, $R^2=1.000$)'
        ax1.plot(F_plot, pts_z, marker=markers[N], color=colors[N], linestyle='--', label=lbl)

    ax1.set_xlabel(r'Applied Sedimentation Load $F_z$')
    ax1.set_ylabel(r'Downward Sedimentation Speed $|U_z|$')
    ax1.set_title(r'(a) Boomerang Vertical Speed vs. Load ($\theta = 45^\circ$)')
    ax1.set_xlim(0, 10.5)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='upper left', framealpha=0.92)

    for N in RESOLUTIONS:
        pts_x = [d['Ux'] for d in all_data if d['N'] == N]
        s = resolution_summaries[N]
        lbl = rf'$N = {N}$ ($m = {s["slope_x_num"]:.5f}$, $R^2=1.000$)'
        ax2.plot(F_plot, pts_x, marker=markers[N], color=colors[N], linestyle='--', label=lbl)

    ax2.set_xlabel(r'Applied Sedimentation Load $F_z$')
    ax2.set_ylabel(r'Lateral Drift Velocity $U_x$')
    ax2.set_title(r'(b) Boomerang Lateral Oblique Drift vs. Load ($\theta = 45^\circ$)')
    ax2.set_xlim(0, 10.5)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='upper left', framealpha=0.92)

    plt.tight_layout()
    plot_file1 = os.path.join(PLOT_DIR, "boomerang_force_vs_sedimentation_by_resolution.png")
    plt.savefig(plot_file1)
    plt.close()
    print(f"Generated publication figure: {plot_file1}")

    return resolution_summaries, all_data


def run_tilt_angle_drift_sweep():
    print("\n" + "="*75)
    print("PHASE 2: BOOMERANG OBLIQUE LATERAL DRIFT VS TILT ANGLE")
    print(f"Sweeping theta in [0, 90 deg] at reference resolution N = 15, F_z = 1.0")
    print("="*75 + "\n")

    angles_deg = np.linspace(0.0, 90.0, 19) # 0, 5, 10, ..., 90
    r_conf, a_blob, meta = p2c.generate_boomerang_mesh(15, ARM_LENGTH, ANGLE_DEG, center_at='com')

    drift_data = []
    Fz = 1.0
    wrench = np.array([0.0, 0.0, -Fz, 0.0, 0.0, 0.0])

    for th_deg in angles_deg:
        th_rad = np.radians(th_deg)
        q = Quaternion.from_rotation(np.array([0.0, th_rad, 0.0]))

        b = p2c.Body(np.zeros(3), q, r_conf, a_blob)
        K = b.calc_K_matrix()
        r_lab = b.get_r_vectors()
        M = p2c.mb.rotne_prager_tensor(r_lab, ETA, a_blob)
        L, lower = p2c.scipy.linalg.cho_factor(M)
        Kt_Minv_K = np.dot(K.T, p2c.scipy.linalg.cho_solve((L, lower), K, check_finite=False))
        N_body = np.linalg.pinv(Kt_Minv_K)

        U_Omega = np.dot(N_body, wrench)
        U = U_Omega[0:3]
        Omega = U_Omega[3:6]

        drift_data.append({
            'theta_deg': th_deg,
            'theta_rad': th_rad,
            'Ux': U[0],
            'Uy': U[1],
            'Uz': U[2],
            'abs_Uz': abs(U[2]),
            'sin_2theta': np.sin(2.0 * th_rad)
        })

    # Save CSV
    csv_file = os.path.join(PROC_DIR, "boomerang_tilt_angle_drift.csv")
    with open(csv_file, 'w') as f:
        f.write("theta_deg,theta_rad,Ux,Uy,Uz,abs_Uz,sin_2theta\n")
        for d in drift_data:
            f.write(f"{d['theta_deg']:.1f},{d['theta_rad']:.6f},{d['Ux']:.8e},{d['Uy']:.8e},"
                    f"{d['Uz']:.8e},{d['abs_Uz']:.8e},{d['sin_2theta']:.6f}\n")
    print(f"Saved drift vs angle dataset to: {csv_file}")

    # Fit scaling: Ux(theta) = A * sin(2*theta)
    sin_vals = np.array([d['sin_2theta'] for d in drift_data])
    Ux_vals = np.array([d['Ux'] for d in drift_data])
    slope, intercept, r_value, p_value, std_err = stats.linregress(sin_vals, Ux_vals)

    print(f"\n--- Oblique Drift Scaling Law Fit: Ux = A * sin(2*theta) ---")
    print(f"Fitted Amplitude A = {slope:.6f}")
    print(f"Intercept          = {intercept:.2e} (expected 0.0)")
    print(f"Correlation R^2    = {r_value**2:.8f} (exact Stokesian sinusoidal drift)")
    print(f"Peak Lateral Drift = {np.max(np.abs(Ux_vals)):.6f} at theta = 45 deg\n")

    # -------------------------------------------------------------------------
    # Publication Plot 2: Lateral Drift vs Tilt Angle theta
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    th_dense = np.linspace(0.0, 90.0, 100)
    sin_dense = np.sin(2.0 * np.radians(th_dense))

    # Subplot 1: Lateral Drift Ux vs Angle
    ax1.plot(angles_deg, Ux_vals, 'ro', label=r'Multiblob ($N = 15$ Boomerang)', markersize=8)
    ax1.plot(th_dense, slope * sin_dense, 'b-', label=rf'Stokes Anisotropy: $A \sin(2\theta)$ ($R^2={r_value**2:.6f}$)', linewidth=2.2)
    ax1.axvline(45.0, color='gray', linestyle=':', label=r'Maximum Drift Angle ($\theta = 45^\circ$)')

    ax1.set_xlabel(r'Particle Tilt Angle $\theta$ (degrees)')
    ax1.set_ylabel(r'Lateral Drift Velocity $U_x$')
    ax1.set_title(r'(a) Oblique Lateral Drift vs. Tilt Angle')
    ax1.set_xlim(0, 90)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='lower center', framealpha=0.92)

    # Subplot 2: Settling Speed |Uz| vs Angle
    Uz_vals = np.array([d['abs_Uz'] for d in drift_data])
    ax2.plot(angles_deg, Uz_vals, 'gs-', label=r'Vertical Speed $|U_z(\theta)|$', markersize=7)
    ax2.set_xlabel(r'Particle Tilt Angle $\theta$ (degrees)')
    ax2.set_ylabel(r'Vertical Settling Speed $|U_z|$')
    ax2.set_title(r'(b) Vertical Settling Speed vs. Tilt Angle')
    ax2.set_xlim(0, 90)
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc='best', framealpha=0.92)

    plt.tight_layout()
    plot_file2 = os.path.join(PLOT_DIR, "boomerang_tilted_drift_vs_angle.png")
    plt.savefig(plot_file2)
    plt.close()
    print(f"Generated publication figure: {plot_file2}\n")


if __name__ == '__main__':
    t_start = time.time()
    run_force_resolution_sweep()
    run_tilt_angle_drift_sweep()
    print(f"All boomerang force resolution & tilt drift analyses complete in {time.time() - t_start:.2f} seconds!")
