"""
MTP_Experiments/scripts/phase2a_disc_trajectory.py

Dynamic time-stepping trajectory simulations of single rigid disc sedimentation.
1. Unbounded Stokes flow: demonstrates constant terminal velocity, horizontal drift, and orientation conservation.
2. Near-wall Stokes flow (Swan-Brady wall mobility): demonstrates wall drag enhancement and wall-induced reorientation torque.
3. Timestep sensitivity study: validates numerical convergence across timestep sizes.
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

PHASE2_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_shapes')
STRUCTURES_DIR = os.path.join(PHASE2_DIR, 'structures')
RAW_DATA_DIR = os.path.join(PHASE2_DIR, 'raw_data')
PROCESSED_DATA_DIR = os.path.join(PHASE2_DIR, 'processed_data')
PLOTS_DIR = os.path.join(PHASE2_DIR, 'plots')
REPORTS_DIR = os.path.join(PHASE2_DIR, 'reports')
LOGS_DIR = os.path.join(PHASE2_DIR, 'logs')

def solve_instantaneous_velocity(body_obj, r_conf, a_blob, wrench, eta=1.0, domain='unbounded'):
    """
    Computes instantaneous rigid-body velocity U and Omega for current body state.
    Supports domain='unbounded' (RPY) or domain='single_wall' (Swan & Brady).
    """
    r_lab = body_obj.get_r_vectors()
    K = body_obj.calc_K_matrix()
    
    if domain == 'unbounded':
        M = mb.rotne_prager_tensor(r_lab, eta, a_blob)
    elif domain == 'single_wall':
        M = mb.single_wall_fluid_mobility(r_lab, eta, a_blob)
    else:
        raise ValueError(f"Unknown domain: {domain}")
        
    L, lower = scipy.linalg.cho_factor(M)
    Kt_Minv_K = np.dot(K.T, scipy.linalg.cho_solve((L, lower), K, check_finite=False))
    N_body = np.linalg.pinv(Kt_Minv_K)
    
    U_Omega = np.dot(N_body, wrench)
    U = U_Omega[0:3]
    Omega = U_Omega[3:6]
    return U, Omega

def simulate_trajectory(r_conf, a_blob, initial_pos, initial_q, force_vec,
                        dt=0.1, n_steps=200, eta=1.0, domain='unbounded', scheme='midpoint'):
    """
    Simulates the deterministic trajectory of a single rigid disc over time.
    Uses 2nd-order Runge-Kutta / midpoint or forward Euler integration.
    """
    loc = np.array(initial_pos, dtype=np.float64)
    q = Quaternion(initial_q)
    wrench = np.zeros(6, dtype=np.float64)
    wrench[0:3] = force_vec
    
    trajectory = []
    
    for step in range(n_steps + 1):
        t = step * dt
        b = Body(loc, q, r_conf, a_blob)
        U, Omega = solve_instantaneous_velocity(b, r_conf, a_blob, wrench, eta=eta, domain=domain)
        
        # Calculate current tilt angle relative to vertical z
        rot_mat = q.rotation_matrix()
        normal_lab = np.dot(rot_mat, np.array([0.0, 0.0, 1.0]))
        tilt_angle_deg = np.degrees(np.arccos(np.clip(normal_lab[2], -1.0, 1.0)))
        
        record = {
            'step': step,
            'time': float(t),
            'x': float(loc[0]),
            'y': float(loc[1]),
            'z': float(loc[2]),
            'Ux': float(U[0]),
            'Uy': float(U[1]),
            'Uz': float(U[2]),
            'U_mag': float(np.linalg.norm(U)),
            'Omega_x': float(Omega[0]),
            'Omega_y': float(Omega[1]),
            'Omega_z': float(Omega[2]),
            'Omega_mag': float(np.linalg.norm(Omega)),
            'q0': float(q.entries[0]),
            'q1': float(q.entries[1]),
            'q2': float(q.entries[2]),
            'q3': float(q.entries[3]),
            'normal_x': float(normal_lab[0]),
            'normal_y': float(normal_lab[1]),
            'normal_z': float(normal_lab[2]),
            'tilt_angle_deg': float(tilt_angle_deg)
        }
        trajectory.append(record)
        
        if step == n_steps:
            break
            
        # Time integration
        if scheme == 'euler':
            loc = loc + U * dt
            q_dt = Quaternion.from_rotation(Omega * dt)
            q = q_dt * q
            q = Quaternion(q.entries / np.linalg.norm(q.entries))
        elif scheme == 'midpoint':
            # Half-step prediction
            loc_mid = loc + 0.5 * U * dt
            q_dt_half = Quaternion.from_rotation(0.5 * Omega * dt)
            q_mid = q_dt_half * q
            q_mid = Quaternion(q_mid.entries / np.linalg.norm(q_mid.entries))
            
            b_mid = Body(loc_mid, q_mid, r_conf, a_blob)
            U_mid, Omega_mid = solve_instantaneous_velocity(b_mid, r_conf, a_blob, wrench, eta=eta, domain=domain)
            
            loc = loc + U_mid * dt
            q_dt_full = Quaternion.from_rotation(Omega_mid * dt)
            q = q_dt_full * q
            q = Quaternion(q.entries / np.linalg.norm(q.entries))
            
        # Termination if particle hits wall (z <= a_blob)
        if domain == 'single_wall' and (loc[2] - np.max(np.abs(r_conf[:, 0])) * np.sin(np.radians(tilt_angle_deg)) <= a_blob):
            print(f"  Disc reached near-wall cutoff at step {step}, t={t:.2f} s")
            break
            
    return trajectory

def run_trajectory_studies():
    log_lines = []
    def log(msg):
        print(msg)
        log_lines.append(msg)
        
    log("=" * 80)
    log("PHASE 2A: DYNAMIC TRAJECTORY & WALL EFFECT SIMULATIONS")
    log("=" * 80)
    
    # Use standard medium disc (N = 95 blobs) for efficient, high-accuracy trajectory dynamics
    disc_file = os.path.join(STRUCTURES_DIR, "disc_N_95_R_1.0_rings_5.vertex")
    r_conf = read_vertex_file(disc_file)
    with open(disc_file, 'r') as fp:
        for line in fp:
            if not line.startswith('#'):
                a_blob = float(line.split()[1]) / 2.0
                break
                
    R_geom = 1.0
    eta = 1.0
    F_mag = 1.0
    force_vec = np.array([0.0, 0.0, -F_mag])
    
    # -------------------------------------------------------------
    # 1. Unbounded Trajectories: Face-on, Edge-on, and Tilted (45 deg)
    # -------------------------------------------------------------
    log("\n>>> Simulating Unbounded Stokes Trajectories (N=95, dt=0.1, t_max=10.0 s) <<<")
    
    # Face-on
    q_face = [1.0, 0.0, 0.0, 0.0]
    traj_face = simulate_trajectory(r_conf, a_blob, [0, 0, 10], q_face, force_vec, dt=0.1, n_steps=100, domain='unbounded')
    log(f"  Face-on final z = {traj_face[-1]['z']:.4f} | Terminal Uz = {traj_face[-1]['Uz']:.6f} | Omega = {traj_face[-1]['Omega_mag']:.2e}")
    
    # Edge-on
    q_edge = [np.cos(np.pi/4), 0.0, np.sin(np.pi/4), 0.0]
    traj_edge = simulate_trajectory(r_conf, a_blob, [0, 0, 10], q_edge, force_vec, dt=0.1, n_steps=100, domain='unbounded')
    log(f"  Edge-on final z = {traj_edge[-1]['z']:.4f} | Terminal Uz = {traj_edge[-1]['Uz']:.6f} | Omega = {traj_edge[-1]['Omega_mag']:.2e}")
    
    # Tilted 45 deg
    theta_tilt = np.radians(45.0)
    q_tilt = [np.cos(theta_tilt/2), 0.0, np.sin(theta_tilt/2), 0.0]
    traj_tilt = simulate_trajectory(r_conf, a_blob, [0, 0, 10], q_tilt, force_vec, dt=0.1, n_steps=100, domain='unbounded')
    log(f"  Tilted final x = {traj_tilt[-1]['x']:.4f} (drift) | final z = {traj_tilt[-1]['z']:.4f} | final tilt = {traj_tilt[-1]['tilt_angle_deg']:.2f}°")
    
    # -------------------------------------------------------------
    # 2. Wall Interaction Trajectory: Tilted Disc Approaching Wall
    # -------------------------------------------------------------
    log("\n>>> Simulating Near-Wall Sedimentation with Swan-Brady Mobility <<<")
    # Start closer to wall to highlight hydrodynamic wall torque: H0 = 4.0
    traj_wall = simulate_trajectory(r_conf, a_blob, [0, 0, 3.5], q_tilt, force_vec, dt=0.1, n_steps=120, domain='single_wall')
    log(f"  Wall Trajectory: initial tilt = {traj_wall[0]['tilt_angle_deg']:.2f}° -> final tilt = {traj_wall[-1]['tilt_angle_deg']:.2f}°")
    log(f"  Wall Trajectory: initial Uz = {traj_wall[0]['Uz']:.5f} -> final Uz = {traj_wall[-1]['Uz']:.5f} (wall deceleration)")
    
    # -------------------------------------------------------------
    # 3. Timestep Sensitivity Study (dt in {0.2, 0.1, 0.05})
    # -------------------------------------------------------------
    log("\n>>> Running Timestep Sensitivity Study (dt = 0.2, 0.1, 0.05) <<<")
    timesteps = [0.2, 0.1, 0.05]
    dt_trajectories = {}
    for dt_val in timesteps:
        n_st = int(round(6.0 / dt_val))
        t_data = simulate_trajectory(r_conf, a_blob, [0, 0, 10], q_tilt, force_vec, dt=dt_val, n_steps=n_st, domain='unbounded')
        dt_trajectories[dt_val] = t_data
        log(f"  dt = {dt_val:4.2f} s | Steps = {n_st:3d} | Final x = {t_data[-1]['x']:.6f} | Final z = {t_data[-1]['z']:.6f} | Final Uz = {t_data[-1]['Uz']:.6f}")

    # Save CSVs
    csv_unbounded = os.path.join(RAW_DATA_DIR, 'phase2a_trajectory_unbounded.csv')
    with open(csv_unbounded, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(traj_tilt[0].keys()))
        writer.writeheader()
        writer.writerows(traj_tilt)
    log(f"\n[Unbounded Trajectory Saved] -> {csv_unbounded}")
    
    csv_wall = os.path.join(RAW_DATA_DIR, 'phase2a_trajectory_wall.csv')
    with open(csv_wall, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(traj_wall[0].keys()))
        writer.writeheader()
        writer.writerows(traj_wall)
    log(f"[Wall Trajectory Saved]      -> {csv_wall}")
    
    # Generate Trajectory Plots
    generate_trajectory_plots(traj_face, traj_edge, traj_tilt, traj_wall, dt_trajectories, PLOTS_DIR)
    log(f"[Trajectory Plots Saved]     -> {PLOTS_DIR}")
    
    # Save Log
    log_path = os.path.join(LOGS_DIR, 'phase2a_disc_trajectory.log')
    with open(log_path, 'w', encoding='utf-8') as fp:
        fp.write("\n".join(log_lines))
    log(f"[Log Saved]                  -> {log_path}")
    log("=" * 80)

def generate_trajectory_plots(traj_face, traj_edge, traj_tilt, traj_wall, dt_trajs, output_dir):
    """Plots trajectory comparisons, horizontal drift, wall reorientation, and timestep checks."""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # 1. z(t) Comparison: Face-on vs Edge-on vs Tilted (Unbounded)
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    t_f = [r['time'] for r in traj_face]
    z_f = [r['z'] for r in traj_face]
    t_e = [r['time'] for r in traj_edge]
    z_e = [r['z'] for r in traj_edge]
    t_t = [r['time'] for r in traj_tilt]
    z_t = [r['z'] for r in traj_tilt]
    
    ax.plot(t_f, z_f, '-', color='#1f77b4', lw=2.2, label='Face-on (Slowest Settling)')
    ax.plot(t_e, z_e, '-', color='#d62728', lw=2.2, label='Edge-on (Fastest Settling)')
    ax.plot(t_t, z_t, '--', color='#2ca02c', lw=2.2, label='Tilted 45° (Intermediate Settling)')
    ax.set_xlabel('Time $t$ (seconds)', fontsize=12)
    ax.set_ylabel('Vertical Position $z(t)$', fontsize=12)
    ax.set_title('Settling Trajectories Across Orientations (Unbounded Stokes Flow)', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=11)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'trajectory_z_vs_time.png'))
    plt.close(fig)
    
    # 2. Horizontal Drift (x-z Plane Trajectory for Tilted Disc)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    x_t = [r['x'] for r in traj_tilt]
    ax.plot(x_t, z_t, 'o-', color='#9467bd', lw=2, ms=4, markevery=5, label='Tilted Disc Path (45°)')
    ax.plot([0, 0], [z_t[0], z_t[-1]], 'k:', lw=1.5, label='Gravity Direction')
    ax.set_xlabel('Lateral Position $x(t)$ (Drift)', fontsize=12)
    ax.set_ylabel('Vertical Position $z(t)$', fontsize=12)
    ax.set_title('Lateral Drift Trajectory of a Tilted Disc', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=11)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'trajectory_lateral_drift_xz.png'))
    plt.close(fig)
    
    # 3. Wall Effect: Deceleration and Wall-Induced Reorientation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), dpi=300)
    t_w = [r['time'] for r in traj_wall]
    z_w = [r['z'] for r in traj_wall]
    Uz_w = [abs(r['Uz']) for r in traj_wall]
    tilt_w = [r['tilt_angle_deg'] for r in traj_wall]
    
    ax1.plot(z_w, Uz_w, 'o-', color='#d62728', lw=2, ms=3)
    ax1.set_xlabel('Height Above Wall $z$', fontsize=12)
    ax1.set_ylabel('Settling Speed $|U_z|$', fontsize=12)
    ax1.set_title('Wall Lubrication: Settling Velocity Deceleration', fontsize=11, fontweight='bold')
    ax1.grid(True, ls=':', alpha=0.6)
    
    ax2.plot(t_w, tilt_w, 's-', color='#1f77b4', lw=2, ms=3)
    ax2.set_xlabel('Time $t$ (seconds)', fontsize=12)
    ax2.set_ylabel('Disc Tilt Angle $\\theta$ (degrees)', fontsize=12)
    ax2.set_title('Wall-Induced Hydrodynamic Reorientation', fontsize=11, fontweight='bold')
    ax2.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'trajectory_wall_effect_dynamics.png'))
    plt.close(fig)
    
    # 4. Timestep Convergence Plot
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    for idx, (dt_val, tr) in enumerate(dt_trajs.items()):
        t_arr = [r['time'] for r in tr]
        z_arr = [r['z'] for r in tr]
        ax.plot(t_arr, z_arr, label=f'$\\Delta t = {dt_val}$ s', lw=2.0 - idx*0.4, ls=['-', '--', ':'][idx], color=colors[idx])
    ax.set_xlabel('Time $t$ (seconds)', fontsize=12)
    ax.set_ylabel('Vertical Position $z(t)$', fontsize=12)
    ax.set_title('Timestep Sensitivity & Convergence ($z(t)$ across $\\Delta t$)', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=11)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'trajectory_timestep_convergence.png'))
    plt.close(fig)

if __name__ == '__main__':
    run_trajectory_studies()
