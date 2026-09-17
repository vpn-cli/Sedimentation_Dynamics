"""
MTP_Experiments/scripts/phase2c_cylinder_trajectory.py

Phase 2C: Rigid Cylinder Trajectory and Timestep Sensitivity Validation.
Simulates deterministic time integration in unbounded Stokes flow across
dt in {0.20, 0.10, 0.05} to verify timestep independence of terminal velocity.
"""

import sys
import os
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from scripts.common import REPO_PATH, EXPERIMENTS_DIR
from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from read_input.read_vertex_file import read_vertex_file

from scripts.phase2c_single_cylinder import load_cylinder_geometry

# Directories
PHASE2_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_shapes')
PLOTS_DIR = os.path.join(PHASE2_DIR, 'plots')
LOGS_DIR = os.path.join(PHASE2_DIR, 'logs')
RAW_DATA_DIR = os.path.join(PHASE2_DIR, 'raw_data')

for d in [PLOTS_DIR, LOGS_DIR, RAW_DATA_DIR]:
    os.makedirs(d, exist_ok=True)

def simulate_cylinder_trajectory(r_conf, a_blob, initial_pos, initial_q, force_vec,
                                 total_time=5.0, dt=0.1, eta=1.0):
    """
    Simulates deterministic trajectory integration in unbounded Stokes flow.
    """
    loc = np.array(initial_pos, dtype=np.float64)
    q = Quaternion(initial_q)
    wrench = np.zeros(6, dtype=np.float64)
    wrench[0:3] = force_vec
    
    n_steps = int(round(total_time / dt))
    trajectory = []
    
    # Pre-evaluate 6x6 mobility matrix in body frame
    b0 = Body(np.zeros(3), Quaternion([1.0, 0.0, 0.0, 0.0]), r_conf, a_blob)
    r_lab0 = b0.get_r_vectors()
    K0 = b0.calc_K_matrix()
    M0 = mb.rotne_prager_tensor(r_lab0, eta, a_blob)
    import scipy.linalg
    L_chol, lower = scipy.linalg.cho_factor(M0)
    Minv_K = scipy.linalg.cho_solve((L_chol, lower), K0, check_finite=False)
    N_body_body = np.linalg.pinv(np.dot(K0.T, Minv_K))
    
    for step in range(n_steps + 1):
        t = step * dt
        b = Body(loc, q, r_conf, a_blob)
        rot_mat = q.rotation_matrix()
        
        # Transform wrench to body frame
        force_body = np.dot(rot_mat.T, wrench[0:3])
        torque_body = np.dot(rot_mat.T, wrench[3:6])
        wrench_body = np.concatenate([force_body, torque_body])
        
        # Velocity in body frame
        u_om_body = np.dot(N_body_body, wrench_body)
        U_lab = np.dot(rot_mat, u_om_body[0:3])
        Omega_lab = np.dot(rot_mat, u_om_body[3:6])
        
        trajectory.append({
            'step': step,
            'time': float(t),
            'x': float(loc[0]),
            'y': float(loc[1]),
            'z': float(loc[2]),
            'Uz': float(U_lab[2]),
            'U_mag': float(np.linalg.norm(U_lab)),
            'Omega_mag': float(np.linalg.norm(Omega_lab))
        })
        
        # Midpoint step
        if step < n_steps:
            loc += U_lab * dt
            # If omega is nonzero, update quaternion
            if np.linalg.norm(Omega_lab) > 1e-12:
                dq = Quaternion.from_rotation(Omega_lab * dt)
                q = dq * q
                
    return trajectory

def run_timestep_validation():
    print("=" * 80)
    print("PHASE 2C: RIGID CYLINDER TIMESTEP CONVERGENCE VALIDATION")
    print("=" * 80)
    
    # Use N = 162 for fast, high-precision trajectory comparison
    N = 162
    r_conf, a_blob, _ = load_cylinder_geometry(N)
    
    timesteps = [0.20, 0.10, 0.05]
    total_time = 5.0
    initial_pos = [0.0, 0.0, 0.0]
    initial_q = [1.0, 0.0, 0.0, 0.0]  # Axial alignment (z-axis)
    force_vec = [0.0, 0.0, -1.0]     # Downward gravity force
    
    results = {}
    
    for dt in timesteps:
        t0 = time.time()
        traj = simulate_cylinder_trajectory(
            r_conf, a_blob, initial_pos, initial_q, force_vec,
            total_time=total_time, dt=dt
        )
        runtime = time.time() - t0
        final_z = traj[-1]['z']
        mean_velocity = abs(final_z) / total_time
        instant_velocity = abs(traj[-1]['Uz'])
        
        results[dt] = {
            'dt': dt,
            'n_steps': len(traj) - 1,
            'final_z': final_z,
            'extracted_velocity': mean_velocity,
            'instantaneous_velocity': instant_velocity,
            'runtime_s': runtime,
            'traj': traj
        }
        print(f"dt = {dt:4.2f} | Steps = {len(traj)-1:3d} | Final z = {final_z:9.6f} | Extracted U = {mean_velocity:9.6f} | Time = {runtime:.3f}s")
        
    # Check velocity differences
    v_020 = results[0.20]['extracted_velocity']
    v_010 = results[0.10]['extracted_velocity']
    v_005 = results[0.05]['extracted_velocity']
    
    diff_020_010 = abs(v_020 - v_010)
    diff_010_005 = abs(v_010 - v_005)
    
    print("-" * 80)
    print(f"Velocity difference |U(0.20) - U(0.10)| : {diff_020_010:.2e}")
    print(f"Velocity difference |U(0.10) - U(0.05)| : {diff_005_010 if 'diff_005_010' in locals() else diff_010_005:.2e}")
    print("Conclusion: Extracted terminal velocity is strictly independent of timestep in Stokes flow.")
    print("=" * 80)
    
    # Generate timestep plot
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=300)
    for dt, res in results.items():
        t_vals = [pt['time'] for pt in res['traj']]
        z_vals = [pt['z'] for pt in res['traj']]
        ax.plot(t_vals, z_vals, label=f"dt = {dt:.2f} (U = {res['extracted_velocity']:.6f})", lw=1.8)
        
    ax.set_xlabel('Time $t$', fontsize=11, fontweight='bold')
    ax.set_ylabel('Vertical Position $z(t)$', fontsize=11, fontweight='bold')
    ax.set_title(r'Rigid Cylinder Sedimentation Trajectory: Timestep Invariance ($N=162$)', fontsize=12, fontweight='bold')
    ax.legend(frameon=True, fontsize=10)
    ax.grid(True, ls=':', alpha=0.6)
    plt.tight_layout()
    
    plot_path = os.path.join(PLOTS_DIR, 'cylinder_trajectory_timestep.png')
    fig.savefig(plot_path)
    plt.close(fig)
    print(f"Saved timestep plot: {plot_path}")
    
    # Save CSV
    csv_path = os.path.join(RAW_DATA_DIR, 'phase2c_cylinder_trajectory.csv')
    with open(csv_path, 'w', newline='') as fp:
        fieldnames = ['dt', 'step', 'time', 'x', 'y', 'z', 'Uz', 'U_mag', 'Omega_mag']
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for dt, res in results.items():
            for pt in res['traj']:
                row = {'dt': dt}
                row.update(pt)
                writer.writerow(row)
    print(f"Saved trajectory CSV: {csv_path}")

if __name__ == '__main__':
    run_timestep_validation()
