"""
MTP_Experiments/scripts/phase2_boomerang_trajectory_sim.py

Dynamic Time-Stepping Sedimentation Simulation of a Boomerang Particle in Stokes Flow.
Demonstrates:
  1. Steady Edge-On Descent (Apex Down, aligned with gravity).
  2. Pitching Reorientation from Flat Pose toward Stable Terminal Pose.
  3. Oblique Lateral Drift at 45° tilt.
  4. Chiral / Dihedral Boomerang (15° out-of-plane twist) producing 3D Helical Spiral Descent.
  5. Exports trajectory JSON for the interactive 3D WebGL viewer.
  6. Generates 300 DPI publication trajectory comparison figures.
"""

import sys
import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
from body.body import Body
from quaternion_integrator.quaternion import Quaternion

# Directory configuration
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'boomerang')
RAW_DIR = os.path.join(OUT_DIR, 'raw_data')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')

for d in [RAW_DIR, PROC_DIR, PLOT_DIR]:
    os.makedirs(d, exist_ok=True)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})


def simulate_boomerang_trajectory(case_name='flat', initial_q=None, dihedral_deg=0.0,
                                  N_blobs=15, dt=0.5, n_steps=60, Fz=1.0, eta=1.0):
    """
    Simulates time-stepping 6-DOF dynamic trajectory of a boomerang particle under gravity.
    """
    r_conf, a_blob, meta = p2c.generate_boomerang_mesh(
        N_blobs=N_blobs, arm_length=2.1, angle_deg=90.0,
        dihedral_deg=dihedral_deg, center_at='centroid'
    )

    pos = np.array([0.0, 0.0, 0.0])
    q = Quaternion([1.0, 0.0, 0.0, 0.0]) if initial_q is None else initial_q
    force_vec = np.array([0.0, 0.0, -Fz])

    traj = []
    current_time = 0.0

    print(f"Simulating Boomerang: {case_name:18s} (dihedral={dihedral_deg:4.1f}°, N={N_blobs})...")

    for step in range(n_steps + 1):
        # Solve instantaneous Stokes velocities in laboratory frame
        U, Omega, _ = p2c.solve_sedimentation(r_conf, a_blob, eta, force_vec, orientation=q)

        # Compute blob coordinates in laboratory frame for rendering
        rot_mat = q.rotation_matrix()
        r_lab = pos + np.dot(r_conf, rot_mat.T)

        traj.append({
            'step': step,
            'time': float(current_time),
            'x': float(pos[0]),
            'y': float(pos[1]),
            'z': float(pos[2]),
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
            'blobs_lab': r_lab.tolist()
        })

        if step == n_steps:
            break

        # Translational update (Forward Euler)
        pos = pos + U * dt

        # Rotational update via quaternion kinematics: dq = exp(Omega * dt / 2)
        dq = Quaternion.from_rotation(Omega * dt)
        q = dq * q
        # Normalize quaternion
        q_norm = np.linalg.norm(q.entries)
        if q_norm > 0:
            q.entries = q.entries / q_norm
            q.s = np.array(q.entries[0])
            q.p = np.array(q.entries[1:4])

        current_time += dt

    return traj, r_conf, a_blob, meta


def export_multi_resolution_meshes():
    """
    Exports reference blob configurations for resolutions N=7, 15, 29
    for planar and chiral boomerangs.
    """
    mesh_data = {}
    cases = [
        ('planar', 0.0),
        ('chiral', 15.0)
    ]
    resolutions = [7, 15, 29]

    for cname, dih in cases:
        mesh_data[cname] = {}
        for N in resolutions:
            r_conf, a_blob, meta = p2c.generate_boomerang_mesh(
                N_blobs=N, arm_length=2.1, angle_deg=90.0,
                dihedral_deg=dih, center_at='centroid'
            )
            mesh_data[cname][f"res_{N}"] = {
                'N_blobs': len(r_conf),
                'blob_radius': float(a_blob),
                'r_conf': r_conf.tolist(),
                'centroid': meta['centroid'],
                'apex_coord': meta['apex_coord']
            }
    return mesh_data


def plot_trajectory_comparison(traj_edge, traj_flat, traj_tilted, traj_chiral):
    """
    Generates a 300 DPI publication figure comparing boomerang trajectories.
    """
    fig = plt.figure(figsize=(16, 10))

    # Subplot 1: 3D Trajectory paths
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')

    def get_xyz(traj):
        return [d['x'] for d in traj], [d['y'] for d in traj], [d['z'] for d in traj]

    xe, ye, ze = get_xyz(traj_edge)
    xf, yf, zf = get_xyz(traj_flat)
    xt, yt, zt = get_xyz(traj_tilted)
    xc, yc, zc = get_xyz(traj_chiral)

    ax1.plot(xe, ye, ze, 'b-', lw=2.2, label='Apex Down (Edge-On Gliding)')
    ax1.plot(xf, yf, zf, 'g--', lw=2.2, label='Flat Pose (Pitching Reorientation)')
    ax1.plot(xt, yt, zt, 'm-.', lw=2.2, label='Tilted 45° (Oblique Drift)')
    ax1.plot(xc, yc, zc, 'r-', lw=2.5, label='Chiral Boomerang (3D Helical Spiral)')

    ax1.scatter([0], [0], [0], color='black', s=50, label='Start (0,0,0)')
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position')
    ax1.set_zlabel('Z (Settling)')
    ax1.set_title('3D Sedimentation Trajectories', fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.legend(loc='lower left', fontsize=8)

    # Subplot 2: Horizontal Drift (X vs Y)
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(xe, ye, 'bo-', lw=1.8, ms=4, label='Apex Down')
    ax2.plot(xf, yf, 'gs-', lw=1.8, ms=4, label='Flat Pose')
    ax2.plot(xt, yt, 'md-', lw=1.8, ms=4, label='Tilted 45° (Lateral Drift)')
    ax2.plot(xc, yc, 'r^-', lw=2.0, ms=4, label='Chiral Spiral (Orbit)')
    ax2.set_xlabel('Horizontal X')
    ax2.set_ylabel('Horizontal Y')
    ax2.set_title('Horizontal Projection (X-Y Plane)', fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()
    ax2.axis('equal')

    # Subplot 3: Settling Speed |Uz| vs Time
    t_vals = [d['time'] for d in traj_edge]
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_edge], 'b-', lw=2, label='Apex Down')
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_flat], 'g--', lw=2, label='Flat Pose')
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_tilted], 'm-.', lw=2, label='Tilted 45°')
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_chiral], 'r-', lw=2, label='Chiral Spiral')
    ax3.set_xlabel('Dimensionless Time $t / \\tau_c$')
    ax3.set_ylabel('Settling Speed $|U_z|$')
    ax3.set_title('Settling Speed Evolution', fontweight='bold')
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend()

    # Subplot 4: Induced Rotation Rate ||Omega|| vs Time
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_edge], 'b-', lw=2, label='Apex Down: $\\|\\mathbf{\\Omega}\\| \\approx 0$')
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_flat], 'g--', lw=2, label='Flat: In-Plane Pitching')
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_tilted], 'm-.', lw=2, label='Tilted 45°: Reorientation')
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_chiral], 'r-', lw=2, label='Chiral: Continuous Autorotation')
    ax4.set_xlabel('Dimensionless Time $t / \\tau_c$')
    ax4.set_ylabel(r'Angular Velocity Magnitude $\|\mathbf{\Omega}\|$')
    ax4.set_title('Induced Gravitational Rotation', fontweight='bold')
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend()

    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "boomerang_dynamic_simulation_trajectories.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"Saved 3D trajectory comparison figure: {plot_path}")


def main():
    print("\n" + "="*70)
    print("BOOMERANG PARTICLE DYNAMIC SEDIMENTATION TRAJECTORY SIMULATION")
    print("="*70 + "\n")

    n_steps = 60
    dt = 0.5

    # Case 1: Apex Down (bisector aligned with z: rotate -45° about z, then tilt into vertical)
    q_edge = Quaternion.from_rotation(np.array([0.0, np.pi / 2.0, 0.0])) * Quaternion.from_rotation(np.array([0.0, 0.0, -np.pi / 4.0]))
    traj_edge, _, _, _ = simulate_boomerang_trajectory('Apex Down', initial_q=q_edge, dt=dt, n_steps=n_steps)

    # Case 2: Flat Pose (horizontal xy plane)
    traj_flat, _, _, _ = simulate_boomerang_trajectory('Flat Pose', initial_q=None, dt=dt, n_steps=n_steps)

    # Case 3: Tilted 45° (in-plane tilt toward gravity)
    q_tilt = Quaternion.from_rotation(np.array([0.0, np.radians(45.0), 0.0]))
    traj_tilted, _, _, _ = simulate_boomerang_trajectory('Tilted 45°', initial_q=q_tilt, dt=dt, n_steps=n_steps)

    # Case 4: Chiral / Dihedral Boomerang (15° twist)
    traj_chiral, _, _, _ = simulate_boomerang_trajectory('Chiral Spiral', initial_q=None, dihedral_deg=15.0, dt=dt, n_steps=n_steps)

    # Export multi-resolution meshes
    print("\nExporting multi-resolution mesh data (N=7, 15, 29 blobs)...")
    mesh_data = export_multi_resolution_meshes()

    # Save comprehensive simulation JSON
    json_path = os.path.join(PROC_DIR, "boomerang_simulation_data.json")
    payload = {
        'metadata': {
            'N_blobs': 15,
            'arm_length': 2.1,
            'angle_deg': 90.0,
            'dt': dt,
            'n_steps': n_steps,
            'total_time': dt * n_steps,
            'Fz': 1.0,
            'eta': 1.0
        },
        'meshes': mesh_data,
        'trajectory_edge': traj_edge,
        'trajectory_flat': traj_flat,
        'trajectory_tilted': traj_tilted,
        'trajectory_chiral': traj_chiral
    }

    with open(json_path, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"Saved comprehensive simulation JSON: {json_path}")

    # Plot 300 DPI publication trajectory figure
    print("\nGenerating publication figures...")
    plot_trajectory_comparison(traj_edge, traj_flat, traj_tilted, traj_chiral)

    print("\n" + "="*70)
    print("Boomerang Dynamic Simulation Complete!")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
