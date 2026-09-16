"""
MTP_Experiments/scripts/phase2_robotic_arm_trajectory_sim.py

Dynamic Time-Stepping Sedimentation Simulation of a Robotic Arm in Stokes Flow.
Demonstrates:
  1. Steady non-tumbling descent of a Straight Symmetric Arm (Omega = 0).
  2. Pitching reorientation of a Planar Bent Arm (theta_bend = 45 deg) toward equilibrium.
  3. Continuous 3D autorotation and helical spiraling descent of a Chiral Twisted Arm (theta_bend = 45 deg, phi_twist = 45 deg).
  4. Exports trajectory JSON for the interactive 3D WebGL viewer.
  5. Generates 300 DPI publication trajectory comparison figures.
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
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'robotic_arm')
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


def simulate_arm_trajectory(config='straight', bend_deg=0.0, twist_deg=0.0,
                            link_resolution=12, dt=1.0, n_steps=60, Fz=1.0, eta=1.0):
    """
    Simulates time-stepping 6-DOF dynamic trajectory of a rigid robotic arm under gravity.
    """
    r_conf, a_blob, meta = p2c.assemble_robotic_arm_rigid(
        N_links=7, link_resolution=link_resolution, link_spacing=2.5,
        config=config, bend_angle_deg=bend_deg, twist_angle_deg=twist_deg,
        center_at='centroid'
    )
    link_centers_ref = np.array(meta['link_centers'])

    pos = np.array([0.0, 0.0, 0.0])
    q = Quaternion([1.0, 0.0, 0.0, 0.0])
    force_vec = np.array([0.0, 0.0, -Fz])

    traj = []
    current_time = 0.0

    print(f"Simulating Arm: {config:8s} (bend={bend_deg:4.1f}°, twist={twist_deg:4.1f}°, N_blobs={len(r_conf)})...")

    for step in range(n_steps + 1):
        # Solve instantaneous Stokes velocities in laboratory frame
        U, Omega, _ = p2c.solve_sedimentation(r_conf, a_blob, eta, force_vec, orientation=q)

        # Compute current link centers in laboratory frame
        rot_mat = q.rotation_matrix()
        link_centers_lab = pos + np.dot(link_centers_ref, rot_mat.T)

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
            'link_centers': link_centers_lab.tolist()
        })

        if step == n_steps:
            break

        # Translational update (Forward Euler)
        pos = pos + U * dt

        # Rotational update via quaternion kinematics: dq = exp(Omega * dt / 2)
        # Note: In laboratory frame, rotation about Omega is premultiplied
        dq = Quaternion.from_rotation(Omega * dt)
        q = dq * q
        # Normalize quaternion to prevent numerical drift
        q_norm = np.linalg.norm(q.entries)
        if q_norm > 0:
            q.entries = q.entries / q_norm
            q.s = np.array(q.entries[0])
            q.p = np.array(q.entries[1:4])

        current_time += dt

    return traj, r_conf, a_blob, meta


def export_multi_resolution_meshes():
    """
    Exports reference blob configurations and link centers for all 3 configurations
    across resolutions (1-blob N=7, 12-blob N=84, 42-blob N=294).
    """
    mesh_data = {}
    configs = [
        ('straight', 0.0, 0.0),
        ('bent', 45.0, 0.0),
        ('chiral', 45.0, 45.0)
    ]
    resolutions = [1, 12, 42]

    for cfg, bend, twist in configs:
        mesh_data[cfg] = {}
        for res in resolutions:
            r_conf, a_blob, meta = p2c.assemble_robotic_arm_rigid(
                N_links=7, link_resolution=res, link_spacing=2.5,
                config=cfg, bend_angle_deg=bend, twist_angle_deg=twist,
                center_at='centroid'
            )
            mesh_data[cfg][f"res_{res}"] = {
                'N_blobs': len(r_conf),
                'blob_radius': float(a_blob),
                'r_conf': r_conf.tolist(),
                'link_centers': meta['link_centers']
            }
    return mesh_data


def plot_trajectory_comparison(traj_straight, traj_bent, traj_chiral):
    """
    Generates a 300 DPI publication-grade comparison figure showing:
      1. 3D spatial settling trajectories (straight, bent, chiral spiral).
      2. 2D horizontal drift (X-Y paths).
      3. Settling speed |Uz|(t) and induced rotation rate ||Omega||(t).
    """
    fig = plt.figure(figsize=(16, 10))

    # Subplot 1: 3D Trajectory paths
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    
    xs_s = [d['x'] for d in traj_straight]
    ys_s = [d['y'] for d in traj_straight]
    zs_s = [d['z'] for d in traj_straight]

    xs_b = [d['x'] for d in traj_bent]
    ys_b = [d['y'] for d in traj_bent]
    zs_b = [d['z'] for d in traj_bent]

    xs_c = [d['x'] for d in traj_chiral]
    ys_c = [d['y'] for d in traj_chiral]
    zs_c = [d['z'] for d in traj_chiral]

    ax1.plot(xs_s, ys_s, zs_s, 'b-', lw=2.2, label='Straight Arm (Steady Non-Tumbling)')
    ax1.plot(xs_b, ys_b, zs_b, 'g--', lw=2.2, label='Bent Arm 45° (Planar Pitching)')
    ax1.plot(xs_c, ys_c, zs_c, 'r-', lw=2.5, label='Chiral Arm 45° (3D Helical Spiral)')

    # Scatter starting and ending points
    ax1.scatter([0], [0], [0], color='black', s=50, label='Start (0,0,0)')
    ax1.scatter([xs_s[-1]], [ys_s[-1]], [zs_s[-1]], color='blue', s=40)
    ax1.scatter([xs_b[-1]], [ys_b[-1]], [zs_b[-1]], color='green', s=40)
    ax1.scatter([xs_c[-1]], [ys_c[-1]], [zs_c[-1]], color='red', s=40)

    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position')
    ax1.set_zlabel('Z (Settling)')
    ax1.set_title('3D Sedimentation Trajectories', fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.legend(loc='lower left', fontsize=9)

    # Subplot 2: 2D Horizontal Drift (X vs Y)
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(xs_s, ys_s, 'bo-', lw=1.8, ms=4, label='Straight Arm')
    ax2.plot(xs_b, ys_b, 'gs-', lw=1.8, ms=4, label='Bent Arm 45°')
    ax2.plot(xs_c, ys_c, 'r^-', lw=2.0, ms=4, label='Chiral Arm 45° (Spiral Orbit)')
    ax2.set_xlabel('Horizontal X')
    ax2.set_ylabel('Horizontal Y')
    ax2.set_title('Horizontal Projection (Lateral Drift & Orbit)', fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()
    ax2.axis('equal')

    # Subplot 3: Settling Velocity |Uz| vs Time
    t_vals = [d['time'] for d in traj_straight]
    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_straight], 'b-', lw=2, label='Straight Arm')
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_bent], 'g--', lw=2, label='Bent Arm 45°')
    ax3.plot(t_vals, [abs(d['Uz']) for d in traj_chiral], 'r-', lw=2, label='Chiral Arm 45°')
    ax3.set_xlabel('Dimensionless Time $t / \\tau_c$')
    ax3.set_ylabel('Settling Speed $|U_z|$')
    ax3.set_title('Terminal Settling Speed Evolution', fontweight='bold')
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend()

    # Subplot 4: Induced Angular Velocity ||Omega|| vs Time
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_straight], 'b-', lw=2, label='Straight: $\\|\\mathbf{\\Omega}\\| \\equiv 0$')
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_bent], 'g--', lw=2, label='Bent: Pitching Reorientation')
    ax4.plot(t_vals, [d['Omega_mag'] for d in traj_chiral], 'r-', lw=2, label='Chiral: Continuous Autorotation')
    ax4.set_xlabel('Dimensionless Time $t / \\tau_c$')
    ax4.set_ylabel(r'Angular Velocity Magnitude $\|\mathbf{\Omega}\|$')
    ax4.set_title('Induced Gravitational Autorotation', fontweight='bold')
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend()

    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "robotic_arm_dynamic_simulation_trajectories.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"Saved 3D trajectory comparison figure: {plot_path}")


def main():
    print("\n" + "="*70)
    print("ROBOTIC ARM DYNAMIC SEDIMENTATION TRAJECTORY SIMULATION")
    print("="*70 + "\n")

    n_steps = 60
    dt = 1.0

    # 1. Simulate Straight Arm
    traj_straight, _, _, _ = simulate_arm_trajectory(config='straight', bend_deg=0.0, twist_deg=0.0,
                                                     link_resolution=12, dt=dt, n_steps=n_steps)

    # 2. Simulate Planar Bent Arm (45 deg)
    traj_bent, _, _, _ = simulate_arm_trajectory(config='bent', bend_deg=45.0, twist_deg=0.0,
                                                 link_resolution=12, dt=dt, n_steps=n_steps)

    # 3. Simulate Chiral Twisted Arm (45 deg bend, 45 deg twist)
    traj_chiral, _, _, _ = simulate_arm_trajectory(config='chiral', bend_deg=45.0, twist_deg=45.0,
                                                   link_resolution=12, dt=dt, n_steps=n_steps)

    # Export multi-resolution meshes for 3D viewer
    print("\nExporting multi-resolution mesh data (1, 12, 42 blobs/link)...")
    mesh_data = export_multi_resolution_meshes()

    # Save complete simulation dataset
    json_path = os.path.join(PROC_DIR, "robotic_arm_simulation_data.json")
    payload = {
        'metadata': {
            'N_links': 7,
            'link_spacing': 2.5,
            'link_radius': 1.0,
            'dt': dt,
            'n_steps': n_steps,
            'total_time': dt * n_steps,
            'Fz': 1.0,
            'eta': 1.0
        },
        'meshes': mesh_data,
        'trajectory_straight': traj_straight,
        'trajectory_bent': traj_bent,
        'trajectory_chiral': traj_chiral
    }

    with open(json_path, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"Saved comprehensive simulation JSON: {json_path}")

    # Plot 300 DPI publication trajectory figure
    print("\nGenerating publication figures...")
    plot_trajectory_comparison(traj_straight, traj_bent, traj_chiral)

    print("\n" + "="*70)
    print("Robotic Arm Dynamic Simulation Complete!")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
