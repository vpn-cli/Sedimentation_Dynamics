"""
MTP_Experiments/scripts/phase2_ellipsoid_trajectory_sim.py

Dynamic Time-Stepping Sedimentation Simulation of an Ellipsoid in Stokes Flow.
Demonstrates:
  - Anisotropic sedimentation of a prolate spheroid (lambda = 2.0).
  - Tilted descent at theta = 45 deg (maximum oblique lateral drift).
  - Comparison with theta = 0 deg (streamlined vertical) and theta = 90 deg (broadside-on).
  - 3D particle pose tracking and hydrodynamic velocity vector visualization.
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
OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'ellipsoid')
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


def simulate_trajectory(N=162, a=2.0, b=1.0, eta=1.0, Fz=1.0, theta_deg=45.0, dt=0.5, n_steps=60):
    """
    Simulates time-stepping trajectory of an ellipsoid sedimenting under gravity.
    """
    r_conf, a_blob, meta = p2c.generate_ellipsoid_mesh(N, a, b)
    theta_rad = np.radians(theta_deg)

    # Initial state
    pos = np.array([0.0, 0.0, 0.0])
    # Orientation: tilt symmetry axis (x-axis) toward z-axis in xz-plane
    q = Quaternion.from_rotation(np.array([0.0, theta_rad, 0.0]))
    force_vec = np.array([0.0, 0.0, -Fz])

    traj = []
    current_time = 0.0

    print(f"Simulating Ellipsoid (lambda={a/b:.1f}, theta={theta_deg} deg, Fz={Fz}, N={N})...")

    for step in range(n_steps + 1):
        # Solve instantaneous Stokes sedimentation velocities
        U, Omega, _ = p2c.solve_sedimentation(r_conf, a_blob, eta, force_vec, orientation=q)

        traj.append({
            'step': step,
            'time': current_time,
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
            'q3': float(q.entries[3])
        })

        if step == n_steps:
            break

        # Time integration: Forward Euler for translation
        pos = pos + U * dt

        # Orientation update via quaternion kinematics (dphi = Omega * dt)
        # Note: In unbounded flow, Omega is zero within machine precision.
        phi = Omega * dt
        q_rot = Quaternion.from_rotation(phi)
        q = q * q_rot
        current_time += dt

    return traj, r_conf, a_blob, meta


def main():
    print("\n" + "="*70)
    print("ELLIPSOID DYNAMIC SEDIMENTATION SIMULATION")
    print("="*70 + "\n")

    # Run primary simulation at theta = 45 deg (maximum lateral drift)
    traj_45, r_conf, a_blob, meta = simulate_trajectory(N=162, a=2.0, b=1.0, theta_deg=45.0, dt=1.0, n_steps=40)

    # Run comparisons at theta = 0 deg (vertical streamlined) and theta = 90 deg (broadside-on)
    traj_0, _, _, _ = simulate_trajectory(N=162, a=2.0, b=1.0, theta_deg=0.0, dt=1.0, n_steps=40)
    traj_90, _, _, _ = simulate_trajectory(N=162, a=2.0, b=1.0, theta_deg=90.0, dt=1.0, n_steps=40)

    # Save trajectory CSV
    csv_file = os.path.join(RAW_DIR, "ellipsoid_trajectory_sim.csv")
    header = "step,time,x,y,z,Ux,Uy,Uz,U_mag,Omega_mag,q0,q1,q2,q3"
    lines = [header]
    for d in traj_45:
        lines.append(f"{d['step']},{d['time']:.2f},{d['x']:.6f},{d['y']:.6f},{d['z']:.6f},"
                     f"{d['Ux']:.8e},{d['Uy']:.8e},{d['Uz']:.8e},{d['U_mag']:.8e},{d['Omega_mag']:.8e},"
                     f"{d['q0']:.6f},{d['q1']:.6f},{d['q2']:.6f},{d['q3']:.6f}")
    with open(csv_file, 'w') as f:
        f.write("\n".join(lines))
    print(f"\nSaved simulation CSV: {csv_file}")

    # Export trajectory JSON for interactive 3D web simulation
    json_payload = {
        'metadata': meta,
        'trajectory_45': traj_45,
        'trajectory_0': traj_0,
        'trajectory_90': traj_90,
        'r_conf_ref': r_conf.tolist(),
        'a_blob': a_blob
    }
    json_file = os.path.join(PROC_DIR, "ellipsoid_simulation_data.json")
    with open(json_file, 'w') as f:
        json.dump(json_payload, f, indent=2)
    print(f"Saved simulation JSON: {json_file}")

    # =========================================================================
    # MULTI-PANEL PUBLICATION VISUALIZATION
    # =========================================================================
    print("Generating comprehensive simulation visualization figure...")
    fig = plt.figure(figsize=(15, 11))

    # Panel 1: 3D Trajectory & Particle Pose Tracking
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    x_45 = [d['x'] for d in traj_45]
    y_45 = [d['y'] for d in traj_45]
    z_45 = [d['z'] for d in traj_45]

    ax1.plot(x_45, y_45, z_45, 'r-', linewidth=2.5, label=r'Centroid Path ($\theta=45^\circ$)')

    # Draw particle multiblob at 5 snapshot intervals along descent
    snapshot_indices = [0, 10, 20, 30, 40]
    colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#d62728']

    for idx, col in zip(snapshot_indices, colors):
        pt = traj_45[idx]
        pos_snap = np.array([pt['x'], pt['y'], pt['z']])
        q_snap = Quaternion([pt['q0'], pt['q1'], pt['q2'], pt['q3']])
        R_snap = q_snap.rotation_matrix()
        r_snap = np.dot(r_conf, R_snap.T) + pos_snap

        # Plot blobs of the ellipsoid
        ax1.scatter(r_snap[:, 0], r_snap[:, 1], r_snap[:, 2], color=col, s=15, alpha=0.75)
        # Draw velocity arrow
        scale_vec = 40.0
        ax1.quiver(pos_snap[0], pos_snap[1], pos_snap[2],
                   pt['Ux']*scale_vec, pt['Uy']*scale_vec, pt['Uz']*scale_vec,
                   color='black', arrow_length_ratio=0.3, linewidth=1.5)

    ax1.set_xlabel('X (Lateral Drift)')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z (Vertical Depth)')
    ax1.set_title(r'(A) 3D Sedimentation Trajectory & Poses ($\theta=45^\circ$)', fontweight='bold')
    ax1.view_init(elev=18, azim=40)
    ax1.legend(loc='upper right')

    # Panel 2: X-Z Side View showing Oblique Drift Angle
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(x_45, z_45, 'r-', linewidth=2.5, label=r'$\theta = 45^\circ$ (Max Lateral Drift)')
    ax2.plot([d['x'] for d in traj_0], [d['z'] for d in traj_0], 'b--', linewidth=2, label=r'$\theta = 0^\circ$ (Streamlined Vertical)')
    ax2.plot([d['x'] for d in traj_90], [d['z'] for d in traj_90], 'g-.', linewidth=2, label=r'$\theta = 90^\circ$ (Broadside Vertical)')

    # Annotate drift angle on 45 deg trajectory
    drift_slope = traj_45[-1]['x'] / traj_45[-1]['z']
    drift_angle_deg = np.degrees(np.arctan(abs(traj_45[-1]['x'] / traj_45[-1]['z'])))
    ax2.annotate(rf'Oblique Drift Angle $\alpha \approx {drift_angle_deg:.2f}^\circ$',
                 xy=(traj_45[25]['x'], traj_45[25]['z']),
                 xytext=(traj_45[25]['x'] + 0.02, traj_45[25]['z'] + 0.3),
                 arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
                 fontsize=10, fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", fc="yellow", alpha=0.3))

    ax2.set_xlabel('Lateral Drift Displacement X')
    ax2.set_ylabel('Vertical Settlement Z')
    ax2.set_title(r'(B) 2D Glide Path in $X-Z$ Plane', fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()

    # Panel 3: Time Series of Displacements Z(t) and X(t)
    ax3 = fig.add_subplot(2, 2, 3)
    t_vals = [d['time'] for d in traj_45]
    ax3.plot(t_vals, [d['z'] for d in traj_45], 'r-', linewidth=2, label=r'$Z(t)$ ($\theta=45^\circ$)')
    ax3.plot(t_vals, [d['z'] for d in traj_0], 'b--', linewidth=2, label=r'$Z(t)$ ($\theta=0^\circ$)')
    ax3.plot(t_vals, [d['z'] for d in traj_90], 'g-.', linewidth=2, label=r'$Z(t)$ ($\theta=90^\circ$)')
    ax3.plot(t_vals, [d['x'] for d in traj_45], 'm:', linewidth=2.5, label=r'$X(t)$ Lateral Drift ($\theta=45^\circ$)')

    ax3.set_xlabel('Time t')
    ax3.set_ylabel('Position Coordinate')
    ax3.set_title('(C) Kinematic Displacements vs Time', fontweight='bold')
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend()

    # Panel 4: Orientation Stability & Zero Angular Velocity
    ax4 = fig.add_subplot(2, 2, 4)
    omega_vals = [d['Omega_mag'] for d in traj_45]
    ax4.plot(t_vals, omega_vals, 'k.-', linewidth=1.5, label=r'$\|\mathbf{\Omega}\|(t)$')
    ax4.set_xlabel('Time t')
    ax4.set_ylabel(r'Angular Velocity Magnitude $\|\mathbf{\Omega}\|$')
    ax4.set_title(r'(D) Angular Velocity Stability ($\|\mathbf{\Omega}\| \equiv 0$)', fontweight='bold')
    ax4.set_ylim(-1e-18, 5e-18)
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend()

    plt.tight_layout()
    plot_file = os.path.join(PLOT_DIR, "ellipsoid_dynamic_simulation_trajectories.png")
    plt.savefig(plot_file)
    plt.close()
    print(f"Saved 4-panel simulation figure: {plot_file}")

    print("\n" + "="*70)
    print("Ellipsoid Simulation Complete!")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
