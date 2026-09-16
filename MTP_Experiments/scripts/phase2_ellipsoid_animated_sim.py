"""
MTP_Experiments/scripts/phase2_ellipsoid_animated_sim.py

Generates an animated GIF of the dynamic low-Reynolds-number sedimentation
of a prolate ellipsoid (lambda = 2.0) in Stokes flow.
Visualizes:
  - 3D multiblob settling with non-tumbling pose stability (Omega = 0).
  - Oblique lateral drift (Ux > 0, Uz < 0).
  - 3-body race comparing theta = 0 deg (broadside), theta = 45 deg (gliding), and theta = 90 deg (streamlined).
  - Real-time telemetry HUD.
"""

import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import shutil

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

OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'ellipsoid')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')

def generate_animation():
    print("Loading precomputed simulation data...")
    json_path = os.path.join(PROC_DIR, "ellipsoid_simulation_data.json")
    with open(json_path, 'r') as f:
        data = json.load(f)

    traj_45 = data['trajectory_45']
    traj_0 = data['trajectory_0']
    traj_90 = data['trajectory_90']
    r_conf = np.array(data['r_conf_ref'])
    a_blob = data['a_blob']

    n_frames = len(traj_45)
    print(f"Preparing animation with {n_frames} frames...")

    fig = plt.figure(figsize=(16, 8), dpi=120)
    fig.patch.set_facecolor('#0d1117')

    # Grid layout: left is 3D view, right top is 2D XZ race, right bottom is telemetry
    ax_3d = fig.add_subplot(1, 2, 1, projection='3d', facecolor='#161b22')
    ax_xz = fig.add_subplot(2, 2, 2, facecolor='#161b22')
    ax_hud = fig.add_subplot(2, 2, 4, facecolor='#161b22')

    # Styling
    for ax in [ax_xz, ax_hud]:
        ax.tick_params(colors='#c9d1d9')
        for spine in ax.spines.values():
            spine.set_color('#30363d')
    ax_3d.tick_params(colors='#c9d1d9')
    ax_3d.xaxis.pane.set_edgecolor('#30363d')
    ax_3d.yaxis.pane.set_edgecolor('#30363d')
    ax_3d.zaxis.pane.set_edgecolor('#30363d')
    ax_3d.xaxis.pane.fill = False
    ax_3d.yaxis.pane.fill = False
    ax_3d.zaxis.pane.fill = False

    # Title
    fig.suptitle("MTP: Stokes Flow Sedimentation of Prolate Ellipsoid (λ = 2.0)",
                 color='#58a6ff', fontsize=16, fontweight='bold', y=0.96)

    # 45 deg orientation
    q_45 = Quaternion([traj_45[0]['q0'], traj_45[0]['q1'], traj_45[0]['q2'], traj_45[0]['q3']])
    R_45 = q_45.rotation_matrix()
    r_body_rot = np.dot(r_conf, R_45.T)

    # Animation update function
    def update(frame):
        ax_3d.cla()
        ax_xz.cla()
        ax_hud.cla()

        pt45 = traj_45[frame]
        pt0 = traj_0[frame]
        pt90 = traj_90[frame]

        pos45 = np.array([pt45['x'], pt45['y'], pt45['z']])
        r_current = r_body_rot + pos45

        # ----------------------------------------------------
        # 1. 3D View
        # ----------------------------------------------------
        # Draw particle blobs
        ax_3d.scatter(r_current[:, 0], r_current[:, 1], r_current[:, 2],
                      color='#388bfd', edgecolors='#58a6ff', s=25, alpha=0.9, label='Prolate Blobs (N=162)')

        # Draw trajectory history trail
        past_x = [traj_45[i]['x'] for i in range(frame + 1)]
        past_y = [traj_45[i]['y'] for i in range(frame + 1)]
        past_z = [traj_45[i]['z'] for i in range(frame + 1)]
        ax_3d.plot(past_x, past_y, past_z, color='#f78166', linewidth=2.5, label='Drift Path')

        # Draw Force arrow (Gravity, downwards)
        scale_f = 0.5
        ax_3d.quiver(pos45[0], pos45[1], pos45[2],
                     0.0, 0.0, -scale_f,
                     color='#ff7b72', arrow_length_ratio=0.25, linewidth=2.2, label='Gravity F_z')

        # Draw Velocity arrow (Oblique drift)
        scale_v = 12.0
        ax_3d.quiver(pos45[0], pos45[1], pos45[2],
                     pt45['Ux']*scale_v, pt45['Uy']*scale_v, pt45['Uz']*scale_v,
                     color='#7ee787', arrow_length_ratio=0.25, linewidth=2.2, label='Velocity U')

        ax_3d.set_xlim(-0.5, 0.5)
        ax_3d.set_ylim(-0.5, 0.5)
        ax_3d.set_zlim(-1.8, 0.2)
        ax_3d.set_xlabel('X (Lateral Drift)', color='#c9d1d9', labelpad=8)
        ax_3d.set_ylabel('Y', color='#c9d1d9', labelpad=8)
        ax_3d.set_zlabel('Z (Depth)', color='#c9d1d9', labelpad=8)
        ax_3d.set_title(f"3D Multiblob Sedimentation (t = {pt45['time']:.1f} s)", color='#f0f6fc', fontsize=12)
        ax_3d.view_init(elev=20, azim=45)
        ax_3d.legend(loc='upper right', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # ----------------------------------------------------
        # 2. X-Z Glide Path Race
        # ----------------------------------------------------
        # Full trajectories as faint background lines
        ax_xz.plot([d['x'] for d in traj_0], [d['z'] for d in traj_0], color='#8b949e', linestyle='--', alpha=0.5)
        ax_xz.plot([d['x'] for d in traj_45], [d['z'] for d in traj_45], color='#8b949e', linestyle='-', alpha=0.5)
        ax_xz.plot([d['x'] for d in traj_90], [d['z'] for d in traj_90], color='#8b949e', linestyle='-.', alpha=0.5)

        # Current positions of the 3 orientations
        ax_xz.plot([traj_0[i]['x'] for i in range(frame+1)], [traj_0[i]['z'] for i in range(frame+1)],
                   color='#79c0ff', linewidth=2, label='θ = 0° (Broadside)')
        ax_xz.scatter(pt0['x'], pt0['z'], color='#79c0ff', s=70, marker='^', zorder=5)

        ax_xz.plot(past_x, past_z,
                   color='#f78166', linewidth=2.5, label='θ = 45° (Gliding Drift)')
        ax_xz.scatter(pt45['x'], pt45['z'], color='#f78166', s=70, marker='o', zorder=5)

        ax_xz.plot([traj_90[i]['x'] for i in range(frame+1)], [traj_90[i]['z'] for i in range(frame+1)],
                   color='#d2a8ff', linewidth=2, label='θ = 90° (Streamlined)')
        ax_xz.scatter(pt90['x'], pt90['z'], color='#d2a8ff', s=70, marker='s', zorder=5)

        ax_xz.set_xlim(-0.02, 0.12)
        ax_xz.set_ylim(-1.8, 0.1)
        ax_xz.set_xlabel('Lateral Drift Displacement X', color='#c9d1d9', fontsize=10)
        ax_xz.set_ylabel('Vertical Settlement Z', color='#c9d1d9', fontsize=10)
        ax_xz.set_title("Glide Path Comparison (X-Z Plane)", color='#f0f6fc', fontsize=11)
        ax_xz.grid(True, linestyle=':', alpha=0.3, color='#8b949e')
        ax_xz.legend(loc='lower left', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # ----------------------------------------------------
        # 3. Real-Time Telemetry & Hydrodynamic HUD
        # ----------------------------------------------------
        ax_hud.axis('off')
        drift_angle = np.degrees(np.arctan(abs(pt45['Ux'] / pt45['Uz']))) if pt45['Uz'] != 0 else 0.0

        hud_text = (
            "HYDRODYNAMIC TELEMETRY HUD\n"
            "--------------------------------------------------\n"
            f"Time Elapsed       :  t = {pt45['time']:.2f} s\n"
            f"Vertical Altitude  :  Z = {pt45['z']:.4f}\n"
            f"Lateral Drift      :  X = +{pt45['x']:.4f}\n"
            f"Settling Speed     :  |Uz| = {abs(pt45['Uz']):.5f}\n"
            f"Lateral Velocity   :  Ux  = +{pt45['Ux']:.5f}\n"
            f"Glide Angle (α)    :  α  = {drift_angle:.2f}°\n"
            "--------------------------------------------------\n"
            f"Angular Speed |Ω|  :  {pt45['Omega_mag']:.2e} rad/s\n"
            "Pose Stability     :  NON-TUMBLING (Ω ≡ 0)\n"
            "Reynolds Number    :  Re < 10⁻⁴ (Stokes Flow)\n"
            "Hydrodynamic Regime:  Force-Velocity Linear\n"
            "Stokes-Perrin Model:  Coupled 6x6 RPY Kernel"
        )
        ax_hud.text(0.05, 0.95, hud_text, transform=ax_hud.transAxes,
                    fontsize=10.5, fontfamily='monospace', color='#58a6ff',
                    verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.8', facecolor='#0d1117', edgecolor='#30363d', alpha=0.9))

        return ax_3d, ax_xz, ax_hud

    print("Rendering animation frames...")
    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=80, blit=False)

    gif_path = os.path.join(PLOT_DIR, "ellipsoid_sedimentation_simulation.gif")
    anim.save(gif_path, writer='pillow', fps=12)
    plt.close()
    print(f"Saved animated simulation GIF: {gif_path}")

    # Copy to artifact directory so that Antigravity IDE can render it directly
    artifact_dir = r"C:\Users\Asus\.gemini\antigravity-ide\brain\0c831f51-dcfa-45e2-afe6-607ffdc616ff"
    artifact_gif = os.path.join(artifact_dir, "ellipsoid_sedimentation_simulation.gif")
    artifact_png = os.path.join(artifact_dir, "ellipsoid_dynamic_simulation_trajectories.png")
    
    shutil.copyfile(gif_path, artifact_gif)
    src_png = os.path.join(PLOT_DIR, "ellipsoid_dynamic_simulation_trajectories.png")
    if os.path.exists(src_png):
        shutil.copyfile(src_png, artifact_png)
    print(f"Copied visualization artifacts to {artifact_dir}")

if __name__ == '__main__':
    generate_animation()
