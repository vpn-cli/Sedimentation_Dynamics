"""
MTP_Experiments/scripts/phase2_robotic_arm_animated_sim.py

Generates an animated GIF of the dynamic low-Reynolds-number sedimentation
of a 7-segment robotic arm in Stokes flow.
Visualizes:
  - 3D multiblob settling and articulated link geometry.
  - Side-by-side dynamic race:
      1. Straight Arm (steady vertical settling, Omega = 0).
      2. Planar Bent Arm (pitching reorientation).
      3. Chiral Arm (continuous 3D helical spiral autorotation).
  - Real-time telemetry HUD displaying velocities and rotation rates.
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

OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'robotic_arm')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')


def generate_animation():
    print("Loading precomputed robotic arm simulation data...")
    json_path = os.path.join(PROC_DIR, "robotic_arm_simulation_data.json")
    with open(json_path, 'r') as f:
        data = json.load(f)

    traj_straight = data['trajectory_straight']
    traj_bent = data['trajectory_bent']
    traj_chiral = data['trajectory_chiral']

    n_frames = len(traj_chiral)
    print(f"Preparing animation with {n_frames} frames...")

    fig = plt.figure(figsize=(16, 8), dpi=120)
    fig.patch.set_facecolor('#0d1117')

    ax_3d = fig.add_subplot(1, 2, 1, projection='3d', facecolor='#161b22')
    ax_xy = fig.add_subplot(2, 2, 2, facecolor='#161b22')
    ax_hud = fig.add_subplot(2, 2, 4, facecolor='#161b22')

    # Style 3D axis
    ax_3d.xaxis.pane.fill = False
    ax_3d.yaxis.pane.fill = False
    ax_3d.zaxis.pane.fill = False
    ax_3d.xaxis.pane.set_edgecolor('#30363d')
    ax_3d.yaxis.pane.set_edgecolor('#30363d')
    ax_3d.zaxis.pane.set_edgecolor('#30363d')
    ax_3d.tick_params(colors='#8b949e', labelsize=8)

    # Subplot 2: 2D XY drift & descent
    ax_xy.tick_params(colors='#8b949e', labelsize=8)
    ax_xy.grid(True, linestyle=':', color='#30363d', alpha=0.7)
    for spine in ax_xy.spines.values():
        spine.set_color('#30363d')

    # Subplot 3: HUD
    ax_hud.axis('off')

    def update(frame):
        ax_3d.cla()
        ax_xy.cla()
        ax_hud.cla()
        ax_hud.axis('off')

        # Current steps
        d_s = traj_straight[frame]
        d_b = traj_bent[frame]
        d_c = traj_chiral[frame]

        current_time = d_c['time']

        # 1. Plot 3D trajectories
        xs_s = [d['x'] for d in traj_straight[:frame+1]]
        ys_s = [d['y'] for d in traj_straight[:frame+1]]
        zs_s = [d['z'] for d in traj_straight[:frame+1]]

        xs_b = [d['x'] for d in traj_bent[:frame+1]]
        ys_b = [d['y'] for d in traj_bent[:frame+1]]
        zs_b = [d['z'] for d in traj_bent[:frame+1]]

        xs_c = [d['x'] for d in traj_chiral[:frame+1]]
        ys_c = [d['y'] for d in traj_chiral[:frame+1]]
        zs_c = [d['z'] for d in traj_chiral[:frame+1]]

        ax_3d.plot(xs_s, ys_s, zs_s, color='#58a6ff', lw=1.8, label='Straight Arm')
        ax_3d.plot(xs_b, ys_b, zs_b, color='#3fb950', lw=1.8, label='Bent Arm 45°')
        ax_3d.plot(xs_c, ys_c, zs_c, color='#f85149', lw=2.2, label='Chiral 3D Spiral')

        # Draw arm skeletons (links)
        links_s = np.array(d_s['link_centers'])
        links_b = np.array(d_b['link_centers'])
        links_c = np.array(d_c['link_centers'])

        ax_3d.plot(links_s[:, 0], links_s[:, 1], links_s[:, 2], 'o-', color='#58a6ff', lw=3.0, ms=6)
        ax_3d.plot(links_b[:, 0], links_b[:, 1], links_b[:, 2], 's-', color='#3fb950', lw=3.0, ms=6)
        ax_3d.plot(links_c[:, 0], links_c[:, 1], links_c[:, 2], '^-', color='#f85149', lw=3.5, ms=7)

        ax_3d.set_title(f"3D Sedimentation Dynamics (t = {current_time:.1f})", color='#c9d1d9', fontsize=12, pad=12)
        ax_3d.set_xlabel("X (Drift)", color='#8b949e', fontsize=9)
        ax_3d.set_ylabel("Y (Drift)", color='#8b949e', fontsize=9)
        ax_3d.set_zlabel("Z (Settling)", color='#8b949e', fontsize=9)
        ax_3d.set_xlim([-10, 10])
        ax_3d.set_ylim([-5, 5])
        ax_3d.set_zlim([-1.2, 0.2])
        ax_3d.legend(loc='lower left', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # 2. Plot 2D Horizontal Drift & Trajectory
        ax_xy.plot(xs_s, ys_s, color='#58a6ff', lw=1.5, label='Straight')
        ax_xy.plot(xs_b, ys_b, color='#3fb950', lw=1.5, label='Bent 45°')
        ax_xy.plot(xs_c, ys_c, color='#f85149', lw=2.0, label='Chiral Spiral')
        ax_xy.scatter([d_s['x']], [d_s['y']], color='#58a6ff', s=35)
        ax_xy.scatter([d_b['x']], [d_b['y']], color='#3fb950', s=35)
        ax_xy.scatter([d_c['x']], [d_c['y']], color='#f85149', s=45)

        ax_xy.set_title("Horizontal Trajectory Projection (X vs Y)", color='#c9d1d9', fontsize=11)
        ax_xy.set_xlabel("X Position", color='#8b949e', fontsize=9)
        ax_xy.set_ylabel("Y Position", color='#8b949e', fontsize=9)
        ax_xy.tick_params(colors='#8b949e', labelsize=8)
        ax_xy.grid(True, linestyle=':', color='#30363d', alpha=0.7)
        ax_xy.legend(loc='upper right', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # 3. Telemetry HUD
        hud_text = (
            f"STOKESIAN TELEMETRY HUD\n"
            f"--------------------------------------------------\n"
            f"Step: {d_c['step']:02d} / {n_frames-1:02d}    |  Sim Time: {current_time:5.1f} tau_c\n\n"
            f"[1] STRAIGHT ARM (Symmetric Rod):\n"
            f"    Altitude Z : {d_s['z']:+7.4f}  |  Speed |Uz|: {abs(d_s['Uz']):.5f}\n"
            f"    Rotation   : ||Omega|| = {d_s['Omega_mag']:.2e} rad/s (NON-TUMBLING)\n\n"
            f"[2] BENT ARM 45° (Planar Asymmetric):\n"
            f"    Altitude Z : {d_b['z']:+7.4f}  |  Speed |Uz|: {abs(d_b['Uz']):.5f}\n"
            f"    Rotation   : ||Omega|| = {d_b['Omega_mag']:.4e} rad/s (PITCHING)\n\n"
            f"[3] CHIRAL ARM 45° (3D Twisted Chain):\n"
            f"    Altitude Z : {d_c['z']:+7.4f}  |  Drift X,Y: [{d_c['x']:+.3f}, {d_c['y']:+.3f}]\n"
            f"    Rotation   : ||Omega|| = {d_c['Omega_mag']:.4e} rad/s (3D AUTOROTATION)\n"
            f"    Descent    : HELICAL / SPIRAL TRAJECTORY\n"
        )
        ax_hud.text(0.05, 0.95, hud_text, transform=ax_hud.transAxes,
                    fontfamily='monospace', fontsize=9, verticalalignment='top',
                    color='#58a6ff',
                    bbox=dict(boxstyle='round,pad=0.6', facecolor='#161b22', edgecolor='#30363d', alpha=0.9))

        return ax_3d, ax_xy, ax_hud

    print("Rendering animation frames...")
    ani = animation.FuncAnimation(fig, update, frames=range(0, n_frames, 2), interval=100, blit=False)

    gif_path = os.path.join(PLOT_DIR, "robotic_arm_sedimentation_simulation.gif")
    print(f"Saving GIF to {gif_path} (this may take a few moments)...")
    ani.save(gif_path, writer='pillow', fps=12)
    plt.close()
    print(f"Animated simulation successfully saved to: {gif_path}")

    # Copy to artifacts directory
    artifact_dir = r"C:\Users\Asus\.gemini\antigravity-ide\brain\3b9ac79f-de7a-4c76-8fb6-545b94769604"
    if os.path.exists(artifact_dir):
        dest_gif = os.path.join(artifact_dir, "robotic_arm_sedimentation_simulation.gif")
        shutil.copy(gif_path, dest_gif)
        print(f"Copied animated GIF to artifact directory: {dest_gif}")


if __name__ == '__main__':
    generate_animation()
