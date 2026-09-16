"""
MTP_Experiments/scripts/phase2_boomerang_animated_sim.py

Generates an animated GIF of the dynamic low-Reynolds-number sedimentation
of a boomerang particle in Stokes flow.
Visualizes:
  - 3D multiblob settling and articulated L-shape geometry.
  - Side-by-side dynamic race:
      1. Apex Down (steady edge-on glide, Omega ≈ 0).
      2. Flat Pose (pitching reorientation).
      3. Tilted 45° (oblique lateral drift).
      4. Chiral Boomerang (continuous 3D helical spiral autorotation).
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

OUT_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical', 'boomerang')
PLOT_DIR = os.path.join(OUT_DIR, 'plots')
PROC_DIR = os.path.join(OUT_DIR, 'processed_data')


def generate_animation():
    print("Loading precomputed boomerang simulation data...")
    json_path = os.path.join(PROC_DIR, "boomerang_simulation_data.json")
    with open(json_path, 'r') as f:
        data = json.load(f)

    traj_edge = data['trajectory_edge']
    traj_flat = data['trajectory_flat']
    traj_tilted = data['trajectory_tilted']
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

    # Style 2D axis
    ax_xy.tick_params(colors='#8b949e', labelsize=8)
    ax_xy.grid(True, linestyle=':', color='#30363d', alpha=0.7)
    for spine in ax_xy.spines.values():
        spine.set_color('#30363d')

    ax_hud.axis('off')

    def update(frame):
        ax_3d.cla()
        ax_xy.cla()
        ax_hud.cla()
        ax_hud.axis('off')

        d_e = traj_edge[frame]
        d_f = traj_flat[frame]
        d_t = traj_tilted[frame]
        d_c = traj_chiral[frame]

        current_time = d_c['time']

        # 1. 3D Trajectories
        def get_sub_xyz(traj):
            return [d['x'] for d in traj[:frame+1]], [d['y'] for d in traj[:frame+1]], [d['z'] for d in traj[:frame+1]]

        xe, ye, ze = get_sub_xyz(traj_edge)
        xf, yf, zf = get_sub_xyz(traj_flat)
        xt, yt, zt = get_sub_xyz(traj_tilted)
        xc, yc, zc = get_sub_xyz(traj_chiral)

        ax_3d.plot(xe, ye, ze, color='#58a6ff', lw=1.8, label='Apex Down (Gliding)')
        ax_3d.plot(xf, yf, zf, color='#3fb950', lw=1.8, label='Flat (Pitching)')
        ax_3d.plot(xt, yt, zt, color='#d29922', lw=1.8, label='Tilted 45° (Drift)')
        ax_3d.plot(xc, yc, zc, color='#f85149', lw=2.2, label='Chiral Spiral')

        # Draw current particle blobs
        blobs_c = np.array(d_c['blobs_lab'])
        ax_3d.scatter(blobs_c[:, 0], blobs_c[:, 1], blobs_c[:, 2], color='#f85149', s=30, alpha=0.85)

        blobs_t = np.array(d_t['blobs_lab'])
        ax_3d.scatter(blobs_t[:, 0], blobs_t[:, 1], blobs_t[:, 2], color='#d29922', s=25, alpha=0.6)

        ax_3d.set_title(f"3D Sedimentation Dynamics (t = {current_time:.1f} τ_c)", color='#c9d1d9', fontsize=12, pad=12)
        ax_3d.set_xlabel("X (Drift)", color='#8b949e', fontsize=9)
        ax_3d.set_ylabel("Y (Drift)", color='#8b949e', fontsize=9)
        ax_3d.set_zlabel("Z (Settling)", color='#8b949e', fontsize=9)
        ax_3d.set_xlim([-4, 4])
        ax_3d.set_ylim([-4, 4])
        ax_3d.set_zlim([-4.0, 0.5])
        ax_3d.legend(loc='lower left', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # 2. Horizontal Projection
        ax_xy.plot(xe, ye, color='#58a6ff', lw=1.5, label='Apex Down')
        ax_xy.plot(xf, yf, color='#3fb950', lw=1.5, label='Flat')
        ax_xy.plot(xt, yt, color='#d29922', lw=1.5, label='Tilted 45°')
        ax_xy.plot(xc, yc, color='#f85149', lw=2.0, label='Chiral Spiral')
        ax_xy.scatter([d_e['x']], [d_e['y']], color='#58a6ff', s=35)
        ax_xy.scatter([d_f['x']], [d_f['y']], color='#3fb950', s=35)
        ax_xy.scatter([d_t['x']], [d_t['y']], color='#d29922', s=35)
        ax_xy.scatter([d_c['x']], [d_c['y']], color='#f85149', s=45)

        ax_xy.set_title("Horizontal Trajectory Projection (X vs Y)", color='#c9d1d9', fontsize=11)
        ax_xy.set_xlabel("X Position", color='#8b949e', fontsize=9)
        ax_xy.set_ylabel("Y Position", color='#8b949e', fontsize=9)
        ax_xy.tick_params(colors='#8b949e', labelsize=8)
        ax_xy.grid(True, linestyle=':', color='#30363d', alpha=0.7)
        ax_xy.legend(loc='upper right', facecolor='#161b22', edgecolor='#30363d', labelcolor='#c9d1d9', fontsize=8)

        # 3. Telemetry HUD
        hud_text = (
            f"BOOMERANG STOKESIAN TELEMETRY HUD\n"
            f"--------------------------------------------------\n"
            f"Step: {d_c['step']:02d} / {n_frames-1:02d}    |  Sim Time: {current_time:5.1f} tau_c\n\n"
            f"[1] APEX DOWN (Edge-On Gliding):\n"
            f"    Altitude Z : {d_e['z']:+7.4f}  |  Speed |Uz|: {abs(d_e['Uz']):.5f}\n"
            f"    Rotation   : ||Omega|| = {d_e['Omega_mag']:.2e} rad/s (STEADY GLIDE)\n\n"
            f"[2] FLAT POSE (Horizontal XY Plane):\n"
            f"    Altitude Z : {d_f['z']:+7.4f}  |  Speed |Uz|: {abs(d_f['Uz']):.5f}\n"
            f"    Rotation   : ||Omega|| = {d_f['Omega_mag']:.4e} rad/s (PITCHING)\n\n"
            f"[3] TILTED 45° (Oblique Descent):\n"
            f"    Altitude Z : {d_t['z']:+7.4f}  |  Drift X,Y: [{d_t['x']:+.3f}, {d_t['y']:+.3f}]\n"
            f"    Descent    : LATERAL DRIFT + REORIENTATION\n\n"
            f"[4] CHIRAL BOOMERANG (15° Dihedral Twist):\n"
            f"    Altitude Z : {d_c['z']:+7.4f}  |  Rotation: {d_c['Omega_mag']:.4e} rad/s\n"
            f"    Descent    : 3D HELICAL SPIRAL AUTOROTATION\n"
        )
        ax_hud.text(0.05, 0.95, hud_text, transform=ax_hud.transAxes,
                    fontfamily='monospace', fontsize=9, verticalalignment='top',
                    color='#58a6ff',
                    bbox=dict(boxstyle='round,pad=0.6', facecolor='#161b22', edgecolor='#30363d', alpha=0.9))

        return ax_3d, ax_xy, ax_hud

    print("Rendering animation frames...")
    ani = animation.FuncAnimation(fig, update, frames=range(0, n_frames, 2), interval=100, blit=False)

    gif_path = os.path.join(PLOT_DIR, "boomerang_sedimentation_simulation.gif")
    print(f"Saving GIF to {gif_path} (this may take a few moments)...")
    ani.save(gif_path, writer='pillow', fps=12)
    plt.close()
    print(f"Animated simulation successfully saved to: {gif_path}")

    # Copy to artifacts directory
    artifact_dir = r"C:\Users\Asus\.gemini\antigravity-ide\brain\3b9ac79f-de7a-4c76-8fb6-545b94769604"
    if os.path.exists(artifact_dir):
        dest_gif = os.path.join(artifact_dir, "boomerang_sedimentation_simulation.gif")
        shutil.copy(gif_path, dest_gif)
        print(f"Copied animated GIF to artifact directory: {dest_gif}")


if __name__ == '__main__':
    generate_animation()
