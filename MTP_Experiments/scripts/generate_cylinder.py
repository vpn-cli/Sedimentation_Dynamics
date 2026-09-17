"""
MTP_Experiments/scripts/generate_cylinder.py

Parametric generator for rigid circular cylinder multiblob configurations.
Generates an axisymmetric, closed circular cylinder surface mesh with exact
center-of-mass at (0, 0, 0) and symmetry axis aligned with the z-axis.

Dimensions:
- Radius: R = 1.0
- Length: L = 2.0 (Aspect ratio L / (2R) = 1.0, equidimensional cylinder)

Multiblob Surface Representation:
- Total surface area: A_total = 2*pi*R*L + 2*pi*R^2 = 4*pi + 2*pi = 6*pi
- Lateral barrel surface (r = R, z in [-L/2, L/2]): Area = 4*pi (2/3 of total)
- Top end-cap disk (z = +L/2, r <= R): Area = pi (1/6 of total)
- Bottom end-cap disk (z = -L/2, r <= R): Area = pi (1/6 of total)

Refinement Hierarchy:
- Exact MTP sequence: N in {12, 42, 162, 642, 2562}
- Since each N is divisible by 6, the blob distribution maintains identical surface density:
    N_lat = 2*N // 3
    N_cap_top = N // 6
    N_cap_bot = N // 6
"""

import os
import sys
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

def generate_cylinder_multiblob(N, R_geom=1.0, L_geom=2.0, overlap_factor=0.50):
    """
    Generates an exact N-blob configuration on the closed surface of a circular cylinder.
    
    Parameters:
    -----------
    N : int
        Total number of blobs (must be divisible by 6 for exact area partitioning).
    R_geom : float
        Cylinder radius (default 1.0).
    L_geom : float
        Cylinder total length (default 2.0).
    overlap_factor : float
        Ratio of blob radius to characteristic grid spacing h = sqrt(6*pi / N).
        
    Returns:
    --------
    coords : np.ndarray, shape (N, 3)
    a_blob : float
    metadata : dict
    """
    if N % 6 != 0:
        raise ValueError(f"N must be divisible by 6 for exact 2/3 : 1/6 : 1/6 partitioning, got N={N}")
        
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    N_lat = int(round(2 * N / 3))
    N_cap = int(round(N / 6))
    
    # 1. Lateral cylindrical surface (barrel)
    # Cylindrical Fibonacci phyllotaxis spiral
    j = np.arange(N_lat, dtype=np.float64)
    z_lat = -L_geom / 2.0 + (L_geom / N_lat) * (j + 0.5)
    theta_lat = j * golden_angle
    x_lat = R_geom * np.cos(theta_lat)
    y_lat = R_geom * np.sin(theta_lat)
    
    # 2. Top flat circular end-cap disk at z = +L/2
    # Vogel golden-angle spiral
    k_top = np.arange(N_cap, dtype=np.float64)
    r_top = R_geom * np.sqrt((k_top + 0.5) / N_cap)
    theta_top = k_top * golden_angle
    x_top = r_top * np.cos(theta_top)
    y_top = r_top * np.sin(theta_top)
    z_top = np.full(N_cap, L_geom / 2.0)
    
    # 3. Bottom flat circular end-cap disk at z = -L/2
    # Vogel golden-angle spiral rotated by pi for geometric reflection symmetry
    k_bot = np.arange(N_cap, dtype=np.float64)
    r_bot = R_geom * np.sqrt((k_bot + 0.5) / N_cap)
    theta_bot = k_bot * golden_angle + np.pi
    x_bot = r_bot * np.cos(theta_bot)
    y_bot = r_bot * np.sin(theta_bot)
    z_bot = np.full(N_cap, -L_geom / 2.0)
    
    # Assemble complete configuration
    x = np.concatenate([x_lat, x_top, x_bot])
    y = np.concatenate([y_lat, y_top, y_bot])
    z = np.concatenate([z_lat, z_top, z_bot])
    
    # Center of mass correction to exact origin (machine precision)
    com_raw = np.array([np.mean(x), np.mean(y), np.mean(z)])
    x -= com_raw[0]
    y -= com_raw[1]
    z -= com_raw[2]
    
    coords = np.column_stack([x, y, z])
    
    # Characteristic spacing based on uniform surface area density
    # Total surface area A = 2*pi*R*L + 2*pi*R^2 = 6*pi for R=1, L=2
    A_total = 2.0 * np.pi * R_geom * L_geom + 2.0 * np.pi * (R_geom**2)
    h_char = np.sqrt(A_total / N)
    a_blob = float(overlap_factor * h_char)
    
    # Inter-blob distances
    d = pdist(coords)
    min_d = float(np.min(d))
    mean_d = float(np.mean(d))
    
    # Radius of gyration
    Rg = float(np.sqrt(np.mean(np.sum(coords**2, axis=1))))
    # Theoretical Rg for surface shell of circular cylinder with flat caps:
    # <r^2> = (4*pi * (1 + 1/3) + 2*pi * (1/2 + 1)) / (6*pi) = (16/3 + 3)/6 = 25/18
    # Rg_th = sqrt(25/18) = 5*sqrt(2)/6 approx 1.1785113
    Rg_th = float(np.sqrt(25.0 / 18.0) * R_geom)
    
    # Inertia tensor components (assuming equal point masses)
    Ixx = float(np.sum(y**2 + z**2))
    Iyy = float(np.sum(x**2 + z**2))
    Izz = float(np.sum(x**2 + y**2))
    Ixy = float(-np.sum(x * y))
    Ixz = float(-np.sum(x * z))
    Iyz = float(-np.sum(y * z))
    
    metadata = {
        'shape': 'cylinder',
        'discretization': 'surface_shell_spiral',
        'N_blobs': int(N),
        'N_lat': int(N_lat),
        'N_cap_top': int(N_cap),
        'N_cap_bot': int(N_cap),
        'R_geom': float(R_geom),
        'L_geom': float(L_geom),
        'aspect_ratio': float(L_geom / (2.0 * R_geom)),
        'blob_radius': float(a_blob),
        'characteristic_h': float(h_char),
        'min_distance': float(min_d),
        'mean_distance': float(mean_d),
        'R_g': float(Rg),
        'R_g_theory': float(Rg_th),
        'R_g_error_pct': float(abs(Rg - Rg_th) / Rg_th * 100.0),
        'Ixx': float(Ixx),
        'Iyy': float(Iyy),
        'Izz': float(Izz),
        'Ixx_over_Iyy': float(Ixx / Iyy),
        'Ixy_over_Ixx': float(Ixy / Ixx),
        'Ixz_over_Ixx': float(Ixz / Ixx),
        'Iyz_over_Ixx': float(Iyz / Ixx),
        'center_of_mass': coords.mean(axis=0).tolist(),
        'com_offset_norm': float(np.linalg.norm(coords.mean(axis=0)))
    }
    
    return coords, a_blob, metadata

def save_cylinder_vertex_file(filepath, coords, a_blob, metadata=None):
    """
    Saves coordinates to standard .vertex file format.
    Line 1: N_blobs  blob_spacing (2 * a_blob)
    Remaining lines: x  y  z
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    N = len(coords)
    spacing = 2.0 * a_blob
    with open(filepath, 'w') as fp:
        if metadata:
            fp.write(f"# Circular Cylinder Multiblob Configuration: N={N}, R={metadata.get('R_geom', 1.0):.4f}, L={metadata.get('L_geom', 2.0):.4f}\n")
            fp.write(f"# a_blob={a_blob:.6f}, spacing={spacing:.6f}, R_g={metadata.get('R_g', 0.0):.6f} (Theory {metadata.get('R_g_theory', 0.0):.6f})\n")
            fp.write(f"# Discretization: N_lat={metadata.get('N_lat', 0)}, N_top={metadata.get('N_cap_top', 0)}, N_bot={metadata.get('N_cap_bot', 0)}\n")
            fp.write(f"# Ixx/Iyy={metadata.get('Ixx_over_Iyy', 1.0):.6f}, COM_norm={metadata.get('com_offset_norm', 0.0):.2e}\n")
        fp.write(f"{N}  {spacing:.16e}\n")
        for pt in coords:
            fp.write(f"{pt[0]:.16e}\t{pt[1]:.16e}\t{pt[2]:.16e}\n")
    print(f"Saved cylinder vertex file: {filepath} (N={N}, a_blob={a_blob:.6f})")

def generate_mtp_cylinder_hierarchy(output_dir, R_geom=1.0, L_geom=2.0):
    """Generates the required MTP refinement hierarchy: N in {12, 42, 162, 642, 2562}."""
    RESOLUTIONS = [12, 42, 162, 642, 2562]
    configs = {}
    
    for N in RESOLUTIONS:
        coords, a_blob, meta = generate_cylinder_multiblob(N, R_geom=R_geom, L_geom=L_geom)
        filename = f"cylinder_N_{N}_R_{R_geom:.1f}_L_{L_geom:.1f}.vertex"
        filepath = os.path.join(output_dir, filename)
        save_cylinder_vertex_file(filepath, coords, a_blob, meta)
        meta['vertex_file'] = filepath
        configs[N] = meta
        
    return configs

def plot_cylinder_geometry_validation(configs, output_path):
    """Generates a 5-panel figure validating cylinder geometries across resolutions."""
    fig = plt.figure(figsize=(18, 4), dpi=300)
    resolutions = [12, 42, 162, 642, 2562]
    
    for idx, N in enumerate(resolutions, 1):
        ax = fig.add_subplot(1, 5, idx, projection='3d')
        coords, a_blob, meta = generate_cylinder_multiblob(N, R_geom=1.0, L_geom=2.0)
        
        # Color by section: lateral (blue), top cap (red), bottom cap (green)
        N_lat = meta['N_lat']
        N_cap = meta['N_cap_top']
        
        ax.scatter(coords[:N_lat, 0], coords[:N_lat, 1], coords[:N_lat, 2],
                   c='#1f77b4', s=max(400.0/np.sqrt(N), 4), alpha=0.8, label='Barrel')
        ax.scatter(coords[N_lat:N_lat+N_cap, 0], coords[N_lat:N_lat+N_cap, 1], coords[N_lat:N_lat+N_cap, 2],
                   c='#d62728', s=max(400.0/np.sqrt(N), 4), alpha=0.8, label='Top Cap')
        ax.scatter(coords[N_lat+N_cap:, 0], coords[N_lat+N_cap:, 1], coords[N_lat+N_cap:, 2],
                   c='#2ca02c', s=max(400.0/np.sqrt(N), 4), alpha=0.8, label='Bottom Cap')
        
        ax.set_title(f"N = {N}\n$a_{{\\rm blob}}={a_blob:.3f}$, $R_g={meta['R_g']:.3f}$", fontsize=10, fontweight='bold')
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        ax.set_xlabel('X', fontsize=8)
        ax.set_ylabel('Y', fontsize=8)
        ax.set_zlabel('Z', fontsize=8)
        ax.tick_params(labelsize=6)
        if idx == 1:
            ax.legend(loc='upper right', fontsize=7)
            
    plt.suptitle("Rigid Circular Cylinder Multiblob Surface Discretization Hierarchy ($R=1.0, L=2.0$)",
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved geometry validation plot: {output_path}")

if __name__ == '__main__':
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'phase2_shapes', 'structures'))
    plots_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'phase2_shapes', 'plots'))
    configs = generate_mtp_cylinder_hierarchy(base_dir, R_geom=1.0, L_geom=2.0)
    
    print("\n" + "=" * 105)
    print("MTP RIGID CYLINDER MULTIBLOB REFINEMENT HIERARCHY SUMMARY")
    print("=" * 105)
    print(f"{'N':>6} | {'N_lat':>6} | {'N_cap':>6} | {'a_blob':>9} | {'min_d':>8} | {'R_g':>8} | {'R_g Err%':>9} | {'Ixx/Iyy':>9} | {'Ixy/Ixx':>9} | {'COM Norm':>10}")
    print("-" * 105)
    for N, c in configs.items():
        print(f"{N:6d} | {c['N_lat']:6d} | {c['N_cap_top']:6d} | {c['blob_radius']:9.5f} | {c['min_distance']:8.5f} | "
              f"{c['R_g']:8.5f} | {c['R_g_error_pct']:8.4f}% | {c['Ixx_over_Iyy']:9.5f} | {c['Ixy_over_Ixx']:9.2e} | {c['com_offset_norm']:10.2e}")
    print("=" * 105)
    
    plot_path = os.path.join(plots_dir, 'cylinder_geometry_validation.png')
    plot_cylinder_geometry_validation(configs, plot_path)
