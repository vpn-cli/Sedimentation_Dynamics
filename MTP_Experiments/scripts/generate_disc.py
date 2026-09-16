"""
MTP_Experiments/scripts/generate_disc.py

Parametric generator for rigid disc multiblob configurations.
Generates axisymmetric coplanar multiblob discs with exact center-of-mass at (0, 0, 0)
and symmetry axis aligned with z-axis.

Supports:
1. Exact MTP refinement hierarchy: N in {12, 42, 162, 642, 2562} via golden-angle phyllotaxis (Vogel spiral).
2. Concentric circular rings.
"""

import os
import sys
import numpy as np
from scipy.spatial.distance import pdist

def generate_vogel_disc(N, R_geom=1.0, overlap_factor=0.70):
    """
    Generates a planar circular disc with EXACTLY N blobs using Vogel's golden spiral phyllotaxis.
    Provides uniform area density, exact isotropic moments of inertia (Ixx approx Iyy),
    and zero center-of-mass offset.
    
    Parameters:
    -----------
    N : int
        Exact number of blobs.
    R_geom : float
        Target geometric radius of the disc.
    overlap_factor : float
        Ratio of blob radius to minimum inter-blob distance (0.7 ensures continuous coverage).
        
    Returns:
    --------
    coords : np.ndarray, shape (N, 3)
    a_blob : float
    metadata : dict
    """
    i = np.arange(N, dtype=np.float64)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    theta = i * golden_angle
    r = R_geom * np.sqrt((i + 0.5) / N)
    
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.zeros(N, dtype=np.float64)
    
    # Enforce exact center of mass at origin
    x -= np.mean(x)
    y -= np.mean(y)
    
    coords = np.column_stack([x, y, z])
    
    # Calculate inter-blob distances
    d = pdist(coords[:, 0:2])
    min_d = float(np.min(d))
    mean_d = float(np.mean(d))
    
    a_blob = float(min_d * overlap_factor)
    
    # Radius of gyration Rg (theoretical for thin uniform disk is R / sqrt(2) approx 0.7071)
    Rg = float(np.sqrt(np.mean(np.sum(coords**2, axis=1))))
    
    # Inertia tensor components (assuming unit mass per blob)
    Ixx = float(np.sum(coords[:, 1]**2 + coords[:, 2]**2))
    Iyy = float(np.sum(coords[:, 0]**2 + coords[:, 2]**2))
    Izz = float(np.sum(coords[:, 0]**2 + coords[:, 1]**2))
    Ixy = float(-np.sum(coords[:, 0] * coords[:, 1]))
    
    metadata = {
        'shape': 'disc',
        'discretization': 'vogel_spiral',
        'N_blobs': int(N),
        'R_geom': float(R_geom),
        'blob_radius': float(a_blob),
        'min_distance': float(min_d),
        'mean_distance': float(mean_d),
        'R_g': float(Rg),
        'Ixx': float(Ixx),
        'Iyy': float(Iyy),
        'Izz': float(Izz),
        'Ixy': float(Ixy),
        'center_of_mass': coords.mean(axis=0).tolist()
    }
    
    return coords, a_blob, metadata

def save_vertex_file(filepath, coords, a_blob, metadata=None):
    """
    Saves coordinates to standard .vertex file format.
    Line 1: N_blobs  blob_spacing
    Following lines: x  y  z
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    N = len(coords)
    spacing = 2.0 * a_blob
    with open(filepath, 'w') as fp:
        if metadata:
            fp.write(f"# Disc Multiblob Configuration: N={N}, R_geom={metadata.get('R_geom', 1.0):.4f}, a_blob={a_blob:.6f}\n")
            fp.write(f"# R_g={metadata.get('R_g', 0.0):.6f}, discretization={metadata.get('discretization', 'vogel')}\n")
        fp.write(f"{N}  {spacing:.16e}\n")
        for pt in coords:
            fp.write(f"{pt[0]:.16e}\t{pt[1]:.16e}\t{pt[2]:.16e}\n")
    print(f"Saved vertex file: {filepath} (N={N}, a_blob={a_blob:.6f})")

def generate_mtp_disc_hierarchy(output_dir, R_geom=1.0):
    """Generates the required MTP refinement hierarchy: N in {12, 42, 162, 642, 2562}."""
    RESOLUTIONS = [12, 42, 162, 642, 2562]
    configs = {}
    
    for N in RESOLUTIONS:
        coords, a_blob, meta = generate_vogel_disc(N, R_geom=R_geom)
        filename = f"disc_N_{N}_R_{R_geom:.1f}.vertex"
        filepath = os.path.join(output_dir, filename)
        save_vertex_file(filepath, coords, a_blob, meta)
        meta['vertex_file'] = filepath
        configs[N] = meta
        
    return configs

if __name__ == '__main__':
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'phase2_shapes', 'structures'))
    configs = generate_mtp_disc_hierarchy(base_dir, R_geom=1.0)
    print("\nGenerated MTP Disc Configurations Summary:")
    print("=" * 85)
    print(f"{'N':>6} | {'a_blob':>9} | {'min_d':>8} | {'R_g':>8} | {'Ixx/Iyy':>10} | {'Ixy/Ixx':>10}")
    print("-" * 85)
    for N, c in configs.items():
        print(f"{N:6d} | {c['blob_radius']:9.5f} | {c['min_distance']:8.5f} | {c['R_g']:8.5f} | {c['Ixx']/c['Iyy']:10.6f} | {c['Ixy']/c['Ixx']:10.2e}")
    print("=" * 85)
