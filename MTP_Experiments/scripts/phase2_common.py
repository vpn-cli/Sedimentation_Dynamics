"""
MTP_Experiments/scripts/phase2_common.py

Phase 2: Common utilities, geometry generators, analytical solutions,
and hydrodynamic solver routines for non-spherical rigid bodies (Ellipsoid & Robotic Arm).
Strictly interfaces with ../RigidMultiblobsWall without modifying repository files.
"""

import sys
import os
import types
import importlib.util
import time
import numpy as np
import scipy.linalg
import scipy.sparse.linalg as spla

# Python 3.12+ compatibility shim for legacy 'imp' imports in RigidMultiblobsWall
if 'imp' not in sys.modules:
    imp = types.ModuleType('imp')
    def find_module(name):
        spec = importlib.util.find_spec(name)
        if spec is None:
            raise ImportError(f"No module named '{name}'")
        return spec
    imp.find_module = find_module
    sys.modules['imp'] = imp

# Path resolution: dynamic relative paths
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
EXPERIMENTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
WORKSPACE_DIR = os.path.abspath(os.path.join(EXPERIMENTS_DIR, '..'))
REPO_PATH = os.path.abspath(os.path.join(WORKSPACE_DIR, 'RigidMultiblobsWall'))

if REPO_PATH not in sys.path:
    sys.path.insert(0, REPO_PATH)

from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from mobility import mobility_numba
from read_input.read_vertex_file import read_vertex_file

# Import baseline SPHERE_CONFIGS from common.py
from common import SPHERE_CONFIGS, STRUCTURES_DIR

# Phase 2 Directory Structure
PHASE2_DIR = os.path.join(EXPERIMENTS_DIR, 'phase2_nonspherical')
ELLIPSOID_DIR = os.path.join(PHASE2_DIR, 'ellipsoid')
ROBOTIC_ARM_DIR = os.path.join(PHASE2_DIR, 'robotic_arm')
PHASE2_PROCESSED_DATA_DIR = os.path.join(PHASE2_DIR, 'processed_data')
PHASE2_PLOTS_DIR = os.path.join(PHASE2_DIR, 'plots')
PHASE2_REPORTS_DIR = os.path.join(PHASE2_DIR, 'reports')
PHASE2_LOGS_DIR = os.path.join(PHASE2_DIR, 'logs')


# =============================================================================
# 1. ANALYTICAL SOLUTIONS: PERRIN FORMULAS FOR PROLATE SPHEROID
# =============================================================================

def perrin_analytical_prolate(a, b, eta=1.0):
    """
    Computes exact analytical Perrin (1934, 1936) drag and mobility coefficients
    for a prolate spheroid with semi-major axis 'a' (along symmetry axis) and
    semi-minor axis 'b' (b = c, aspect ratio lambda = a/b > 1) in unbounded Stokes flow.

    Returns dictionary containing:
      - 'aspect_ratio': a/b
      - 'eccentricity': e = sqrt(1 - b^2/a^2)
      - 'Req': (a * b^2)**(1/3) equivalent volume sphere radius
      - 'K_parallel': shape factor for motion along symmetry axis
      - 'K_perp': shape factor for motion perpendicular to symmetry axis
      - 'mu_parallel': translational mobility parallel to symmetry axis
      - 'mu_perp': translational mobility perpendicular to symmetry axis
      - 'anisotropy_ratio': mu_parallel / mu_perp
      - 'C_parallel': rotational shape factor parallel to symmetry axis
      - 'C_perp': rotational shape factor perpendicular to symmetry axis
      - 'mu_r_parallel': rotational mobility about symmetry axis
      - 'mu_r_perp': rotational mobility about perpendicular axis
    """
    if a < b:
        raise ValueError("For prolate spheroid, semi-major axis 'a' must be >= 'b'")
    
    if np.isclose(a, b, atol=1e-12):
        # Spherical limit
        Req = a
        mu_trans = 1.0 / (6.0 * np.pi * eta * a)
        mu_rot = 1.0 / (8.0 * np.pi * eta * (a**3))
        return {
            'aspect_ratio': 1.0,
            'eccentricity': 0.0,
            'Req': Req,
            'K_parallel': 1.0,
            'K_perp': 1.0,
            'mu_parallel': mu_trans,
            'mu_perp': mu_trans,
            'anisotropy_ratio': 1.0,
            'C_parallel': 1.0,
            'C_perp': 1.0,
            'mu_r_parallel': mu_rot,
            'mu_r_perp': mu_rot
        }

    lam = a / b
    e = np.sqrt(1.0 - (b / a)**2)
    L = np.log((1.0 + e) / (1.0 - e))
    Req = (a * (b**2))**(1.0 / 3.0)

    # Translational shape factors (Perrin 1934; Happel & Brenner 1965, p. 147)
    denom_par = -2.0 * e + (1.0 + e**2) * L
    K_par = (8.0 / 3.0) * (e**3) / denom_par

    denom_perp = 2.0 * e + (3.0 * e**2 - 1.0) * L
    K_perp = (16.0 / 3.0) * (e**3) / denom_perp

    mu_par = 1.0 / (6.0 * np.pi * eta * a * K_par)
    mu_perp = 1.0 / (6.0 * np.pi * eta * a * K_perp)

    # Rotational shape factors (Perrin 1934; Kim & Karrila 1991)
    denom_r_par = 2.0 * e - (1.0 - e**2) * L
    C_par = (2.0 / 3.0) * (e**3) * (1.0 - e**2) / denom_r_par

    denom_r_perp = -2.0 * e + (1.0 + e**2) * L
    C_perp = (4.0 / 3.0) * (e**3) * (2.0 - e**2) / denom_r_perp

    mu_r_par = 1.0 / (8.0 * np.pi * eta * a * (b**2) * C_par)
    mu_r_perp = 1.0 / (8.0 * np.pi * eta * a * (b**2) * C_perp)

    return {
        'aspect_ratio': lam,
        'eccentricity': e,
        'Req': Req,
        'K_parallel': K_par,
        'K_perp': K_perp,
        'mu_parallel': mu_par,
        'mu_perp': mu_perp,
        'anisotropy_ratio': mu_par / mu_perp,
        'C_parallel': C_par,
        'C_perp': C_perp,
        'mu_r_parallel': mu_r_par,
        'mu_r_perp': mu_r_perp
    }


# =============================================================================
# 2. GEOMETRY GENERATION: ELLIPSOID MESH
# =============================================================================

def generate_ellipsoid_mesh(N, a, b, c=None, sphere_type='Rh_calibrated', custom_blob_radius=None):
    """
    Generates an ellipsoidal shell multiblob configuration by affine scaling of
    the validated unit spherical shell triangulations from the repository.

    Parameters:
      - N: Blob resolution (12, 42, 162, 642, 2562)
      - a: Semi-major axis (along x-axis)
      - b: Semi-minor axis (along y-axis)
      - c: Semi-minor axis (along z-axis, defaults to b for prolate spheroid)
      - sphere_type: 'Rh_calibrated' (default) or 'Rg_geometric'
      - custom_blob_radius: Optional override for blob hydrodynamic radius

    Returns:
      - r_conf: (N, 3) numpy array of blob coordinates in body frame, centered at (0, 0, 0)
      - a_blob: Blob hydrodynamic radius
      - meta: Metadata dictionary
    """
    if c is None:
        c = b
    
    if N not in SPHERE_CONFIGS[sphere_type]:
        raise ValueError(f"Resolution N={N} not in SPHERE_CONFIGS[{sphere_type}]")
    
    sphere_meta = SPHERE_CONFIGS[sphere_type][N]
    vfile = sphere_meta['vertex_file']
    r_sphere = read_vertex_file(vfile)

    # Unitize sphere vertices if needed (so geometric radius = 1.0)
    Rg_sphere = sphere_meta['Rg']
    r_unit = r_sphere / Rg_sphere

    # Affine scaling: stretch along principal axes
    r_ellip = np.zeros_like(r_unit)
    r_ellip[:, 0] = a * r_unit[:, 0]
    r_ellip[:, 1] = b * r_unit[:, 1]
    r_ellip[:, 2] = c * r_unit[:, 2]

    # Effective volume radius
    Req = (a * b * c)**(1.0 / 3.0)

    # Scale blob radius proportionally to equivalent radius
    if custom_blob_radius is not None:
        a_blob = float(custom_blob_radius)
    else:
        # Scale blob radius with equivalent radius Req
        a_blob = sphere_meta['blob_radius'] * (Req / Rg_sphere)

    meta = {
        'N': N,
        'a': a,
        'b': b,
        'c': c,
        'Req': Req,
        'aspect_ratio': a / b,
        'blob_radius': a_blob,
        'sphere_type': sphere_type,
        'source_vertex_file': vfile
    }

    return r_ellip, a_blob, meta


# =============================================================================
# 3. GEOMETRY GENERATION: ROBOTIC ARM / COMPLEX RIGID BODY
# =============================================================================

def assemble_robotic_arm_rigid(N_links=7, link_resolution=12, link_spacing=2.5,
                                link_radius=1.0, config='straight', bend_angle_deg=0.0,
                                center_at='centroid'):
    """
    Assembles a multi-segment robotic arm into a SINGLE RIGID BODY.
    Uses spherical multiblob shells for each segment as specified in the repository's
    'robot_arm_N_7.clones' and 'Generate_list_vertex_clones_const_files.m'.

    Parameters:
      - N_links: Number of spherical segments in chain (e.g. 7 or 15)
      - link_resolution: 1 (single blob) or 12 (shell_N_12) or 42 (shell_N_42)
      - link_spacing: Center-to-center distance between adjacent links (default 2.5)
      - link_radius: Geometric radius of each spherical link (default 1.0)
      - config: 'straight' or 'bent' (elbow bend at middle link)
      - bend_angle_deg: Bend angle in degrees for 'bent' config (e.g. 45 or 90)
      - center_at: 'centroid' (centers body frame at geometric centroid) or 'root' (first link at 0)

    Returns:
      - r_conf: (N_total, 3) coordinates of all blobs in the rigid body frame
      - a_blob: Hydrodynamic radius of individual blobs
      - meta: Metadata dictionary
    """
    if link_resolution == 1:
        r_link = np.zeros((1, 3))
        a_blob = link_radius
    elif link_resolution == 12:
        vfile = os.path.join(STRUCTURES_DIR, 'shell_N_12_Rg_0_7921_Rh_1.vertex')
        r_link = read_vertex_file(vfile)
        a_blob = SPHERE_CONFIGS['Rh_calibrated'][12]['blob_radius'] * link_radius
    elif link_resolution == 42:
        vfile = os.path.join(STRUCTURES_DIR, 'shell_N_42_Rg_0_8913_Rh_1.vertex')
        r_link = read_vertex_file(vfile)
        a_blob = SPHERE_CONFIGS['Rh_calibrated'][42]['blob_radius'] * link_radius
    else:
        raise ValueError(f"Unsupported link_resolution: {link_resolution}. Use 1, 12, or 42.")

    # Determine link center locations
    link_centers = np.zeros((N_links, 3))
    mid_index = N_links // 2

    cur_pos = np.zeros(3)
    cur_dir = np.array([1.0, 0.0, 0.0])

    for i in range(N_links):
        if i == 0:
            link_centers[i] = cur_pos
        else:
            if config == 'bent' and i == mid_index:
                theta = np.radians(bend_angle_deg)
                # Rotate direction vector about z-axis
                rot_z = np.array([
                    [np.cos(theta), -np.sin(theta), 0.0],
                    [np.sin(theta),  np.cos(theta), 0.0],
                    [0.0,            0.0,           1.0]
                ])
                cur_dir = np.dot(rot_z, cur_dir)
            cur_pos = cur_pos + link_spacing * cur_dir
            link_centers[i] = cur_pos

    # Center coordinates according to center_at
    centroid = np.mean(link_centers, axis=0)
    if center_at == 'centroid':
        link_centers -= centroid
    elif center_at == 'root':
        pass  # Root link already at (0, 0, 0)
    else:
        raise ValueError(f"Unknown center_at option: {center_at}")

    # Assemble all blobs into a unified rigid reference configuration
    all_blobs = []
    for center in link_centers:
        blobs_in_link = r_link + center
        all_blobs.append(blobs_in_link)
    
    r_conf = np.concatenate(all_blobs, axis=0)

    meta = {
        'N_links': N_links,
        'link_resolution': link_resolution,
        'N_blobs_total': len(r_conf),
        'link_spacing': link_spacing,
        'link_radius': link_radius,
        'blob_radius': a_blob,
        'config': config,
        'bend_angle_deg': bend_angle_deg,
        'center_at': center_at,
        'centroid': centroid
    }

    return r_conf, a_blob, meta


# =============================================================================
# 4. SOLVER: FULL 6x6 MOBILITY TENSOR & SYMMETRIES
# =============================================================================

def solve_rigid_mobility_tensor(r_conf, a_blob, eta=1.0, location=None, orientation=None):
    """
    Computes the full 6x6 rigid body mobility tensor:
      N_body = (K^T * M^-1 * K)^-1 = [ [M_tt, M_tr],
                                       [M_rt, M_rr] ]
    mapping generalized wrench [F; T] to kinematic velocity [U; Omega].

    Also performs rigorous verification of:
      1. Onsager reciprocal symmetry: ||M_tr - M_rt^T||_inf
      2. Symmetric positive-definiteness: all 6 eigenvalues > 0
    """
    if location is None:
        location = np.zeros(3)
    if orientation is None:
        orientation = Quaternion([1.0, 0.0, 0.0, 0.0])

    b = Body(location, orientation, r_conf, a_blob)
    K = b.calc_K_matrix()  # (3*N, 6)
    r_lab = b.get_r_vectors()

    # Compute RPY blob-blob mobility matrix using lab-frame coordinates
    M = mb.rotne_prager_tensor(r_lab, eta, a_blob)

    # Cholesky factorization of M for numerical stability
    L, lower = scipy.linalg.cho_factor(M)
    Minv_K = scipy.linalg.cho_solve((L, lower), K, check_finite=False)
    Kt_Minv_K = np.dot(K.T, Minv_K)

    # Invert (K^T M^-1 K) to obtain 6x6 body mobility
    N_body = np.linalg.pinv(Kt_Minv_K)

    M_tt = N_body[0:3, 0:3]
    M_tr = N_body[0:3, 3:6]
    M_rt = N_body[3:6, 0:3]
    M_rr = N_body[3:6, 3:6]

    # Symmetry checks
    sym_error_tt = np.max(np.abs(M_tt - M_tt.T))
    sym_error_rr = np.max(np.abs(M_rr - M_rr.T))
    onsager_error = np.max(np.abs(M_tr - M_rt.T))
    eigenvalues = np.linalg.eigvalsh(N_body)

    metrics = {
        'N_body': N_body,
        'M_tt': M_tt,
        'M_tr': M_tr,
        'M_rt': M_rt,
        'M_rr': M_rr,
        'sym_error_tt': sym_error_tt,
        'sym_error_rr': sym_error_rr,
        'onsager_error': onsager_error,
        'eigenvalues': eigenvalues,
        'is_spd': bool(np.all(eigenvalues > 0)),
        'min_eigenvalue': np.min(eigenvalues)
    }

    return N_body, metrics


# =============================================================================
# 5. SOLVER: STOKES SEDIMENTATION VELOCITIES [U, OMEGA]
# =============================================================================

def solve_sedimentation(r_conf, a_blob, eta, force_vec, torque_vec=None, orientation=None, method='auto'):
    """
    Solves unbounded rigid multiblob Stokes sedimentation.
    Uses Cholesky factorization for N <= 642, and preconditioned GMRES for N >= 2562.
    Returns (U, Omega, runtime_seconds).
    """
    t0 = time.time()
    Nblobs = len(r_conf)
    force = np.asarray(force_vec, dtype=np.float64).flatten()
    torque = np.zeros(3) if torque_vec is None else np.asarray(torque_vec, dtype=np.float64).flatten()
    wrench = np.concatenate([force, torque])

    if orientation is None:
        orientation = Quaternion([1.0, 0.0, 0.0, 0.0])

    b = Body(np.zeros(3), orientation, r_conf, a_blob)
    K = b.calc_K_matrix()
    r_lab = b.get_r_vectors()

    if method == 'auto':
        method = 'gmres' if Nblobs >= 2562 else 'cholesky'

    if method == 'cholesky':
        M = mb.rotne_prager_tensor(r_lab, eta, a_blob)
        L, lower = scipy.linalg.cho_factor(M)
        Kt_Minv_K = np.dot(K.T, scipy.linalg.cho_solve((L, lower), K, check_finite=False))
        N_body = np.linalg.pinv(Kt_Minv_K)
        U_Omega = np.dot(N_body, wrench)
        U = U_Omega[0:3]
        Omega = U_Omega[3:6]
    elif method == 'gmres':
        L_box = np.zeros(3)
        def matvec(v):
            f_blob = v[:3*Nblobs]
            u_body = v[3*Nblobs:3*Nblobs+6]
            Mf = mobility_numba.no_wall_mobility_trans_times_force_numba(r_lab, f_blob, eta, a_blob, L_box).flatten()
            Ku = np.dot(K, u_body)
            res_blob = Mf - Ku
            res_body = -np.dot(K.T, f_blob)
            return np.concatenate([res_blob, res_body])

        sys_size = 3*Nblobs + 6
        A_op = spla.LinearOperator((sys_size, sys_size), matvec=matvec, dtype=np.float64)

        m_self_inv = 6.0 * np.pi * eta * a_blob
        Kt_M0inv_K = m_self_inv * np.dot(K.T, K)
        S_inv = np.linalg.inv(Kt_M0inv_K)

        def pc_matvec(r):
            r_blob = r[:3*Nblobs]
            r_body = r[3*Nblobs:3*Nblobs+6]
            lambda_0 = m_self_inv * r_blob
            rhs_s = -r_body - np.dot(K.T, lambda_0)
            u_approx = np.dot(S_inv, rhs_s)
            lam_approx = m_self_inv * (r_blob + np.dot(K, u_approx))
            return np.concatenate([lam_approx, u_approx])

        PC_op = spla.LinearOperator((sys_size, sys_size), matvec=pc_matvec, dtype=np.float64)
        RHS = np.concatenate([np.zeros(3*Nblobs), -wrench])

        sol, info = spla.gmres(A_op, RHS, M=PC_op, rtol=1e-8, atol=1e-10, restart=60, maxiter=300)
        if info != 0:
            raise RuntimeError(f"GMRES did not converge (info={info})")

        U = sol[3*Nblobs:3*Nblobs+3]
        Omega = sol[3*Nblobs+3:3*Nblobs+6]
    else:
        raise ValueError(f"Unknown solver method: {method}")

    runtime = time.time() - t0
    return U, Omega, runtime


# =============================================================================
# 6. CENTER OF MOBILITY (CoM) EVALUATION
# =============================================================================

def compute_center_of_mobility(r_conf, a_blob, eta=1.0, initial_reference_point=None):
    """
    Computes the Center of Mobility (CoM / CoH) of an arbitrary rigid body
    using Delong's equation (as implemented in boomerang/evaluate_com_mobility.py).
    
    The condition is that at the CoM, M_wF is symmetric: (M_wF)^T = M_wF.
    Solves for vector r_shift from initial_reference_point to CoM.
    """
    if initial_reference_point is None:
        initial_reference_point = np.zeros(3)

    eijk = np.zeros((3, 3, 3))
    eijk[0, 1, 2] = 1;  eijk[1, 2, 0] = 1;  eijk[2, 0, 1] = 1
    eijk[1, 0, 2] = -1; eijk[0, 2, 1] = -1; eijk[2, 1, 0] = -1

    N_body, _ = solve_rigid_mobility_tensor(r_conf, a_blob, eta=eta, location=initial_reference_point)
    M_wF = N_body[3:6, 0:3]  # M_rt
    M_wT = N_body[3:6, 3:6]  # M_rr

    A = M_wF - M_wF.T

    B = np.zeros((3, 3, 3))
    for l in range(3):
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    B[l, i, j] += eijk[i, k, l] * M_wT[j, k] - eijk[j, k, l] * M_wT[i, k]

    RHS = np.array([A[0, 1], A[0, 2], A[1, 2]])
    C = np.zeros((3, 3))
    for i in range(3):
        C[0, i] = B[i, 0, 1]
        C[1, i] = B[i, 0, 2]
        C[2, i] = B[i, 1, 2]

    try:
        r_shift = np.linalg.solve(C, RHS)
    except np.linalg.LinAlgError:
        r_shift = np.linalg.lstsq(C, RHS, rcond=None)[0]

    com_location = initial_reference_point + r_shift

    # Verify at new location
    r_conf_com = r_conf - r_shift
    N_body_com, _ = solve_rigid_mobility_tensor(r_conf_com, a_blob, eta=eta, location=np.zeros(3))
    M_wF_com = N_body_com[3:6, 0:3]
    com_asymmetry = np.max(np.abs(M_wF_com - M_wF_com.T))

    return {
        'initial_reference_point': initial_reference_point,
        'r_shift_to_com': r_shift,
        'com_location': com_location,
        'com_asymmetry': com_asymmetry,
        'N_body_at_com': N_body_com
    }
