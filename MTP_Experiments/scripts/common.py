"""
MTP_Experiments/scripts/common.py

Shared utilities, metadata tracking, and solver wrappers for MTP experiments.
Imports directly from ../RigidMultiblobsWall without modifying any original code.
"""

import sys
import os
import types
import importlib.util
import time

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

# Path resolution: configurable and relative to this script
EXPERIMENTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WORKSPACE_DIR = os.path.abspath(os.path.join(EXPERIMENTS_DIR, '..'))
REPO_PATH = os.path.abspath(os.path.join(WORKSPACE_DIR, 'RigidMultiblobsWall'))

if REPO_PATH not in sys.path:
    sys.path.insert(0, REPO_PATH)

import numpy as np
import scipy.linalg
import scipy.sparse.linalg as spla
from body.body import Body
from quaternion_integrator.quaternion import Quaternion
from mobility import mobility as mb
from mobility import mobility_numba
from read_input.read_vertex_file import read_vertex_file

# Directory configuration for Phase 1
PHASE1_DIR = os.path.join(EXPERIMENTS_DIR, 'phase1_sphere')
PHASE1_BASELINE_DIR = os.path.join(PHASE1_DIR, 'baseline')
PHASE1_FORCE_DIR = os.path.join(PHASE1_DIR, 'force_sweep')
PHASE1_RES_DIR = os.path.join(PHASE1_DIR, 'resolution_convergence')
PHASE1_RAW_DATA_DIR = os.path.join(PHASE1_DIR, 'raw_data')
PHASE1_PROCESSED_DATA_DIR = os.path.join(PHASE1_DIR, 'processed_data')
PHASE1_PLOTS_DIR = os.path.join(PHASE1_DIR, 'plots')
PHASE1_REPORTS_DIR = os.path.join(PHASE1_DIR, 'reports')
PHASE1_LOGS_DIR = os.path.join(PHASE1_DIR, 'logs')

STRUCTURES_DIR = os.path.join(REPO_PATH, 'multi_bodies', 'Structures')

# Verified sphere metadata from repository audit
SPHERE_CONFIGS = {
    'Rh_calibrated': {
        12: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_12_Rg_0_7921_Rh_1.vertex'),
            'Rg': 0.7921,
            'Rh': 1.0000,
            'blob_radius': 0.41642068286675,
            'Dx': 0.8328413657335
        },
        42: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_42_Rg_0_8913_Rh_1.vertex'),
            'Rg': 0.8913,
            'Rh': 1.0000,
            'blob_radius': 0.243553056072,
            'Dx': 0.487106112144
        },
        162: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_162_Rg_0_9497_Rh_1.vertex'),
            'Rg': 0.9497,
            'Rh': 1.0000,
            'blob_radius': 0.13100877695,
            'Dx': 0.2620175539
        },
        642: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_642_Rg_0_9767_Rh_1.vertex'),
            'Rg': 0.9767,
            'Rh': 1.0000,
            'blob_radius': 0.067527675333,
            'Dx': 0.135055350666
        },
        2562: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_2562_Rg_0_9888_Rh_1.vertex'),
            'Rg': 0.9888,
            'Rh': 1.0000,
            'blob_radius': 0.034204978919,
            'Dx': 0.068409957838
        }
    },
    'Rg_geometric': {
        12: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_12_Rg_1_Rh_1_2625.vertex'),
            'Rg': 1.0000,
            'Rh': 1.2625,
            'blob_radius': 0.525731112119,
            'Dx': 1.051462224238
        },
        42: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_42_Rg_1_Rh_1_1220.vertex'),
            'Rg': 1.0000,
            'Rh': 1.1220,
            'blob_radius': 0.273266528913,
            'Dx': 0.546533057825
        },
        162: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_162_Rg_1_Rh_1_0530.vertex'),
            'Rg': 1.0000,
            'Rh': 1.0530,
            'blob_radius': 0.137952242128,
            'Dx': 0.275904484255
        },
        642: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_642_Rg_1_Rh_1_0239.vertex'),
            'Rg': 1.0000,
            'Rh': 1.0239,
            'blob_radius': 0.069141586774,
            'Dx': 0.138283173547
        },
        2562: {
            'vertex_file': os.path.join(STRUCTURES_DIR, 'shell_N_2562_Rg_1_Rh_1_0113.vertex'),
            'Rg': 1.0000,
            'Rh': 1.0113,
            'blob_radius': 0.034591495181,
            'Dx': 0.069182990361
        }
    }
}

def load_sphere_structure(sphere_type, N):
    """Loads multiblob coordinates and metadata from vertex file."""
    meta = SPHERE_CONFIGS[sphere_type][N]
    r_conf = read_vertex_file(meta['vertex_file'])
    return r_conf, meta

def solve_sphere_unbounded(r_conf, a_blob, eta, force_vec, torque_vec=None, method='auto'):
    """
    Solves unbounded rigid multiblob Stokes sedimentation.
    Uses Cholesky factorization for N <= 642, and preconditioned GMRES with native
    Numba parallel matrix-vector product for N = 2562 (avoiding dense inversion).
    Returns (U, Omega, runtime_seconds).
    """
    t0 = time.time()
    Nblobs = len(r_conf)
    force = np.asarray(force_vec, dtype=np.float64).flatten()
    torque = np.zeros(3) if torque_vec is None else np.asarray(torque_vec, dtype=np.float64).flatten()
    wrench = np.concatenate([force, torque])
    
    q0 = Quaternion([1.0, 0.0, 0.0, 0.0])
    b = Body(np.zeros(3), q0, r_conf, a_blob)
    K = b.calc_K_matrix()
    
    if method == 'auto':
        method = 'gmres' if Nblobs >= 2562 else 'cholesky'
        
    if method == 'cholesky':
        M = mb.rotne_prager_tensor(r_conf, eta, a_blob)
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
            Mf = mobility_numba.no_wall_mobility_trans_times_force_numba(r_conf, f_blob, eta, a_blob, L_box).flatten()
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

def stokes_velocity(F_mag, eta, R):
    """Analytical Stokes velocity: U = F / (6 * pi * eta * R)"""
    return F_mag / (6.0 * np.pi * eta * R)

def angle_between_vectors(v1, v2):
    """Angle in degrees between two 3D vectors."""
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 0.0
    cos_val = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return np.degrees(np.arccos(cos_val))
