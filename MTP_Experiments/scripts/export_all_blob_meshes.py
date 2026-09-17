import sys
import os
import json

SCRIPT_DIR = r"c:\MTP\Sedimentation_Dynamics\MTP_Experiments\scripts"
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import phase2_common as p2c

resolutions = [12, 42, 162, 642, 2562]
a = 2.0
b = 1.0

# Numerical mobility data from resolution convergence study
mobility_data = {
    12: {'mu_par': 0.03524, 'mu_perp': 0.03104, 'err_par': 20.03, 'err_perp': 19.31},
    42: {'mu_par': 0.03897, 'mu_perp': 0.03438, 'err_par': 11.56, 'err_perp': 10.63},
    162: {'mu_par': 0.04161, 'mu_perp': 0.03660, 'err_par': 5.56, 'err_perp': 4.87},
    642: {'mu_par': 0.04293, 'mu_perp': 0.03762, 'err_par': 2.57, 'err_perp': 2.21},
    2562: {'mu_par': 0.04353, 'mu_perp': 0.03807, 'err_par': 1.21, 'err_perp': 1.04}
}

meshes = {}
for N in resolutions:
    r_conf, a_blob, meta = p2c.generate_ellipsoid_mesh(N, a, b, sphere_type='Rh_calibrated')
    # Round coordinates to 5 decimals for compact JS file size
    r_list = [[round(coord, 5) for coord in pt] for pt in r_conf.tolist()]
    meshes[str(N)] = {
        'N': N,
        'a_blob': round(float(a_blob), 5),
        'mu_par': mobility_data[N]['mu_par'],
        'mu_perp': mobility_data[N]['mu_perp'],
        'err_par': mobility_data[N]['err_par'],
        'err_perp': mobility_data[N]['err_perp'],
        'r_conf': r_list
    }
    print(f"N={N}: {len(r_list)} blobs, a_blob={a_blob:.5f}")

out_path = r"c:\MTP\Sedimentation_Dynamics\MTP_Experiments\phase2_nonspherical\ellipsoid\processed_data\all_blob_meshes.json"
with open(out_path, 'w') as f:
    json.dump(meshes, f)

print(f"Saved all meshes to {out_path}")
