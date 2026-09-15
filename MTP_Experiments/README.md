# MTP Experiments Workspace

This directory contains all experiments, data, plots, logs, and reports for the Major Thesis Project (MTP) on particle sedimentation in low-Reynolds-number flows.

## Strict Separation Architecture

* `../RigidMultiblobsWall/`: Reference hydrodynamic computing framework (read-only, untouched).
* `MTP_Experiments/`: All experiment drivers, generated CSVs, plots, logs, and reports.

```
MTP_Experiments/
├── scripts/                      # Reusable experiment drivers and common wrappers
│   ├── common.py
│   ├── phase1a_baseline.py
│   └── ...
├── phase1_sphere/                # Phase 1: Single Sphere Sedimentation (Unbounded)
│   ├── baseline/                 # Phase 1A baseline outputs & tolerance checks
│   ├── force_sweep/              # Phase 1B force sweep results & fits
│   ├── resolution_convergence/   # Phase 1C resolution study & convergence rates
│   ├── raw_data/                 # Machine-readable raw CSV outputs
│   ├── processed_data/           # Aggregated & normalized tables
│   ├── plots/                    # High-resolution publication plots
│   ├── reports/                  # Markdown validation reports
│   └── logs/                     # Console run logs and solver timing
├── phase2_shapes/                # Phase 2: Shape Validation (Disc, Cylinder, Ellipsoid, Boomerang)
│   ├── raw_data/
│   ├── processed_data/
│   ├── plots/
│   ├── reports/
│   └── logs/
└── README.md
```

All generated datasets include complete metadata:
`shape`, `N_blobs`, `blob_radius`, `R_g`, `R_h`, `force`, `viscosity`, `orientation`, `simulation velocity`, `angular velocity`, `theoretical velocity`, `error`, and `runtime`.
