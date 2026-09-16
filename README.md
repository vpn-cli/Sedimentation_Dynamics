# Sedimentation Dynamics in Low-Reynolds-Number Stokes Flows

Numerical investigation of particle shape effects on low-Reynolds-number sedimentation and sedimentation stability, utilizing the multiblob method (`RigidMultiblobsWall`).

---

## 🔬 Research Overview

In creeping flows ($\mathrm{Re} \ll 1$), fluid motion is governed by the linear Stokes equations:
$$\eta \nabla^2 \mathbf{u} - \nabla p = \mathbf{0}, \quad \nabla \cdot \mathbf{u} = 0$$

Under this regime:
1. **Particle Anisotropy & Oblique Drift**: Non-spherical particles (e.g. prolate spheroids, ellipsoids, discs) experience anisotropic translational mobility ($\mu_\parallel > \mu_\perp$). When inclined at an angle $\theta$, gravitational loading induces an oblique lateral drift velocity:
   $$U_x(\theta) = +\frac{1}{2} F_z (\mu_\parallel - \mu_\perp) \sin(2\theta)$$
   reaching maximum lateral drift at $\theta = 45^\circ$.
2. **Non-Tumbling Pose Stability**: Due to orthotropic reflection symmetry ($D_{2h}$), translation-rotation coupling vanishes at the centroid ($M_{tr} \equiv 0$). Isolated ellipsoids sediment without tumbling ($\mathbf{\Omega} \equiv \mathbf{0}$).
3. **Multi-Resolution Convergence**: Particle surfaces are discretized into multiblob meshes using geodesic icosahedral subdivision ($N \in \{12, 42, 162, 642\}$), converging monotonically toward analytical continuum benchmarks (Stokes and Perrin theory).

---

## 📐 Dimensionless Formulation & Scaling

All simulations are formulated in **characteristic dimensionless (reduced) Stokesian units**:
* **Characteristic Length ($L_c$)**: Semi-minor axis $b = 1.0$ (semi-major axis $a = 2.0$, aspect ratio $\lambda = 2.0$).
* **Fluid Viscosity ($\eta_c$)**: Dynamic viscosity $\eta = 1.0$.
* **Sedimentation Force ($F_c$)**: Net gravitational buoyant force magnitude $F_z = 1.0$.

### Mapping to Physical SI Units
Because Stokes flow is scale-invariant, dimensionless simulation values map directly to any physical experiment:
$$\mathbf{r}_{\text{physical}} = \mathbf{r}^* \cdot L_c, \qquad \mathbf{U}_{\text{physical}} = \mathbf{U}^* \cdot \left(\frac{F_c}{\eta_c L_c}\right), \qquad t_{\text{physical}} = t^* \cdot \left(\frac{\eta_c L_c^2}{F_c}\right)$$

---

## 📂 Repository Structure

```text
Sedimentation_Dynamics/
├── MTP_Experiments/
│   ├── phase1_sphere/             # Sphere benchmarks, force sweeps, and baseline calibration
│   ├── phase2_nonspherical/       # Non-spherical body hydrodynamics
│   │   ├── ellipsoid/             # Prolate ellipsoid benchmarks & 3D viewer
│   │   ├── robotic_arm/           # 7-segment articulated robotic arm & 3D viewer
│   │   └── boomerang/             # Bent 2-arm boomerang particle & 3D viewer
│   │       ├── boomerang_simulation_viewer.html  # Interactive 3D WebGL Multi-Resolution Viewer
│   │       ├── plots/                           # 300 DPI publication-ready figures & animated simulation GIF
│   │       ├── processed_data/                  # Discretized meshes & trajectory JSON data
│   │       └── reports/                         # Markdown validation reports
│   └── scripts/                   # Simulation, validation, and analysis scripts
└── RigidMultiblobsWall/           # Reference multiblob Stokesian hydrodynamics library
```

---

## 🪃 Boomerang Colloidal Particle (Phase 2 Focus)

The bent two-arm **boomerang particle** (`RigidMultiblobsWall/multi_bodies/Structures/boomerang_N_15.vertex`) investigates low-Reynolds-number sedimentation of non-axisymmetric, articulated particles with tunable opening angle $\alpha$ and chirality.

### Key Hydrodynamic Insights
1. **Center of Mobility (CoM) vs Centroid Audit**:
   - Discretized at $N=15$ blobs ($L=2.1, a_{\text{blob}}=0.25$).
   - Apex at $(0, 0, 0)$, geometric centroid at $(0.560, 0.560, 0.0)$, and hydrodynamic CoM at $(0.664, 0.664, 0.0)$ ($\Delta \mathbf{r} = +0.104$).
   - Tracking dynamics about the apex introduces spurious lever-arm torques ($\|M_{tr}\| = 0.1110$). Shifting to the CoM minimizes translation-rotation coupling to $\|M_{tr}\| = 0.0055$ ($>95\%$ reduction).
2. **Symmetry & Rigorous Validation**:
   - Onsager reciprocal symmetry error: $\|M_{tr} - M_{rt}^T\|_\infty = 8.67 \times 10^{-18}$ (machine precision).
   - Positive-definiteness: $\lambda_{\min} = 0.02918 > 0$.
   - Force linearity: Exact proportionality across all axes ($R^2 = 1.00000000$).
3. **Four Dynamic Sedimentation Regimes**:
   - **Apex Down (Edge-On Gliding)**: Stable orientation ($\mathbf{\Omega} \approx 0$). Slices vertically through fluid with minimal drag.
   - **Flat Pose (Horizontal)**: High drag settling ($|U_z| = 0.05914$), exhibiting pitching reorientation about the bisecting axis.
   - **Tilted 45° (Oblique Drift)**: Anisotropic lateral velocity ($U_x, U_y \neq 0$) coupled with gradual alignment.
   - **Chiral Spiral (15° Dihedral Twist)**: Breaks planar reflection symmetry ($C_s \to C_1$), coupling vertical sedimentation force $F_z$ directly to vertical hydrodynamic torque $T_z$, producing continuous autorotation and a 3D helical spiral trajectory.

---

## 🚀 Running Simulations & Interactive 3D Viewers

### 1. Boomerang Particle
```bash
# Run validation suite (Onsager, CoM, linearity, angle sweep)
python MTP_Experiments/scripts/phase2_boomerang_validation.py

# Run 6-DOF dynamic trajectory simulation
python MTP_Experiments/scripts/phase2_boomerang_trajectory_sim.py

# Render 300 DPI animated simulation GIF
python MTP_Experiments/scripts/phase2_boomerang_animated_sim.py

# Build interactive 3D WebGL Multi-Resolution Viewer
python MTP_Experiments/scripts/generate_boomerang_html.py
```
**Launch Boomerang 3D Viewer:**
```bash
python -m http.server 8086 --directory MTP_Experiments/phase2_nonspherical/boomerang
```
Open `http://localhost:8086/boomerang_simulation_viewer.html` to inspect real-time 3D settling, toggle resolutions ($N=7, 15, 29$), switch between the 4 sedimentation regimes, inspect velocity vectors, and view real-time telemetry HUD.

### 2. Ellipsoid Benchmark
```bash
python -m http.server 8085 --directory MTP_Experiments/phase2_nonspherical/ellipsoid
```
Open `http://localhost:8085/ellipsoid_simulation_viewer.html` for prolate ellipsoid multi-resolution dynamics ($N=12, 42, 162, 642$).

