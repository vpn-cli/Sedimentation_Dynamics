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
│   ├── phase2_nonspherical/       # Ellipsoid & non-spherical body hydrodynamics
│   │   ├── ellipsoid/
│   │   │   ├── ellipsoid_simulation_viewer.html  # Interactive 3D WebGL Multi-Resolution Viewer
│   │   │   ├── plots/                           # 300 DPI publication-ready figures
│   │   │   ├── processed_data/                  # Discretized meshes & trajectory JSON data
│   │   │   └── reports/                         # Markdown validation reports
│   └── scripts/                   # Simulation, validation, and analysis scripts
└── RigidMultiblobsWall/           # Reference multiblob Stokesian hydrodynamics library
```

---

## 🚀 Running Simulations & Interactive 3D Viewer

### 1. Run Dynamic Sedimentation Trajectories
```bash
python MTP_Experiments/scripts/phase2_ellipsoid_trajectory_sim.py
```

### 2. Run Force Sweep & Resolution Convergence
```bash
python MTP_Experiments/scripts/phase2_ellipsoid_force_resolution.py
```

### 3. Launch Interactive 3D Multi-Resolution Viewer
Start an HTTP server in the ellipsoid experiment directory:
```bash
python -m http.server 8085 --directory MTP_Experiments/phase2_nonspherical/ellipsoid
```
Open your browser at `http://localhost:8085/ellipsoid_simulation_viewer.html` to:
- Dynamically toggle blob resolutions ($N = 12, 42, 162, 642$).
- Inspect real-time 3D settling trajectories, velocity vectors, and telemetry HUD.
- Compare descent dynamics across orientations ($\theta = 0^\circ$ broadside, $\theta = 45^\circ$ gliding, $\theta = 90^\circ$ streamlined).
