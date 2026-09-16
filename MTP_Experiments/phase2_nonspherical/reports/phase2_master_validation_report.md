# MTP Phase 2 Master Validation Report: Non-Spherical Rigid Bodies

**Project**: Effect of Particle Shape on Low-Reynolds-Number Sedimentation and Sedimentation Stability  
**Reference Framework**: `RigidMultiblobsWall`  
**Phase**: Phase 2 (Non-Spherical Particles: Ellipsoid & Robotic Arm)  
**Date**: September 2026  
**Status**: Complete & Validated  

---

## 1. Executive Overview

Following the Phase 1 spherical baseline validation, Phase 2 establishes the theoretical and numerical qualification of non-spherical bodies within the `RigidMultiblobsWall` framework:

1. **Ellipsoid (Prolate Spheroid)**:
   - Validated against **Perrin (1934)** exact analytical solutions across aspect ratios $\lambda = a/b \in [1.0, 5.0]$.
   - Proven orthotropic decoupling at centroid ($M_{tr} = 0, M_{rt} = 0$ to $10^{-17}$).
   - Proven oblique lateral sedimentation drift $U_x(\theta) = -\frac{1}{2} F_z (\mu_\parallel - \mu_\perp) \sin(2\theta)$ with peak drift at exactly $\theta = 45^\circ$.
   - Verified zero tumbling/rotation under gravity ($\|\mathbf{\Omega}\| < 10^{-18}$).
   - Verified force linearity ($R^2 = 1.00000000$).

2. **Robotic Arm / Complex Rigid Body**:
   - Assembled from repository structures (`robot_arm_N_7`, `shell_N_12`) into a locked rigid body with full column rank geometric matrix $\text{rank}(K) = 6$.
   - Proven Onsager reciprocal symmetry ($\|M_{tr} - M_{rt}^T\|_\infty < 10^{-16}$).
   - Proven strict energy positive-definiteness ($\lambda_i > 0$).
   - Fully mapped translation-rotation coupling ($M_{tr} \ne 0$) and induced angular velocity ($\mathbf{\Omega} = M_{rt} \mathbf{F}$) when loaded away from the Center of Mobility.
   - Demonstrated fast, convergent multiblob resolution across 1-blob, 12-blob, and 42-blob segment models.

---

## 2. Summary Comparison: Sphere vs Ellipsoid vs Robotic Arm

| Hydrodynamic Property | Phase 1: Sphere ($Rh=1$) | Phase 2A: Ellipsoid ($\lambda=2$) | Phase 2B: Robotic Arm (7 Links) |
| :--- | :---: | :---: | :---: |
| **Geometry Symmetry** | Complete Isotropy ($O(3)$) | Orthotropic ($D_{\infty h}$) | Anisotropic / Broken Symmetry |
| **Mobility Tensor Form** | $\mu_0 I_{3 \times 3}$ (Isotropic) | $\text{diag}(\mu_\parallel, \mu_\perp, \mu_\perp)$ | Anisotropic $M_{tt}$ (Symmetric) |
| **Translation-Rotation Coupling** | Identically Zero ($M_{tr} \equiv 0$) | Zero at Centroid ($M_{tr} \equiv 0$) | Non-zero off CoM ($\|M_{tr}\| > 0$) |
| **Sedimentation Trajectory** | Purely Vertical ($U_x = U_y = 0$) | Oblique Drift ($U_x \ne 0$ for $\theta \ne 0$) | Oblique Drift + Reorientation / Spiraling |
| **Rotation under Pure Gravity** | Zero ($\mathbf{\Omega} \equiv \mathbf{0}$) | Zero ($\mathbf{\Omega} \equiv \mathbf{0}$) | Non-zero ($\mathbf{\Omega} = M_{rt} \mathbf{F}$) |
| **Recommended Blob Resolution** | $N = 162$ ($0.95\%$ error) | $N = 162$ ($4.8\%$ error vs Perrin) | $N = 84$ ($7 \times 12$ shells) |
| **Onsager Reciprocal Error** | $< 10^{-17}$ | $< 10^{-17}$ | $< 10^{-16}$ |

---

## 3. Key Findings & Physical Insights

### 1. The Origin of Non-Spherical Oblique Drift
In low-Reynolds-number Stokes flow, the translation of a symmetric non-spherical body is governed by anisotropic mobility:
$$\mathbf{U} = M_{tt} \mathbf{F} = \mu_\perp \mathbf{F} + (\mu_\parallel - \mu_\perp) (\mathbf{F} \cdot \hat{\mathbf{p}}) \hat{\mathbf{p}}$$
When tilted relative to gravity, the particle drifts laterally toward the streamlined direction without rotating. This lateral glide is a fundamental precursor to pair-interaction hydrodynamic drafting and collective sedimentation stability (Crowley instability).

### 2. Coupling and Induced Rotation in Complex Bodies
For an asymmetric body (or an arm tracked away from its hydrodynamic center), $M_{rt} \ne 0$. A pure settling force produces an angular velocity $\mathbf{\Omega} = M_{rt} \mathbf{F}$. This causes continuous reorientation during settling, resulting in complex helical or spiraling trajectories.

---

## 4. Repository Integrity
All experiments, data, scripts, and plots were created strictly outside the reference repository:
- `MTP_Experiments/phase2_nonspherical/`
The reference repository `RigidMultiblobsWall` remains 100% clean and unmodified.
