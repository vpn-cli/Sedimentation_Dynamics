# Phase 2 Validation Report: Robotic Arm / Complex Rigid Body Hydrodynamics

**Geometry**: 7-Segment Articulated / Rigid Robotic Arm  
**Solver**: RigidMultiblobsWall RPY / Cholesky  
**Date**: 2026-09-15 19:47:46  

---

## 1. Executive Summary

This report demonstrates the representation and hydrodynamic validation of an arbitrary complex non-spherical body—a multi-segment robotic arm—using the `RigidMultiblobsWall` framework.
Key findings:
1. **Arbitrary Rigid Multiblob Connectivity**: Spherical links are successfully assembled into a unified rigid particle. The kinematic matrix $K$ has full column rank ($\text{rank}(K) = 6$), enforcing exact rigid-body motion across all blobs.
2. **Onsager Reciprocal Symmetry**: The $6 \times 6$ generalized mobility tensor $\mathcal{N}_{body}$ satisfies Onsager reciprocity $\|M_{tr} - M_{rt}^T\|_\infty < 10^{-16}$ to machine precision.
3. **Energy Positive-Definiteness**: All 6 eigenvalues of $\mathcal{N}_{body}$ are strictly positive ($\lambda_i > 0$), confirming thermodynamic consistency.
4. **Center of Mobility & Induced Rotation**: When the reference point is placed away from the Center of Mobility (e.g. at the root joint or for a bent arm), settling under pure gravity induces a nonzero angular velocity $\mathbf{\Omega} = M_{rt} \mathbf{F} \ne \mathbf{0}$, causing the arm to reorient and spiral during sedimentation.
5. **Resolution Verification**: The recommended 12-blob per link discretization ($N_{total} = 84$ blobs) provides accurate multiblob hydrodynamics with fast runtime (< 0.05 s).

---

## 2. Rigid Assembly & Kinematic Properties

| Parameter | Value | Notes |
| :--- | :--- | :--- |
| Number of Links | `7` | Sequentially connected along arm |
| Blobs per Link | `12` | Discretized with `shell_N_12` |
| Total Blobs $N_{total}$ | `84` | $7 \times 12 = 84$ blobs |
| Kinematic Matrix Dimension | `252 x 6` | $(3 N_{blobs}) \times 6$ |
| Rank of $K$ Matrix | `6` | **Full column rank (6)** |
| Condition Number $\kappa(K)$ | `7.7956` | Well-conditioned |
| Total Arm Length | `16.42` | Center-to-center span |

---

## 3. Mobility Tensor Symmetries & Onsager Reciprocity

| Metric | Numerical Value | Target | Status |
| :--- | :--- | :--- | :--- |
| Onsager Reciprocal Error $\|M_{tr} - M_{rt}^T\|_\infty$ | `2.1411e-17` | `< 1e-14` | **PASS (Machine Precision)** |
| Symmetric Positive Definite | `True` | `True` | **PASS** |
| Minimum Eigenvalue $\lambda_{min}$ | `3.6069e-04` | `> 0` | **PASS** |
| Longitudinal Mobility $\mu_{xx}$ | `0.020099` | — | Along arm length |
| Transverse Mobility $\mu_{yy} = \mu_{zz}$ | `0.014922` | — | Broadside-on |

---

## 4. Center of Mobility & Translation-Rotation Coupling

| Configuration | $\|M_{{tr}}\|$ | CoM Shift $[\Delta x, \Delta y]$ | Induced Angular Velocity $\|\mathbf{\Omega}\|$ |
| :--- | :---: | :---: | :---: |
| Straight Arm (Centroid) | `2.1527e-17` | `[+0.00, +0.00]` | `4.7110e-19` |
| Straight Arm (Root) | `3.8262e-03` | `[+7.50, +0.00]` | `2.7052e-03` |
| Bent Arm 45° (Centroid) | `5.2954e-04` | `[-0.02, +0.17]` | `4.7311e-04` |
| Bent Arm 45° (Root) | `6.4805e-03` | `[+6.44, +2.70]` | `5.8332e-03` |
| Bent Arm 90° (Centroid) | `5.5132e-04` | `[-0.14, +0.33]` | `4.0750e-04` |
| Bent Arm 90° (Root) | `7.1264e-03` | `[+3.79, +3.90]` | `6.4928e-03` |

![Coupling Comparison](../plots/robotic_arm_coupling_comparison.png)

> [!IMPORTANT]
> When the reference tracking point is at the Center of Mobility of a straight symmetric arm, $\|M_{tr}\| = 0$ and gravity produces zero rotation.
> When the arm is bent or tracked at the root, non-zero $M_{tr} = M_{rt}^T$ couples sedimentation force to angular rotation $\mathbf{\Omega}$.

---

## 5. Resolution Comparison

| Model | Blobs/Link | Total Blobs | $\mu_{{xx}}$ | $\mu_{{yy}}$ | Anisotropy Ratio | Runtime (s) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Single Blob per Link | 1 | 7 | 0.020484 | 0.015083 | 1.358 | 0.000s |
| 12-Blob Shell per Link (Recommended) | 12 | 84 | 0.020099 | 0.014922 | 1.347 | 0.005s |
| 42-Blob Shell per Link (High-Res) | 42 | 294 | 0.020052 | 0.014919 | 1.344 | 0.041s |

![Resolution Comparison](../plots/robotic_arm_resolution_comparison.png)

---

## 6. Conclusions & Recommendations

1. **Complex Arbitrary Bodies Validated**: The repository's rigid multiblob formulation seamlessly generalizes from simple spheroids to complex articulated/chain structures.
2. **Physical Kinematics Guaranteed**: Exact rigid connectivity is maintained via the $K$ matrix without any spurious deformation.
3. **Coupling Fully Characterized**: Translation-rotation coupling ($M_{tr}$) is quantitatively mapped to geometric asymmetry and tracking point location.