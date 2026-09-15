# Phase 1 Final Validation Report: Directional Isotropy & Horizontal Force

**Study:** Single Sphere Unbounded Sedimentation under Horizontal vs Vertical Forces  
**Date:** September 15, 2026  
**Primary Resolution:** $N = 162$ blobs  
**Robustness Resolutions:** $N = 42, 162, 642$ blobs  
**Fluid Domain:** Unbounded Stokes fluid, $\eta = 1.0$ (no wall)  
**Particle Model:** Calibrated spherical shell ($R_h = 1.0$)  
**Status:** **PASSED — COMPLETE ISOTROPY VERIFIED**

---

## 1. Objective & Physical Setup

This experiment provides the final directional validation of the spherical rigid-multiblob particle before moving to nonspherical bodies:
* An ideal sphere immersed in an unbounded, quiescent Stokes fluid is **rotationally and directionally isotropic**: its hydrodynamic resistance and mobility tensors are scalar multiples of the identity matrix:
  $$ \mathbf{M}_{tt} = \frac{1}{6\pi\eta R_h} \mathbf{I} $$
* Therefore, applying an external force horizontally ($\mathbf{F}_x = [1, 0, 0]$), laterally ($\mathbf{F}_y = [0, 1, 0]$), or vertically ($\mathbf{F}_z = [0, 0, 1]$) must yield:
  1. Identical settling speeds: $|U_x| = |U_y| = |U_z| = \frac{1}{6\pi} \approx 0.0530516477$.
  2. Identical directional mobilities: $M_x = M_y = M_z$.
  3. Strictly zero transverse velocity: $U_\perp = 0$.
  4. Strictly zero angular velocity: $|\mathbf{\Omega}| = 0$.

---

## 2. Primary Case Results: $N = 162$ Blobs ($F = 1.0$)

### Velocity Components & Kinematics Table
| Applied Force Direction | $U_x$ | $U_y$ | $U_z$ | Speed $|U|$ | Theoretical $U_{\rm theory}$ | Relative Error | Alignment Angle | Max $|\mathbf{\Omega}|$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Horizontal: $\mathbf{F}_x = [1, 0, 0]$** | **0.0530534070** | $-4.05 \times 10^{-19}$ | $-1.17 \times 10^{-19}$ | 0.0530534070 | 0.0530516477 | **0.00332%** | $0.0000000000^\circ$ | $1.41 \times 10^{-18}$ |
| **Lateral: $\mathbf{F}_y = [0, 1, 0]$** | $-4.05 \times 10^{-19}$ | **0.0530534070** | $-1.76 \times 10^{-20}$ | 0.0530534070 | 0.0530516477 | **0.00332%** | $0.0000000000^\circ$ | $4.86 \times 10^{-18}$ |
| **Vertical: $\mathbf{F}_z = [0, 0, 1]$** | $1.17 \times 10^{-19}$ | $-1.76 \times 10^{-20}$ | **0.0530534070** | 0.0530534070 | 0.0530516477 | **0.00332%** | $0.0000000000^\circ$ | $8.14 \times 10^{-18}$ |

### Quantitative Comparison: Horizontal vs Vertical
* **Horizontal Speed:** $|U_x| = 0.053053407049681839$
* **Vertical Speed:** $|U_z| = 0.053053407049681853$
* **Absolute Difference:**
  $$ |U_x - U_z| = 1.3878 \times 10^{-17} $$
* **Relative Isotropy Discrepancy:**
  $$ \frac{|U_x - U_z|}{U_z} = 2.6158 \times 10^{-16} $$
  This relative difference is exactly on the order of IEEE 754 double precision machine epsilon ($\epsilon_{\rm mach} \approx 2.22 \times 10^{-16}$).
* **Transverse Velocity:**
  * Transverse magnitude under horizontal force: $U_\perp = \sqrt{U_y^2 + U_z^2} = 4.22 \times 10^{-19} \ll 10^{-12}$.
* **Angular Velocity:**
  * Magnitude under horizontal force: $|\mathbf{\Omega}| = 1.41 \times 10^{-18} \ll 10^{-12}$ rad/s.
* **Alignment:**
  * Colinearity angle between $\mathbf{U}$ and $\mathbf{F}$ is strictly $0.00^\circ$.

---

## 3. Robustness Check: Resolutions $N = 42, 162, 642$

To ensure that spatial isotropy is not an artifact of a particular mesh resolution, the exact directional test was conducted across $N = 42, 162, 642$:

| Resolution $N$ | Horizontal $M_x$ | Lateral $M_y$ | Vertical $M_z$ | Absolute Difference $|M_x - M_z|$ | Relative Discrepancy $\frac{|M_x - M_z|}{M_z}$ | Max $U_\perp$ | Max $|\mathbf{\Omega}|$ | Result |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$N = 42$** | 0.0530506033 | 0.0530506033 | 0.0530506033 | $0.00 \times 10^{0}$ | **$0.00 \times 10^{0}$** | $1.04 \times 10^{-17}$ | $1.15 \times 10^{-17}$ | **PASS** |
| **$N = 162$** | 0.0530534070 | 0.0530534070 | 0.0530534070 | $1.39 \times 10^{-17}$ | **$2.62 \times 10^{-16}$** | $4.22 \times 10^{-19}$ | $8.14 \times 10^{-18}$ | **PASS** |
| **$N = 642$** | 0.0530495835 | 0.0530495835 | 0.0530495835 | $0.00 \times 10^{0}$ | **$0.00 \times 10^{0}$** | $9.82 \times 10^{-18}$ | $1.43 \times 10^{-17}$ | **PASS** |

### Robustness Findings
1. Across all three tested resolutions, the directional mobilities $M_x$, $M_y$, and $M_z$ agree to within **$2.7 \times 10^{-16}$** relative error.
2. Spurious transverse velocities remain bounded below $1.1 \times 10^{-17}$ (well below the $10^{-12}$ threshold).
3. Spurious angular velocities remain bounded below $1.5 \times 10^{-17}$ rad/s (well below the $10^{-12}$ threshold).

---

## 4. Graphical Figures

1. **Horizontal vs Vertical Settling Speed:**  
   `plots/horizontal_vs_vertical_velocity.png`  
   Direct bar comparison illustrating identical velocities along $x, y, z$ matching theoretical Stokes settling across resolutions.
2. **Directional Mobility Isotropy & Error:**  
   `plots/directional_mobility_isotropy.png`  
   Top panel displays $M_x, M_y, M_z$ collapsing onto each other; bottom panel confirms relative discrepancy $|M_x - M_z| / M_z \le 2.62 \times 10^{-16}$.
3. **Spurious Residuals (Machine Zero Check):**  
   `plots/spurious_residuals_check.png`  
   Confirms transverse velocity and rotation magnitudes are strictly at machine noise ($10^{-19} - 10^{-17}$).

---

## 5. Final Decision: Can Sphere Validation Be Considered Complete?

### Acceptance Criteria Evaluation
* [x] **Criterion 1: Horizontal velocity aligned with horizontal force?**  
  **PASSED.** Alignment angle is $0.0000000000^\circ$ (deviation $< 10^{-14}$ deg).
* [x] **Criterion 2: Transverse velocities negligible?**  
  **PASSED.** $|U_\perp| \le 4.22 \times 10^{-19}$ (residual threshold is $10^{-12}$).
* [x] **Criterion 3: Angular velocity negligible?**  
  **PASSED.** $|\mathbf{\Omega}| \le 8.14 \times 10^{-18}$ rad/s (residual threshold is $10^{-12}$).
* [x] **Criterion 4: Horizontal speed agrees with analytical Stokes law?**  
  **PASSED.** Discrepancy is $0.00332\%$ (matches Phase 1 baseline and calibration precision).
* [x] **Criterion 5: Horizontal and vertical mobilities agree within numerical tolerance?**  
  **PASSED.** Relative discrepancy between $M_x$ and $M_z$ is $2.62 \times 10^{-16}$ (machine precision).

### Explicit Verdict
$$ \mathbf{FINAL\ VERDICT:\ PASS} $$

**The single-sphere hydrodynamic validation phase (Phase 1) is formally and completely validated.**  
The rigid-multiblob implementation in `RigidMultiblobsWall` reproduces:
1. Analytical Stokes drag in unbounded flow ($< 0.004\%$ error for calibrated shells).
2. Strict force linearity across two orders of magnitude ($R^2 = 1.0000000000$).
3. Boundary discretization convergence conforming to theoretical $O(N^{-1/2})$ scaling ($p = 0.5493$).
4. Complete directional isotropy ($M_x = M_y = M_z$ to $2.6 \times 10^{-16}$).

**The sphere phase is ready to be closed. The subsequent research phase (Phase 2: Nonspherical Particle Hydrodynamics) can proceed upon approval.**
