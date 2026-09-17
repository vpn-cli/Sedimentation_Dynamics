# Phase 2C Rigid Cylinder Audit Corrections & Caveats

This document specifies the precise corrections, terminology adjustments, and scientific caveats resulting from the Phase 2C independent audit.

---

## 1. Classification Correction: Category B (Not Category A)

* **Original Phrasing / Tendency**: Implying that the multiblob cylinder implementation was "theoretically validated against analytical theory".
* **Correction**: Phase 2C is formally classified as **Category B: Approximate theoretical benchmark + numerical convergence**.
* **Reason**: No closed-form exact analytical solution exists for the 3D Stokes resistance of a finite cylinder with flat end-caps and sharp edges. The comparison against an equivalent-surface-area sphere is an approximate geometric benchmark, not an exact theoretical solution.

---

## 2. Benchmark Clarification: Equivalent-Surface-Area Sphere

* **Original Phrasing**: Labeling $U_{\rm th} = 0.043316$ as "Analytical Theory".
* **Correction**: Replaced with **"Approximate Benchmark (Equivalent-Surface Sphere)"** across all plots, tables, and documentation.
* **Reason**: While a sphere of radius $R_{\rm eff} = \sqrt{1.5} \approx 1.22474$ shares the exact surface area $A = 6\pi$ of the cylinder, its Stokes resistance differs by small shape-dependent edge and form effects ($O(0.1\%-1\%)$).

---

## 3. Secondary Benchmark Corrections: Broersma & Tirado Formulas

* **Broersma (1960) Slender-Body Formula**:
  * **Status**: **INAPPLICABLE FOR $L/D = 1.0$**.
  * **Reason**: The slender-body asymptotic expansion requires $\sigma = \ln(L/R) > 2$ ($L/R > 7.4$). For $L/D = 1.0$ ($L/R = 2$), $\sigma = \ln(2) \approx 0.693$, causing negative resistance coefficients ($\gamma_\parallel = -9.15, \gamma_\perp = -3.75$). It must never be used for $L/D \le 2$.
* **Tirado & García de la Torre (1979) Bead-Shell Formula**:
  * **Status**: **CONTEXTUAL EXTRAPOLATION ONLY**.
  * **Reason**: Tirado's empirical polynomial fit was parameterized strictly for slender cylinders in the range $2 \le L/D \le 30$. Evaluating it at $L/D = 1.0$ yields $U = 0.041221$, which is an extrapolation outside its validated domain.

---

## 4. Force–Velocity Linearity Interpretation

* **Original Interpretation**: Presenting $R^2 = 1.000000$ as empirical confirmation that Stokes flow is linear.
* **Correction**: Explicitly state that $R^2 = 1.000000$ is a **mathematical consequence of linear algebra in the discrete RPY mobility solver** ($\mathbf{U} = \mathbf{M} \mathbf{F}$).
* **Value**: This confirms that the numerical implementation has zero nonlinear drift, perfect floating-point consistency, and zero spurious offset ($b \sim 10^{-18}$).

---

## 5. Timestep Independence Interpretation

* **Original Interpretation**: Using trajectory simulations across $dt \in \{0.20, 0.10, 0.05\}$ as proof of mobility solver accuracy.
* **Correction**: Clarify that in overdamped Stokes flow, the terminal velocity $\mathbf{U}$ is determined instantaneously by the mobility solve. Timestep variation confirms that the kinematic position integrator ($\mathbf{X}_{t+\Delta t} = \mathbf{X}_t + \mathbf{U}\Delta t$) is accurate and stable, but it does not represent an independent test of the hydrodynamic mobility matrix $\mathcal{M}$.

---

## 6. Rotational Residual Reporting

* **Original Description**: Describing the rotational velocity as "zero" or "identically zero".
* **Correction**: Quantified as a **finite numerical discretization residual**:
  $$|\mathbf{\Omega}| = 8.26 \times 10^{-6} \text{ rad/s} \implies \frac{|\mathbf{\Omega}| R}{U} \approx 1.89 \times 10^{-4} \quad (0.0189\%)$$
* **Source**: Minor discrete symmetry breaking in the polygonal spiral packing of blobs on the circular flat end-caps.

---

## 7. Plot Legend and Label Updates

All Phase 2C plots have been regenerated with updated legends:
1. `cylinder_velocity_vs_force_axial.png`: Legend updated to `Approx. Benchmark (Eq. Surface Sphere: slope=0.04332)`.
2. `cylinder_velocity_vs_force_transverse.png`: Legend updated to `Approx. Benchmark (Eq. Surface Sphere: slope=0.04332)`.
3. `cylinder_error_vs_resolution.png`: Y-axis updated to `Difference (%) vs Eq. Surface Sphere Benchmark`.
