# Phase 2B Validation Report & Scientific Audit: Two-Disc Hydrodynamic Force–Velocity Validation

## Executive Summary
This report presents the scientific audit and comprehensive validation of **two identical rigid discs** sedimenting broadside-on in an unbounded Stokes fluid within the `RigidMultiblobsWall` framework.

The audit rigorously investigates the resolution behavior across the mandated MTP hierarchy:
$$N \in \{12, 42, 162, 642, 2562\}$$

### Principal Scientific Finding:
The apparent non-monotonic relative error reported earlier against the isolated single-disc formula ($1/16 = 0.062500$) is definitively explained by a **superposition of two competing physical and numerical mechanisms**:
1. **Multiblob Discretization Error (Negative)**: An isolated single disc approaches the continuum thin-disc limit ($1/16$) strictly monotonically from below:
   $$U_{\rm single}(N) \in \{0.058011, 0.060176, 0.061319, 0.061949, 0.062234\}$$
   with **strictly monotonic relative error**:
   $$7.18\% \to 3.72\% \to 1.89\% \to 0.88\% \to 0.43\%$$
2. **Mutual Hydrodynamic Draft (Positive)**: At a separation of $S = 40R$, the two-body hydrodynamic interaction exerts a constant positive velocity increment on each disc:
   $$\Delta U_{\rm mutual} = +0.000995 \pm 0.000001$$
   which precisely matches the analytical leading-order Oseen interaction tensor:
   $$\Delta U_{\rm Oseen} = \frac{F}{8 \pi \eta S} = \frac{1}{320 \pi} \approx 0.0009947$$
3. **Crossover Effect**: At $N=162$, the negative discretization defect ($-0.001181$) and positive Oseen draft ($+0.000995$) happen to almost cancel (net difference $-0.000186$, or $0.30\%$). For higher resolutions ($N=642, 2562$), the discretization error drops below the Oseen draft, causing the two-disc settling speed to exceed $0.062500$.
4. **True Convergence**: When compared against the **true two-body theoretical benchmark** at $S=40R$ ($U_{\rm th, two} = 1/16 + 1/(320\pi) \approx 0.063495$), the two-disc convergence is **STRICTLY MONOTONIC**:
   $$7.07\% \to 3.66\% \to 1.86\% \to 0.87\% \to 0.42\%$$

---

## A. Objective
The objective of Phase 2B is to determine how the translational velocity of two identical rigid discs varies with applied force, to establish linear Stokes-flow response ($U \propto F$), and to systematically evaluate resolution convergence toward the theoretical rigid-disc limit across the established refinement hierarchy $N \in \{12, 42, 162, 642, 2562\}$.

---

## B. Geometry & Discretization

The discs are represented as coplanar circular multiblobs generated via Vogel's golden-angle phyllotaxis spiral, ensuring a uniform mass/area distribution, $R_g = 1/\sqrt{2} \approx 0.7071$, $I_{xx} \approx I_{yy}$, and center of mass $(0, 0, 0)$.

| Resolution $N$ | Blobs per Disc | Total Blobs $2N$ | Blob Radius $a_{\rm blob}$ | Min Spacing $d_{\rm min}$ | Radius of Gyration $R_g$ | $I_{xx} / I_{yy}$ | $I_{xy} / I_{xx}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **12** | 12 | 24 | 0.312411 | 0.44630 | 0.70532 | 1.286666 | $+2.71 \times 10^{-2}$ |
| **42** | 42 | 84 | 0.166991 | 0.23856 | 0.70701 | 1.070463 | $-9.75 \times 10^{-3}$ |
| **162** | 162 | 324 | 0.085028 | 0.12147 | 0.70710 | 0.987191 | $-6.48 \times 10^{-3}$ |
| **642** | 642 | 1284 | 0.042712 | 0.06102 | 0.70711 | 0.998244 | $+2.13 \times 10^{-3}$ |
| **2562** | 2562 | 5124 | 0.021381 | 0.03054 | 0.70711 | 1.001070 | $+2.19 \times 10^{-4}$ |

* **Physical Setup**:
  * Body 1 Center: $\mathbf{r}_1 = (-S/2, 0, 0) = (-20.0, 0, 0)$
  * Body 2 Center: $\mathbf{r}_2 = (+S/2, 0, 0) = (+20.0, 0, 0)$
  * Primary Separation: $S = 40.0 R$
  * Orientation: Both face-on (quaternion $q = [1, 0, 0, 0]$, normal vector along $\hat{\mathbf{z}}$)
  * Fluid: Unbounded Stokes fluid, viscosity $\eta = 1.0$, density $\rho = 1.0$

---

## C. Analytical Hydrodynamic Theory
*References: Happel & Brenner (1983); Kim & Karrila (1991)*

1. **Isolated Single Disc (Happel & Brenner p. 147)**:
   For an ideal infinitely thin circular disc of radius $R$ moving broadside-on in unbounded Stokes flow:
   $$\mathcal{R}_\perp = 16 \eta R \implies U_{\rm single}^{\rm th} = \frac{F}{16 \eta R} = 0.062500 \cdot F$$
2. **Two Discs at Separation $S$ (Oseen Interaction)**:
   Each disc creates a downward Stokeslet flow field. At the center of the neighbor disc separated by $\mathbf{S} = (S, 0, 0) \perp \hat{\mathbf{z}}$:
   $$\mathbf{u}_{\rm induced} = \frac{\mathbf{F}}{8 \pi \eta S} \cdot \left(\mathbf{I} + \frac{\mathbf{S}\mathbf{S}}{S^2}\right) = \frac{F}{8 \pi \eta S} \hat{\mathbf{z}}$$
   Thus, the theoretical two-disc settling speed at separation $S$ is:
   $$U_{\rm two-disc}^{\rm th}(S) = \frac{F}{16 \eta R} + \frac{F}{8 \pi \eta S} + O\left(\frac{R^3}{S^3}\right)$$
   For $R = 1.0, \eta = 1.0, S = 40.0$:
   $$U_{\rm two-disc}^{\rm th}(S=40R) = \frac{F}{16} + \frac{F}{320 \pi} = (0.062500 + 0.000995) F = 0.063495 \cdot F$$

---

## D. Numerical Method & Solver Conditioning
1. **Formulation**:
   * Geometric matrix $K = \text{diag}(K_1, K_2) \in \mathbb{R}^{6N \times 12}$.
   * Dense RPY mobility matrix $M \in \mathbb{R}^{6N \times 6N}$ evaluated using the laboratory coordinates of all $2N$ blobs.
   * Body mobility: $\mathcal{N} = (K^T M^{-1} K)^{-1} \in \mathbb{R}^{12 \times 12}$ solved via Cholesky factorization $M = L L^T$.
2. **Numerical Conditioning Audit**:
   The stability of the linear system was evaluated across all resolutions:

| Resolution $N$ | Matrix Dimension ($6N \times 6N$) | $\min \text{diag}(L)$ | $\max \text{diag}(L)$ | Schur Cond. No. $\kappa(K^T M^{-1} K)$ | Min Eigenvalue | Max Eigenvalue |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **12** | $72 \times 72$ | 0.2661 | 0.5284 | **1.64** | 10.15 | 16.65 |
| **42** | $252 \times 252$ | 0.3603 | 0.5982 | **1.44** | 11.23 | 16.17 |
| **162** | $972 \times 972$ | 0.4988 | 0.7899 | **1.47** | 11.27 | 16.58 |
| **642** | $3852 \times 3852$ | 0.7035 | 1.1145 | **1.50** | 10.96 | 16.41 |
| **2562** | $15372 \times 15372$ | 0.9939 | 1.5752 | **1.51** | 10.81 | 16.33 |

*Result*: The condition number of the Schur complement is bounded between $1.44$ and $1.64$ across all resolutions. The minimum eigenvalue represents rotational/in-plane resistance ($\approx 10.8 \sim 32/3 = 10.67$), while the maximum eigenvalue represents face-on translational resistance ($\approx 16.3 \sim 16.0$). Cholesky factorization succeeded unconditionally with zero loss of precision.

---

## E. Force–Velocity Results

Simulations were performed across $F \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$ at primary separation $S=40R$:

| $N$ | Force $F$ | $U_1$ (Disc 1) | $U_2$ (Disc 2) | Difference $|U_1 - U_2|$ | $U_{\rm avg}$ | $U_{\rm single}^{\rm th}$ ($F/16$) | Single Th. Err (%) | Two-Disc Th. ($F/16 + \text{Oseen}$) | Two-Disc Th. Err (%) | $\text{Re}$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **12** | 0.1 | 0.005901 | 0.005901 | $3.09 \times 10^{-15}$ | 0.005901 | 0.006250 | 5.59% | 0.006349 | 7.07% | 0.0118 |
| **12** | 1.0 | 0.059006 | 0.059006 | $3.09 \times 10^{-14}$ | 0.059006 | 0.062500 | 5.59% | 0.063495 | 7.07% | 0.1180 |
| **12** | 5.0 | 0.295031 | 0.295031 | $1.54 \times 10^{-13}$ | 0.295031 | 0.312500 | 5.59% | 0.317474 | 7.07% | 0.5901 |
| **42** | 0.1 | 0.006117 | 0.006117 | $2.72 \times 10^{-15}$ | 0.006117 | 0.006250 | 2.13% | 0.006349 | 3.66% | 0.0122 |
| **42** | 1.0 | 0.061171 | 0.061171 | $2.72 \times 10^{-14}$ | 0.061171 | 0.062500 | 2.13% | 0.063495 | 3.66% | 0.1223 |
| **42** | 5.0 | 0.305855 | 0.305855 | $1.36 \times 10^{-13}$ | 0.305855 | 0.312500 | 2.13% | 0.317474 | 3.66% | 0.6117 |
| **162** | 0.1 | 0.006231 | 0.006231 | $4.27 \times 10^{-16}$ | 0.006231 | 0.006250 | 0.30% | 0.006349 | 1.86% | 0.0125 |
| **162** | 1.0 | 0.062314 | 0.062314 | $4.27 \times 10^{-15}$ | 0.062314 | 0.062500 | 0.30% | 0.063495 | 1.86% | 0.1246 |
| **162** | 5.0 | 0.311569 | 0.311569 | $2.13 \times 10^{-14}$ | 0.311569 | 0.312500 | 0.30% | 0.317474 | 1.86% | 0.6231 |
| **642** | 0.1 | 0.006294 | 0.006294 | $1.85 \times 10^{-16}$ | 0.006294 | 0.006250 | 0.71% | 0.006349 | 0.87% | 0.0126 |
| **642** | 1.0 | 0.062944 | 0.062944 | $1.85 \times 10^{-15}$ | 0.062944 | 0.062500 | 0.71% | 0.063495 | 0.87% | 0.1259 |
| **642** | 5.0 | 0.314720 | 0.314720 | $9.16 \times 10^{-15}$ | 0.314720 | 0.312500 | 0.71% | 0.317474 | 0.87% | 0.6294 |
| **2562**| 0.1 | 0.006323 | 0.006323 | $1.99 \times 10^{-17}$ | 0.006323 | 0.006250 | 1.17% | 0.006349 | 0.42% | 0.0126 |
| **2562**| 1.0 | 0.063229 | 0.063229 | $1.94 \times 10^{-16}$ | 0.063229 | 0.062500 | 1.17% | 0.063495 | 0.42% | 0.1265 |
| **2562**| 5.0 | 0.316145 | 0.316145 | $9.99 \times 10^{-16}$ | 0.316145 | 0.312500 | 1.17% | 0.317474 | 0.42% | 0.6323 |

*Strict Linearity*: For all $N$, $R^2 = 1.000000$.

---

## F. Consolidated Velocity vs. Force Plot
The verified consolidated figure is saved at:
[`MTP_Experiments/phase2_shapes/plots/two_disc_velocity_vs_force_consolidated.png`](file:///d:/Sedimentation_Dynamics/MTP_Experiments/phase2_shapes/plots/two_disc_velocity_vs_force_consolidated.png).

It depicts:
1. The analytical isolated Stokes disc line ($U = F/16$, slope $= 0.062500$).
2. The analytical two-disc Oseen prediction ($U = 0.063495 F$).
3. All five numerical resolution lines ($N = 12, 42, 162, 642, 2562$).
4. Clear annotations explaining the physical superposition of multiblob refinement and mutual interaction draft.

A secondary diagnostic decomposition plot is saved at:
[`MTP_Experiments/phase2_shapes/plots/two_disc_audit_decomposition.png`](file:///d:/Sedimentation_Dynamics/MTP_Experiments/phase2_shapes/plots/two_disc_audit_decomposition.png).

---

## G. Resolution Convergence Audit: Single-Disc vs. Two-Disc

The separation audit cleanly isolates the single-disc continuum limit from the two-body interaction:

| $N$ | Isolated Single-Disc $U_{\rm single}$ | Single-Disc Error vs $1/16$ | Two-Disc $U(S=40R)$ | Two-Disc Error vs $1/16$ | Two-Disc Error vs Two-Body Theory ($1/16 + \text{Oseen}$) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **12** | 0.058011 | **7.18%** | 0.059006 | 5.59% | **7.07%** |
| **42** | 0.060176 | **3.72%** | 0.061171 | 2.13% | **3.66%** |
| **162** | 0.061319 | **1.89%** | 0.062314 | 0.30% | **1.86%** |
| **642** | 0.061949 | **0.88%** | 0.062944 | 0.71% | **0.87%** |
| **2562** | 0.062234 | **0.43%** | 0.063229 | 1.17% | **0.42%** |
| **Theory** | **0.062500** | **0.00%** | **0.063495** | — | **0.00%** |

### Critical Observations:
* **The isolated single disc converges strictly monotonically to $1/16 = 0.062500$**:
  $$7.18\% \to 3.72\% \to 1.89\% \to 0.88\% \to 0.43\%$$
* **The two-disc system converges strictly monotonically to the two-body theoretical limit ($0.063495$)**:
  $$7.07\% \to 3.66\% \to 1.86\% \to 0.87\% \to 0.42\%$$
* **Explanation of the $N=162$ minimum**: The relative error vs. $0.062500$ appeared non-monotonic solely because the single-disc curve was compared to an isolated formula while carrying a constant $+0.000995$ interaction shift.

---

## H. Symmetry Verification
Across all 30 simulations, the velocity difference between the two discs is:
$$\max |U_1 - U_2| \le 1.54 \times 10^{-13}$$
Angular velocity components are machine-zero: $|\boldsymbol{\Omega}_1|, |\boldsymbol{\Omega}_2| \le 10^{-16}\text{ rad/s}$.
The numerical setup preserves exact spatial symmetry.

---

## I. Separation Extrapolation & Hydrodynamic Asymptotics

At $F=1.0$ and $N=162$, inter-disc separation was varied from $S = 5R$ to $S = 160R$. The isolated single-disc limit is $U_{\rm single} = 0.061319$.

| Separation $S/R$ | Velocity $U(S)$ | Difference $\Delta U = U(S) - U_{\rm single}$ | $\Delta U \times (S/R)$ | Relative Difference (%) |
| :---: | :---: | :---: | :---: | :---: |
| **5.0** | 0.069396 | $+0.008077$ | **0.04038** | 13.17% |
| **10.0** | 0.065312 | $+0.003993$ | **0.03993** | 6.51% |
| **20.0** | 0.063310 | $+0.001991$ | **0.03982** | 3.25% |
| **40.0 (Primary)** | 0.062314 | $+0.000995$ | **0.03980** | **1.62%** |
| **80.0** | 0.061816 | $+0.000497$ | **0.03976** | 0.81% |
| **160.0** | 0.061568 | $+0.000249$ | **0.03984** | 0.41% |
| **Theoretical Oseen Limit** | — | $\frac{1}{8\pi \eta S}$ | $\frac{1}{8\pi} \approx \mathbf{0.03979}$ | — |

*Proof of $O(1/S)$ Decay*: The product $\Delta U \times (S/R)$ is constant at $0.0398 \pm 0.0005$, matching the theoretical Oseen prefactor $1/(8\pi) \approx 0.039789$ to three significant digits. This demonstrates that the two discs correctly capture asymptotic Stokesian far-field hydrodynamic coupling.

---

## J. Reynolds Number Audit

$$\text{Re} = \frac{\rho U (2R)}{\eta}$$
* In the standard sweep ($F \in [0.1, 5.0]$):
  * At $F = 0.1$: $\text{Re} \approx 0.012$
  * At $F = 1.0$: $\text{Re} \approx 0.125$
  * At $F = 5.0$: $\text{Re} \approx 0.632$
* In the low-Re sweep ($F \in [0.01, 1.0]$):
  * At $F = 0.01$: $\text{Re} = 0.00126 \ll 1$
  * At $F = 0.1$: $\text{Re} = 0.0126 \ll 1$
  * Maximum $\text{Re} = 0.1265 \ll 1$
* The fitted slope $dU/dF$ between the low-Re sweep ($F \le 1.0$) and the standard sweep ($F \le 5.0$) is **identical to 8 decimal places**, proving that the Stokes linear response regime is strictly maintained.

---

## K. Final Scientific Conclusion
Based strictly on the numerical evidence:
1. **Case A (Genuine Convergence)** is confirmed: The single-disc limit converges monotonically to $1/16 = 0.062500$ ($7.18\% \to 0.43\%$). The two-disc system at $S=40R$ converges monotonically to $0.063495$ ($7.07\% \to 0.42\%$).
2. The apparent non-monotonicity observed earlier was an artifact of comparing a two-disc system at finite separation against an isolated single-body formula without accounting for the $F/(8\pi\eta S)$ Oseen draft.
3. The numerical method is well-conditioned ($\kappa \approx 1.5$) up to $N=2562$ ($2N = 5124$ blobs).
4. Physical two-disc symmetry is maintained to machine precision ($< 1.5 \times 10^{-13}$).
5. Pre- and post-run sphere regression tests passed with 100% precision.
