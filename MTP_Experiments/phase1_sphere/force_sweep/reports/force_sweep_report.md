# Phase 1B Validation Report: Force Linearity

**Configuration:** $N = 162$, $R_h = 1.0000$, $R_g = 0.9497$, $\eta = 1.0$, Unbounded Stokes Flow  
**Date:** September 15, 2026  
**Status:** **PASSED**  

## 1. Summary of Results
* **Fitted Slope ($m$):** `0.0530534070`
* **Stokes Law Slope ($m_{theory}$):** `0.0530516477`
* **Relative Slope Error:** **`0.003316%`**
* **Fitted Intercept ($b$):** `2.7384e-17`
* **$R^2$ Metric:** **`1.000000000000`**
* **Max Relative Error across sweep:** `0.003316%`
* **Transverse drift:** $< 10^{-19}$
* **Spurious rotation:** $< 10^{-17}$

## 2. Table of Numerical Values
| $F$ | $U_z$ | $|U|$ | $U_{theory}$ | Relative Error | $|\Omega|$ | $\angle(\mathbf{U}, \mathbf{F})$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.1 | 0.00530534 | 0.00530534 | 0.00530516 | 0.00332% | 8.14e-19 | 0.00e+00° |
| 0.5 | 0.02652670 | 0.02652670 | 0.02652582 | 0.00332% | 4.07e-18 | 0.00e+00° |
| 1.0 | 0.05305341 | 0.05305341 | 0.05305165 | 0.00332% | 8.14e-18 | 0.00e+00° |
| 2.0 | 0.10610681 | 0.10610681 | 0.10610330 | 0.00332% | 1.63e-17 | 0.00e+00° |
| 5.0 | 0.26526704 | 0.26526704 | 0.26525824 | 0.00332% | 4.07e-17 | 0.00e+00° |
| 10.0 | 0.53053407 | 0.53053407 | 0.53051648 | 0.00332% | 8.14e-17 | 0.00e+00° |

## 3. Plots
![Force Velocity Linear](plots/force_velocity_linear.png)
![Force Residuals](plots/force_residuals.png)
