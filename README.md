# Binary Phase Diagram Computation via Common Tangent Construction

## Theoretical Background (Gaskell Section 10)

This code computes a temperature-composition binary phase diagram by implementing the **common tangent construction** method described in **Chapter 10 of Gaskell's *Introduction to the Thermodynamics of Materials***. The system models two components (A and B) that form solid and liquid solutions with non-ideal mixing behavior.

### Governing Thermodynamic Equations

For a binary system at temperature $T$, the **molar Gibbs free energy of mixing** for a phase $\phi$ is expressed as:

$$\Delta G^{\phi}_M = \Delta G^{\phi}_{\text{ref}} + \Delta G^{\phi}_{\text{ideal}} + \Delta G^{\phi}_{\text{excess}}$$

**Reference State Contribution:**

The reference state accounts for the Gibbs energy difference between pure components in their equilibrium states:

$$\Delta G^{\phi}_{\text{ref}} = X_A \cdot (-\Delta G_{m,A}) + X_B \cdot \Delta G_{m,B}$$

where the **Gibbs energy of melting** for each pure component is approximated by:

$$\Delta G_{m,i}(T) = \Delta H_{m,i} \cdot \frac{T_{m,i} - T}{T_{m,i}}, \quad i = A, B$$

**Ideal Mixing Contribution:**

$$\Delta G^{\phi}_{\text{ideal}} = RT \left[ X_A \ln X_A + X_B \ln X_B \right]$$

**Excess Mixing Contribution (Regular Solution Model – Gaskell §10.5):**

$$\Delta G^{\phi}_{\text{excess}} = \Omega_{\phi} \cdot X_A X_B$$

where $\Omega_{\phi}$ is the regular solution parameter (interaction parameter) for phase $\phi$.

### Complete Gibbs Energy Expressions

**Solid Phase ($S$):**

$$G^S_M(X_B, T) = -X_A \cdot \Delta G_{m,A}(T) + RT \left[ X_A \ln X_A + X_B \ln X_B \right] + \Omega_S \cdot X_A X_B$$

**Liquid Phase ($L$):**

$$G^L_M(X_B, T) = X_B \cdot \Delta G_{m,B}(T) + RT \left[ X_A \ln X_A + X_B \ln X_B \right] + \Omega_L \cdot X_A X_B$$

with $X_A = 1 - X_B$.

### Common Tangent Construction (Gaskell §10.4)

Phase coexistence at equilibrium requires:

1. **Equal chemical potentials** (equal slopes of the tangent line):

$$\left. \frac{\partial G^{S}_M}{\partial X_B} \right|_{X_B = X_B^S} = \left. \frac{\partial G^{L}_M}{\partial X_B} \right|_{X_B = X_B^L}$$

2. **Common tangent line** (equal intercepts):

$$\frac{G^L_M(X_B^L) - G^S_M(X_B^S)}{X_B^L - X_B^S} = \left. \frac{\partial G^{S}_M}{\partial X_B} \right|_{X_B = X_B^S}$$

These two conditions form the **system of nonlinear equations** solved at each temperature:

$$
\begin{cases}
\begin{aligned}
\begin{gather}
&\dfrac{\partial G^S_M}{\partial X_B}(X_B^S) - \dfrac{\partial G^L_M}{\partial X_B}(X_B^L) = 0 \tag{1}
&\dfrac{G^L_M(X_B^L) - G^S_M(X_B^S)}{X_B^L - X_B^S} - \dfrac{\partial G^S_M}{\partial X_B}(X_B^S) = 0
\end{gather}
\end{aligned}
\end{cases}
$$


### Miscibility Gaps (Gaskell §10.7)

When $\Omega > 2RT$, the Gibbs energy curve develops **two inflection points**, leading to phase separation within the same phase (solid or liquid miscibility gap). The code detects this condition and computes **solid-solid** or **liquid-liquid** coexistence.

---

## Physical Parameters

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| $T_{m,A}$ | 1135 | K | Melting temperature of pure A |
| $T_{m,B}$ | 1685 | K | Melting temperature of pure B |
| $\Delta H_{m,A}$ | 19.2 × 10³ | J/mol | Enthalpy of fusion of A |
| $\Delta H_{m,B}$ | 22.8 × 10³ | J/mol | Enthalpy of fusion of B |
| $\Omega_S$ | 13.9 × 10³ | J/mol | Solid-phase regular solution parameter |
| $\Omega_L$ | 9.1 × 10³ | J/mol | Liquid-phase regular solution parameter |
| $R$ | 8.3144 | J/(mol·K) | Universal gas constant |

---

## Numerical Implementation


### Key Features

- **Temperature Grid Resolution:** ΔT = 1 K allows precise phase boundary detection
- **Composition Grid:** $X_B \in [10^{-6}, 1 - 10^{-6}]$ avoids log(0) singularities
- **Initial Guess Strategy:** Uses intersections of Gibbs curves and local minima positions
- **Solver Configuration:** Tight tolerances (`FunctionTolerance = 10⁻¹²`, `StepTolerance = 10⁻¹²`) for accurate common tangent positions
- **Phase Detection Logic:** Minimum composition separation of `X_tol = 10⁻⁶` distinguishes true two-phase coexistence from trivial solutions

---

## Output Description

### Phase Diagram Plot
- **Red curves (Solidus):** $X_B^S(T)$ boundaries for solid phases
- **Blue curves (Liquidus):** $X_B^L(T)$ boundaries for liquid phases
- **Phase regions labeled:**
  - $\alpha_1 + \alpha_2$: Solid miscibility gap
  - $S + L$: Two-phase solid-liquid region
  - $S$: Single-phase solid solution
  - $L$: Single-phase liquid solution

### Common Tangent Verification Plot
Six representative temperatures showcase the Gibbs energy curves with superimposed common tangents, visually confirming the numerically computed equilibrium compositions.

---

## Dependencies

- **MATLAB R2020b+** (tested)
- Optimization Toolbox (`fsolve`)
- Signal Processing Toolbox (`findpeaks`)


## Contact
- Name: Mohammad Reza Hadi
- E-mail: mohammadreza.hadi.2002@gmail.com

##  
