# Engine & Injector Sizing — Theory Reference

This document explains the physics and equations behind `Engine.py`, `Chamber.py`, `IPA.py`, and `N2O.py`. It's meant to be read alongside the code: each section below corresponds to one or more functions in those files.

---

## Contents

- [Nomenclature](#nomenclature)
- [1. Engine sizing (`Engine.py`)](#1-engine-sizing-enginepy)
- [2. Chamber & nozzle geometry (`Chamber.py`)](#2-chamber--nozzle-geometry-chamberpy)
- [3. IPA injector — incompressible liquid flow (`IPA.py`)](#3-ipa-injector--incompressible-liquid-flow-ipapy)
- [4. N2O injector — choked compressible flow (`N2O.py`)](#4-n2o-injector--choked-compressible-flow-n2opy)
- [5. Assumptions & limitations](#5-assumptions--limitations)
- [6. References](#6-references)

---

## Nomenclature

| Symbol       | Meaning                                      | Units          |
|--------------|----------------------------------------------|----------------|
| $F$          | Thrust                                       | N              |
| $P_c$        | Chamber pressure                             | Pa             |
| $P_e$        | Nozzle exit pressure                         | Pa             |
| $P_a$        | Ambient pressure                             | Pa             |
| $P_0$        | Upstream (feed) stagnation pressure          | Pa             |
| $T_c$        | Chamber (combustion) temperature             | K              |
| $\gamma$     | Ratio of specific heats, $c_p/c_v$           | –              |
| $R$          | Specific gas constant, $R = \bar{R}/M$       | J/(kg·K)       |
| $\bar{R}$    | Universal gas constant (8314.46)             | J/(kmol·K)     |
| $M$          | Molecular weight                             | kg/kmol (g/mol)|
| $c^*$        | Characteristic velocity                      | m/s            |
| $C_F$        | Thrust coefficient                           | –              |
| $A_t, A_e$   | Throat / exit area                           | m²             |
| $D_t, D_e$   | Throat / exit diameter                       | m              |
| $\varepsilon$| Expansion (area) ratio, $A_e/A_t$            | –              |
| $\dot{m}$    | Mass flow rate                               | kg/s           |
| $I_{sp}$     | Specific impulse                             | s              |
| $g_0$        | Standard gravity (9.80665)                  | m/s²           |
| $L^*$        | Characteristic chamber length                | m              |
| $V_c$        | Chamber volume (injector face to throat)     | m³             |
| $C_R$        | Contraction ratio, $A_c/A_t$                | –              |
| $C_d$        | Discharge coefficient                        | –              |
| $\rho$       | Density                                      | kg/m³          |
| $A$          | Flow (orifice) area                          | m²             |

---

## 1. Engine sizing (`Engine.py`)

### 1.1 Sizing the engine based on target thrust

The engine is sized based on the **target thrust** ($F$) using the following steps:

1. **Throat Area ($A_t$)**:
   The throat area is calculated using the definition of the thrust coefficient ($C_F$):
   $$
   A_t = \frac{F}{C_F \cdot P_c}
   $$
   Here, $C_F$ depends on the nozzle geometry and pressure ratios.

2. **Mass Flow Rate ($\dot{m}$)**:
   The total mass flow rate is determined using the characteristic velocity ($c^*$):
   $$
   \dot{m} = \frac{P_c \cdot A_t}{c^*}
   $$

3. **Specific Impulse ($I_{sp}$)**:
   The specific impulse is calculated as:
   $$
   I_{sp} = \frac{F}{\dot{m} \cdot g_0}
   $$

4. **Exit Area ($A_e$)**:
   The exit area is determined using the expansion ratio ($\varepsilon$):
   $$
   A_e = \varepsilon \cdot A_t
   $$

5. **Throat and Exit Diameters**:
   The throat and exit diameters are derived from their respective areas:
   $$
   D_t = 2 \sqrt{\frac{A_t}{\pi}}, \quad D_e = 2 \sqrt{\frac{A_e}{\pi}}
   $$

These calculations are implemented in the `SizeEngine()` function.

---

### 1.2 Specific gas constant

The combustion gas's specific gas constant comes from the universal gas constant and the gas's molecular weight:

$$
R = \frac{\bar{R}}{M}
$$

Implemented in `SpecificGasConstant()`.

---

### 1.3 The Vandenkerckhove function

Several isentropic relations share a recurring combination of $\gamma$ terms, usually called the Vandenkerckhove function:

$$
\Gamma(\gamma) = \sqrt{\gamma}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{2(\gamma-1)}}
$$

Implemented in `VandenkerckhoveFunction()`.

---

### 1.4 Characteristic velocity, $c^*$

$c^*$ measures how effectively the combustion chamber converts propellant chemical energy into hot, high-pressure gas — it depends only on the propellant combination and chamber conditions, **not** on the nozzle:

$$
c^* = \frac{\sqrt{R T_c}}{\Gamma(\gamma)}
$$

Implemented in `CharacteristicVelocity()`. A real chamber never reaches the ideal value (incomplete mixing, finite reaction rate, heat loss to the walls), so the code scales it by an efficiency:

$$
c^*_{\text{real}} = \eta_{c^*} \, c^*_{\text{ideal}}, \qquad \eta_{c^*} \approx 0.90\text{–}0.98
$$

---

### 1.5 Thrust coefficient, $C_F$

$C_F$ is the "amplification" the nozzle gives to the basic $P_c A_t$ force. It has a momentum term (always positive) and a pressure term (zero only when the nozzle is *ideally expanded*, $P_e = P_a$):

$$
C_F = \underbrace{\sqrt{\frac{2\gamma^2}{\gamma-1}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}\left[1-\left(\frac{P_e}{P_c}\right)^{\frac{\gamma-1}{\gamma}}\right]}}_{\text{momentum term}} \;+\; \underbrace{\varepsilon\,\frac{P_e-P_a}{P_c}}_{\text{pressure term}}
$$

Implemented in `ThrustCoefficient()`.

---

### 1.6 Splitting the Total Mass Flow Rate

Once the total mass flow rate ($\dot{m}$) is calculated, it is split into **oxidizer** and **fuel** mass flow rates based on the oxidizer-to-fuel ratio (O/F):

$$
\dot{m}_{\text{oxidizer}} = \dot{m} \cdot \frac{\text{O/F}}{1 + \text{O/F}}
$$

$$
\dot{m}_{\text{fuel}} = \dot{m} - \dot{m}_{\text{oxidizer}}
$$

Where:
- $\dot{m}_{\text{oxidizer}}$: Oxidizer mass flow rate (e.g., N2O).
- $\dot{m}_{\text{fuel}}$: Fuel mass flow rate (e.g., IPA).
- O/F: Oxidizer-to-fuel ratio, a design parameter based on the propellant combination.

These calculations are implemented in the `SplitMassFlow()` function in `Engine.py`.

## 2. Chamber & nozzle geometry (`Chamber.py`)

`Chamber.py` is responsible for sizing the combustion chamber and nozzle geometry based on the throat area, characteristic length ($L^*$), and contraction ratio. The combustion chamber is where the propellant burns, generating high-pressure, high-temperature gas that expands through the nozzle to produce thrust.

### 2.1 Combustion chamber sizing

The combustion chamber volume is defined by the characteristic length:

$$
V_c = L^* \cdot A_t
$$

Where:
- $L^*$ is the characteristic chamber length, an empirical value depending on the propellant and combustion efficiency (typically 0.5–2 m for small liquid engines).
- $A_t$ is the throat area.

The chamber volume is split into two parts:
1. **Cylindrical section**: The main body of the chamber.
2. **Convergent cone**: The frustum-shaped section leading to the throat.

The convergent cone's volume is calculated as:

$$
V_{\text{cone}} = \frac{\pi}{3} \cdot L_{\text{conv}} \cdot \left(R_c^2 + R_c \cdot R_t + R_t^2\right)
$$

Where:
- $L_{\text{conv}}$ is the length of the convergent cone.
- $R_c$ is the chamber radius.
- $R_t$ is the throat radius.

The cylindrical section's volume is then:

$$
V_{\text{cyl}} = V_c - V_{\text{cone}}
$$

The cylindrical length is:

$$
L_{\text{cyl}} = \frac{V_{\text{cyl}}}{A_c}
$$

Where $A_c$ is the chamber cross-sectional area, determined by the contraction ratio:

$$
A_c = C_R \cdot A_t
$$

The contraction ratio ($C_R$) can either be user-defined or estimated empirically using:

$$
C_R = 8 \cdot D_t^{-0.6} + 1.25
$$

Where $D_t$ is the throat diameter in centimeters.

### 2.2 Nozzle sizing

The nozzle geometry is sized based on the expansion ratio ($\varepsilon = A_e / A_t$) and the half-angle of the divergent cone. The nozzle exit radius is:

$$
R_e = R_t \cdot \sqrt{\varepsilon}
$$

The length of the divergent cone is:

$$
L_{\text{cone}} = \frac{R_e - R_t}{\tan(\alpha)}
$$

Where $\alpha$ is the half-angle of the divergent cone (commonly 15°).

For reference, the 80%-bell-equivalent length is approximated as:

$$
L_{\text{bell, 80}} = 0.8 \cdot \frac{R_e - R_t}{\tan(15^\circ)}
$$

The divergence loss factor, which accounts for the radial component of the exit velocity in a conical nozzle, is:

$$
\lambda = \frac{1 + \cos(\alpha)}{2}
$$

### 2.3 Outputs

The `Chamber.py` module calculates and returns the following:
- **Chamber properties**:
  - Contraction ratio ($C_R$) and its source (user-defined or empirical).
  - Chamber diameter ($D_c$), cylindrical length ($L_{\text{cyl}}$), and total chamber volume ($V_c$).
  - Convergent cone length ($L_{\text{conv}}$) and volume ($V_{\text{cone}}$).
- **Nozzle properties**:
  - Exit diameter ($D_e$), divergent cone length ($L_{\text{cone}}$), and 80%-bell-equivalent length ($L_{\text{bell, 80}}$).
  - Divergence loss factor ($\lambda$).

---

## 3. IPA injector — incompressible liquid flow (`IPA.py`)

The IPA injector uses the **single-phase incompressible (SPI)** orifice equation:

$$
\dot{m} = C_d\, A\, \sqrt{2\rho\,\Delta P}, \qquad \Delta P = P_{feed} - P_c
$$

Rearranged for the required flow area:

$$
A = \frac{\dot{m}}{C_d \sqrt{2\rho\,\Delta P}}
$$

---

## 4. N2O injector — choked compressible flow (`N2O.py`)

The N2O injector uses the **choked (sonic) compressible flow** relation:

$$
\dot{m} = C_d\, A\, P_0 \sqrt{\frac{\gamma}{R\,T_0}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}}
$$

---

## 5. Assumptions & limitations

- **1-D, isentropic, ideal gas.**
- **Combustion gas properties ($T_c$, $\gamma$, $M$) must come from a thermochemistry source.**
- **N2O must actually be gaseous at the injector.**
- **Discharge coefficients ($C_d$) are inputs, subject to user decision.**

---

## 6. References

- Sutton, G. P., and Biblarz, O., *Rocket Propulsion Elements*, Wiley.
- Huzel, D. K., and Huang, D. H., *Modern Engineering for Design of Liquid-Propellant Rocket Engines*, NASA SP-125.
- Humble, R. W., Henry, G. N., and Larson, W. J., *Space Propulsion Analysis and Design*, McGraw-Hill.
- NASA CEA (Chemical Equilibrium with Applications).