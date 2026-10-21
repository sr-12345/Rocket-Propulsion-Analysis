# Engine & Injector Sizing: Theory Reference
 
This document explains the physics behind the code and how each section fits together.
 
The overview of the tool:
 
> You input how much thrust you want and at what chamber pressure. The code calculates how big the throat has to be, how much propellant that burns per second, how big the chamber and nozzle are, and finally how big to drill the injector holes so that exactly that much IPA and N₂O flows into the engine.

The calculations work in a chain. Each stage produces numbers the next stage needs:
 
User inputs:
 Target thrust F, chamber pressure Pc, ambient Pa
 + combustion gas properties (Tc, γ, M) from CEA/RPA

User inputs are fed into Engine.py:
Calculates c* and CF to compute the throat area At
From At and c* we calculate total mass flow ṁ
Using ṁ and the O/F ratio we can then find ṁ_ox, ṁ_fuel

Engine.py feeds At into Chamber.py:
Calculates L* and At to get chamber volume
Chamber diameter/length computed using contraction ratio

ṁ_ox, ṁ_fuel fed to N2O.py and IPA.py:
Required flow areas are calculated from each mass flow rate and the users choice of Cd
From required flow areas we calculate injection orifice sizes based off the fuel flow type

 
Two separate groups of physics are involved, and it helps to keep them apart in your head:
 
1. **The hot side** (`Engine.py`, `Chamber.py`) describes the *combustion products* (hot gas at ~2500 K) flowing through the chamber and nozzle. Its properties (`Tc`, `Gamma_gas`, `MW_gas`) come from a thermochemistry program.
2. **The cold side** (`IPA.py`, `N2O.py`) describes the *propellants before combustion* flowing through the injector orifices. Its properties (`Rho_IPA`, `Gamma_N2O`, `R_N2O`, `N2O_Temp`) describe the liquid IPA and the nitrous oxide.
The two sides only meet through the **mass flow rate** and the **chamber pressure**: the hot side decides how much propellant per second is needed, and the cold side decides what holes deliver it against that chamber pressure.
---
 
## Contents

- [Nomenclature](#nomenclature)
- [1. Engine sizing (`Engine.py`)](#1-engine-sizing-enginepy)
- [2. Chamber and nozzle geometry (`Chamber.py`)](#2-chamber-and-nozzle-geometry-chamberpy)
- [3. Injector sizing: the IPA side (`IPA.py`)](#3-injector-sizing-the-ipa-side-ipapy)
- [4. Injector sizing: the N₂O side (`N2O.py`)](#4-injector-sizing-the-n₂o-side-n2opy)
- [5. Turning flow area into drill sizes](#5-turning-flow-area-into-drill-sizes)
- [6. A worked example](#6-a-worked-example)
- [7. Assumptions and limitations](#7-assumptions-and-limitations)
- [8. References](#8-references)
 
---
 
## Nomenclature
 
| Symbol | Meaning | Units |
|---|---|---|
| $F$ | Thrust | N |
| $P_c$ | Chamber pressure | Pa |
| $P_e$ | Nozzle exit pressure | Pa |
| $P_a$ | Ambient pressure | Pa |
| $P_0$ | Upstream (feed) pressure at an injector | Pa |
| $T_c$ | Chamber (combustion) temperature | K |
| $\gamma$ | Ratio of specific heats, $c_p/c_v$ | – |
| $R$ | Specific gas constant, $\bar R/M$ | J/(kg·K) |
| $\bar R$ | Universal gas constant (8314.46) | J/(kmol·K) |
| $M$ | Molecular weight | kg/kmol (= g/mol) |
| $c^*$ | Characteristic velocity | m/s |
| $C_F$ | Thrust coefficient | – |
| $\eta_{c^*},\ \eta_{C_F}$ | Efficiency factors applied to $c^*$ and $C_F$ | – |
| $A_t,\ A_e,\ A_c$ | Throat / exit / chamber cross-section area | m² |
| $D_t,\ D_e,\ D_c$ | Throat / exit / chamber diameter | m |
| $\varepsilon$ | Expansion ratio, $A_e/A_t$ | – |
| $C_R$ | Contraction ratio, $A_c/A_t$ | – |
| $\dot m$ | Mass flow rate | kg/s |
| $I_{sp}$ | Specific impulse | s |
| $g_0$ | Standard gravity (9.80665) | m/s² |
| $L^*$ | Characteristic chamber length | m |
| $V_c$ | Chamber volume (injector face to throat) | m³ |
| $\text{O/F}$ | Oxidizer-to-fuel mass ratio | – |
| $C_d$ | Discharge coefficient | – |
| $\rho$ | Density | kg/m³ |
| $A$ | Orifice flow area | m² |
 
---
 
## 1. Engine sizing (`Engine.py`)
 
### 1.1 The core idea
 
A rocket nozzle converts the pressure and heat of the chamber gas into a fast jet of propellants. Two constants are key:
 
- **$c^*$ (characteristic velocity)** says how good the *combustion* is: how much chamber pressure you get per unit of propellant flow through a given throat area. It depends on the propellants and chamber conditions, not on the nozzle shape.
- **$C_F$ (thrust coefficient)** says how good the *nozzle* is: how much thrust you get per unit of chamber pressure acting on the throat area. It depends on the nozzle's expansion and the pressures, and only indirectly on the propellants (via $\gamma$).
Their definitions are
 
$$
c^* = \frac{P_c A_t}{\dot m}, \qquad C_F = \frac{F}{P_c A_t}
$$
 
and rearranging them gives the sizing equation directly:
 
$$
\boxed{A_t = \frac{F}{C_F\,P_c}}, \qquad \boxed{\dot m = \frac{P_c A_t}{c^*}}, \qquad \boxed{I_{sp} = \frac{F}{\dot m\, g_0}}
$$
 
So the job of `Engine.py` is to find $c^*$ and $C_F$ from the gas properties; everything else follows. This is done in `SizeEngine()`.
 
### 1.2 Gas constant and the Vandenkerckhove function
 
The combustion gas is assumed to behave as an ideal gas with
 
$$
R = \frac{\bar R}{M}
$$
 
(`SpecificGasConstant()`). A particular combination of $\gamma$ terms appears frequently in isentropic nozzle flow, so it is defined by the Vandenkerckhove function:
 
$$
\Gamma(\gamma) = \sqrt{\gamma}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{2(\gamma-1)}}
$$
 
(`VandenkerckhoveFunction()`). Note this $\Gamma$ is a *function of* $\gamma$ and is not the same thing as $\gamma$.
 
### 1.3 Characteristic velocity, $c^*$
 
For a choked throat, the mass flow through it is fixed by chamber conditions, which leads to
 
$$
c^*_{\text{ideal}} = \frac{\sqrt{R\,T_c}}{\Gamma(\gamma)}
$$
 
(`CharacteristicVelocity()`). Physically, hotter or lighter gas (high $T_c$, low $M$) gives a higher $c^*$.
 
A real chamber never achieves this: propellants don't mix perfectly, reactions aren't instantaneous, and some heat dissipates into the walls. The code applies an efficiency:
 
$$
c^*_{\text{real}} = \eta_{c^*}\,c^*_{\text{ideal}}, \qquad \eta_{c^*}\approx 0.90\text{–}0.98
$$
 
Small engines tend to sit at the low end of that range.
 
### 1.4 Thrust coefficient, $C_F$
 
The thrust coefficient has two parts:
 
$$
C_F = \underbrace{\sqrt{\frac{2\gamma^2}{\gamma-1}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}\left[1-\left(\frac{P_e}{P_c}\right)^{\frac{\gamma-1}{\gamma}}\right]}}_{\text{momentum term}} \;+\; \underbrace{\varepsilon\,\frac{P_e-P_a}{P_c}}_{\text{pressure term}}
$$
 
- The **momentum term** is the thrust from accelerating the gas. It is always positive and grows the more the gas is expanded (lower $P_e/P_c$).
- The **pressure term** is the extra push (or drag) when the exit pressure doesn't match the surroundings. It is zero when $P_e = P_a$, i.e. the nozzle is *ideally expanded*.
`ThrustCoefficient()` evaluates both. By default the code designs for ideal expansion at the ambient pressure you supply ($P_e = P_a$), so the pressure term vanishes and the nozzle is perfect at that one altitude. If you would rather fix the nozzle's area ratio, `EpsilonToPe()` finds the exit pressure that corresponds to a chosen $\varepsilon$ by bisection, and you can pass that in as `P_e`.
 
### 1.5 Finding the expansion ratio
 
$C_F$ needs the area ratio $\varepsilon$, and $\varepsilon$ follows from the pressure ratio via the exit Mach number:
 
1. From the isentropic relation, the pressure ratio fixes the exit Mach number $M_e$ (`ExitMachFromPressureRatio()`):
$$
M_e = \sqrt{\frac{2}{\gamma-1}\left[\left(\frac{P_e}{P_c}\right)^{-\frac{\gamma-1}{\gamma}}-1\right]}
$$
2. The Mach number fixes the area ratio (`AreaRatioFromMach()`):
$$
\varepsilon = \frac{1}{M_e}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M_e^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}
$$
Only the supersonic solution ($M_e>1$, $P_e/P_c<1$) is relevant for a rocket nozzle.
 
### 1.6 Applying nozzle losses and sizing the throat
 
Like $c^*$, the real $C_F$ falls short of the ideal one:
 
$$
C_{F,\text{real}} = \eta_{C_F}\,C_{F,\text{ideal}}
$$
 
`SizeEngine()` then runs the sequence of calculations from 1.1 with the *real* values, so the throat is sized to deliver the target thrust *after* losses. Throat and exit sizes come from the areas:
 
$$
A_e = \varepsilon A_t, \qquad D_t = 2\sqrt{\frac{A_t}{\pi}}, \qquad D_e = 2\sqrt{\frac{A_e}{\pi}}
$$
 
> **Note on $\eta_{C_F}$ and the divergence factor.** A conical nozzle also loses performance because some exit velocity points sideways (see §2.3). `Chamber.py` computes that divergence factor $\lambda$ for reference, but it is *not* automatically multiplied into $C_F$. If you want it accounted for, include it during your choice of `eta_cf`.
 
### 1.7 Splitting the mass flow between propellants
 
The total flow $\dot m$ is shared between oxidizer and fuel according to the O/F ratio ($\text{O/F} = \dot m_{ox}/\dot m_{fuel}$):
 
$$
\dot m_{ox} = \dot m\,\frac{\text{O/F}}{1+\text{O/F}}, \qquad \dot m_{fuel} = \frac{\dot m}{1+\text{O/F}}
$$
 
(`SplitMassFlow()`). These two numbers are fed to the injector sections below. The O/F you choose **must** be the same O/F you used when looking up $T_c$, $\gamma$ and $M$ in CEA, because those gas properties depend on it.
 
---
 
## 2. Chamber and nozzle geometry (`Chamber.py`)
 
The engine section gave a throat area. This section builds the chamber upstream of it and the nozzle downstream of it.
 
### 2.1 How big should the chamber be?
 
Propellant needs time to mix and burn before the gas reaches the throat. The longer it spends in the chamber, the more complete the combustion. Instead of specifying that time directly, we use the **characteristic length**:
 
$$
L^* = \frac{V_c}{A_t} \quad\Longleftrightarrow\quad V_c = L^*\,A_t
$$
 
$L^*$ is the length a chamber of the same volume would have if it were a straight tube with the throat's cross-section. It is empirical and propellant-dependent; 0.5–2 m is typical for small liquid engines, and higher is generally more complete combustion but more weight and more heat-transfer area.
 
As a cross-check, the gas **stay time** in the chamber follows from the volume, the gas density and the flow rate:
 
$$
\tau = \frac{\rho V_c}{\dot m} = \frac{L^*\,c^*}{R\,T_c}
$$
 
(`StayTime()`). Stay times of a few milliseconds are common for small engines.
 
### 2.2 Chamber shape
 
The chamber is a cylinder followed by a conical convergent section that funnels down to the throat. Its size is set in four steps (`SizeChamber()`):
 
1. **Chamber area from the contraction ratio.** Choose $C_R = A_c/A_t$, then
$$
A_c = C_R\,A_t, \qquad R_c = \sqrt{A_c/\pi}
$$
If you don't pick a value, an empirical curve fit for small engines is used ($D_t$ in cm). It is only a starting point; a value you can justify is better.
$$
C_R \approx 8\,D_t^{-0.6} + 1.25
$$
(`ContractionRatioEstimate()`). Be aware that for very small throats this fit can return large values (above 10 for a throat of around 8 mm), which is well outside the usual 3–6 design range.
2. **Convergent cone length** from the chosen half-angle $\beta$ (default 45°):
$$
L_{conv} = \frac{R_c - R_t}{\tan\beta}
$$
The half-angle is the angle between the axis of the chamber and the wall of the convergent section leading to the throat. It determines how quickly the chamber radius decreases to the throat diameter, with smaller half-angles allowing a smoother transition for the flow and avoiding large pressure losses. Larger half-angles can cause flow separation from the chamber walls and turbulence however reduce the engine's length and weight.
3. **Convergent cone volume.** The convergent cone is a conical frustum, so
$$
V_{conv} = \frac{\pi}{3}\,L_{conv}\left(R_c^2 + R_cR_t + R_t^2\right)
$$
4. **Cylinder length.** Whatever volume is left over after the cone goes into the cylinder:
$$
V_{cyl} = V_c - V_{conv}, \qquad L_{cyl} = \frac{V_{cyl}}{A_c}
$$
If $V_{conv}$ already exceeds $V_c$, there is no room for a cylinder and the code raises an error. The fixes are to raise $L^*$, lower $C_R$, or steepen the convergent half-angle.
 
The code reports $L_{cyl}/D_c$ as a sanity check; the GUI flags values outside roughly 0.5–4 (typical designs are around 1–3).
 
### 2.3 Nozzle (divergent) section
 
The exit radius follows directly from the expansion ratio found in §1.5:
 
$$
R_e = R_t\sqrt{\varepsilon}
$$
 
A simple conical nozzle with half-angle $\alpha$ (commonly 15°) has length
 
$$
L_{cone} = \frac{R_e - R_t}{\tan\alpha}
$$
 
A real **bell** nozzle reaches the same area ratio in a shorter length. The common "80% bell" rule of thumb takes 80% of the length of an equivalent 15° cone:
 
$$
L_{bell,80} = 0.8\,\frac{R_e - R_t}{\tan 15^\circ}
$$
 
The price of a cone is **divergence loss**: the exhaust leaves with a sideways component that doesn't contribute to thrust. The correction factor for a conical nozzle is
 
$$
\lambda = \frac{1+\cos\alpha}{2}
$$
 
(for 15° this is about 0.983, a 1.7% loss). `SizeNozzle()` returns all of these.
 
---
 
## 3. Injector sizing: the IPA side (`IPA.py`)
 
Now we have a fuel mass flow, $\dot m_{fuel}$, and a chamber pressure. The injector's job is to deliver that flow with the feed system pressure pushing against the chamber pressure.
 
IPA is a liquid, essentially incompressible, so the **single-phase incompressible (SPI)** orifice equation applies. It comes from Bernoulli's equation: the pressure drop across the orifice is converted into velocity, and the discharge coefficient accounts for real-world losses and flow contraction.
 
$$
\dot m = C_d\,A\,\sqrt{2\rho\,\Delta P}, \qquad \Delta P = P_{feed} - P_c
$$
 
Rearranged for the area we need (`SPI()`):
 
$$
A = \frac{\dot m}{C_d\sqrt{2\rho\,\Delta P}}
$$
 
Things worth noticing:
 
- Flow area scales as $1/C_d$. A sloppier orifice (low $C_d$) needs a bigger hole for the same flow.
- Area scales as $1/\sqrt{\Delta P}$: doubling the pressure drop only shrinks the hole by about 30%.
- $\Delta P$ is the pressure *across the injector*. If feed pressure doesn't comfortably exceed chamber pressure (injector stiffness of 15–30% of $P_c$ is a common target), combustion pressure fluctuations can feed back into the flow.
---
 
## 4. Injector sizing: the N₂O side (`N2O.py`)
 
The N₂O side is treated as a **compressible gas** that flows through the orifice at **sonic (choked) velocity**.
 
### 4.1 Why choking?
 
When the pressure ratio across an orifice is large enough, the gas reaches the speed of sound at the narrowest point. After that, lowering the downstream pressure further doesn't increase the flow; the mass flow depends only on the upstream conditions. This is convenient for injectors because chamber pressure fluctuations can't travel upstream through a choked orifice.
 
The condition is (`ChokedFlowCriteria()`):
 
$$
\frac{P_c}{P_{feed}} \le \left(\frac{2}{\gamma+1}\right)^{\frac{\gamma}{\gamma-1}}
$$
 
For $\gamma \approx 1.3$ the right-hand side is about 0.55, so the feed pressure has to be roughly 1.8× the chamber pressure or more. The code reports whether this is met, and the GUI shows a warning if it isn't.
 
### 4.2 The choked-flow equation
 
For an ideal gas the choked mass flow is
 
$$
\dot m = C_d\,A\,P_0\sqrt{\frac{\gamma}{R\,T_0}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}}
$$
 
where $P_0$ and $T_0$ are the **upstream** (feed-side) pressure and temperature. Solved for the area (`ChokedNitrousAreaRequired()`):
 
$$
A = \frac{\dot m}{C_d\,P_0\sqrt{\dfrac{\gamma}{R\,T_0}\left(\dfrac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}}}
$$
 
> **Which pressure is $P_0$?** It must be the N₂O **feed** pressure, not the chamber pressure. Choked flow does not depend on downstream pressure at all, which is the whole point. The GUIs pass the feed pressure (`P_feedsystem_N2O`). Be careful when calling the function from your own scripts, where `P_c` would silently give an orifice far too large.
 
### 4.3 Limitations specific to N₂O
 
Nitrous oxide is a liquid at typical tank conditions and has a vapour pressure of roughly 50 bar at 290 K. The gas formula above is only appropriate if the N₂O really is gaseous at the injector. If it's liquid or two-phase (which is common for a pressure-fed system with high feed pressures), the flow behaves quite differently and a model such as the homogeneous equilibrium model (HEM) or Dyer's model is more appropriate. Treat the N₂O result as a first estimate in that case.
 
---
 
## 5. Turning flow area into drill sizes
 
Both injectors end up with a required total flow area $A$. Two common layouts are supported (`Coaxial_*` and `Impinging_*` functions).
 
**Coaxial (single hole).** One circular orifice carries all the flow:
 
$$
d = 2\sqrt{\frac{A}{\pi}}
$$
 
**Impinging (multiple holes).** The area is split among $n$ identical holes whose jets collide:
 
$$
d_{each} = 2\sqrt{\frac{A/n}{\pi}}
$$
 
**Coaxial N₂O with an IPA needle.** In the coaxial layout the IPA flows down a needle running through the middle of the N₂O hole, so the N₂O flows through the *annulus* between the needle and the hole wall. The drilled hole must therefore be larger than the flow area alone would suggest. It has to provide the N₂O flow area *plus* the area the needle takes up (`Coaxial_Drill_Diameter_Nitrous()`):
 
$$
d_{drill} = 2\sqrt{\frac{A_{N_2O} + \pi\,(d_{needle}/2)^2}{\pi}}
$$
 
This is why the order of design matters: size the IPA orifice first, choose a needle with a suitable inner diameter (outer diameter from the manufacturer's spec), and then feed that OD into the N₂O calculation.
 
**Sensitivity to $C_d$.** $C_d$ is the least certain input in the injector calculation, so the plots show the drill diameter against a sweep of $C_d$. Since $d \propto 1/\sqrt{C_d}$, a drop in $C_d$ from 0.8 to 0.6 makes the hole about 15% larger. Rough guide values: sharp-edged orifice ≈ 0.6–0.65, drilled hole with L/D of 2–4 ≈ 0.7–0.8, chamfered or rounded entry ≈ 0.85–0.95. Ultimately, flow testing (water flow tests for the IPA side) is the way to pin down real values.
 
---
 
## 6. A worked example
 
Using the values in `UserInputs.py` (50 N, $P_c$ = 8 bar, ambient 1 atm, O/F = 2.5, $T_c$ = 2500 K, $\gamma$ = 1.2, $M$ = 22):
 
| Step | Result |
|---|---|
| Gas constant $R = 8314.46/22$ | 377.9 J/(kg·K) |
| $c^*_{ideal}$ → real (×0.95) | 1499 → 1424 m/s |
| Expansion ratio $\varepsilon$ (ideal expansion to 1 atm) | 1.94 |
| $C_F$ ideal → real (×0.95) | 1.213 → 1.152 |
| Throat area $A_t = F/(C_F P_c)$ | 54.3 mm² ($D_t$ = 8.31 mm) |
| Exit diameter | 11.58 mm |
| Total mass flow $\dot m = P_cA_t/c^*$ | 30.5 g/s |
| Split at O/F = 2.5 | 21.8 g/s N₂O, 8.7 g/s IPA |
| $I_{sp}$ | 167 s |
| Chamber volume $L^*A_t$ (L* = 1 m) | 54.3 cm³ |
| Stay time | 1.5 ms |
 
Those two propellant flows are then what the IPA and N₂O sections turn into hole sizes.
 
---
 
## 7. Assumptions and limitations
 
- **Hot side is 1-D, isentropic, frozen-composition ideal gas.** Real nozzles have boundary layers, finite-rate chemistry and non-uniform exit flow, which is what the efficiency factors approximately account for.
- **Combustion gas properties must be supplied by you.** $T_c$, $\gamma$ and $M$ depend on propellants, $P_c$ and O/F, and should come from CEA, RPA or `rocketcea`. They are *not* the same as the N₂O properties used on the injector side.
- **Mass flow consistency.** If you type the IPA and N₂O flows in by hand rather than letting the engine sizing set them, make sure they add up to the $\dot m$ your thrust target actually requires. In the GUI, enabling engine sizing overwrites the manual flow fields for you.
- **N₂O is modelled as an ideal gas.** See §4.3; liquid or two-phase N₂O needs a different model.
- **IPA is modelled as incompressible and single-phase,** with no allowance for cavitation or flashing.
- **$C_d$ values are inputs.** The tool does not predict them.
- **Geometry is first-order.** The chamber is a cylinder plus cone with no throat radius of curvature, injector-face detail, cooling or structural sizing.
- **Divergence loss $\lambda$ is reported but not applied automatically** (see §1.6).
---
 
## 8. References

- NASA CEA (Chemical Equilibrium with Applications).