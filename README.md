# Mathematical biology project

This project explores patterned vegetation, tipping points, and the rate of climate change, based on the research paper "Patterned vegetation, tipping points, and the rate of climate change" by Chen, Yuxin & Kolokolnikov, Theodore & Tzou, Justin & Gai, Chunyi

## Overview

Major fluctuations in Earth's climate have occurred throughout history, and life has generally managed to adapt. The critical issue is whether a species can withstand the speed at which these changes occur. This project investigates what happens when the precipitation rate decreases at different speeds, using a modified Klausmeier model.

The model describes plant density (`n`) and water concentration (`w`). The main terms are:
- Plant growth: `w n^2`
- Plant death: `-m n`
- Plant diffusion: `delta * d^2n/dx^2`
- Water dynamics: precipitation (`a`), evaporation (`-w`), absorption by plants (`-w n^2`), and water diffusion.

A key aspect is the timescale separation: water dynamics are fast, plant dynamics are slow. The model is simplified by setting the water timescale parameter `b=0` (we neglect the changes in water level according to the timescale of the water) and slope `c=0`. Most theory is developed for `m=1`, but the code allows for other values.

## Analytical Approach

- The homogeneous steady states are computed, including the desert state and two vegetated equilibria (`E_+`, `E_-`).
- Stability is analyzed using the Jacobian.
- The effect of a time-dependent precipitation rate `a(t) = a_0 - ε t + noise` is studied, where `ε` is the rate of precipitation loss and the noise term models environmental variability as a spatio-temporal Wiener process.
- The system is linearized around the stable equilibrium, and a Fourier decomposition is used to analyze pattern formation.
- The Fokker-Planck equation is used to study the variance of perturbations, leading to analytical predictions for the "take-off" time and the delayed pattern onset value `a_d`.

## Numerical Simulations

Finite difference schemes are used for spatial and temporal discretization. The main files provided are:

#### `chen_et_al.m`
Simulates the 1D model with time-dependent precipitation and noise. Computes the maximum plant density and spread, and estimates the delayed pattern onset numerically.

#### `chen_et_al_analysis.m`
Similar to `chen_et_al.m`, but includes the analytical calculation of the delayed take-off value `a_d` using the Laplace method, and compares it to the numerical result.

#### `chen_et_al_over_eps.m`
Studies the effect of varying the rate `ε` (precipitation loss) on pattern formation and extinction thresholds.

#### `chen_et_al_over_sigma.m`
Studies the effect of varying the noise amplitude `σ_0` on the onset of pattern formation.

#### `chen_et_al_2D.m`
Extends the simulation to 2D, showing the emergence of stripes and spots as precipitation decreases.

#### `Klausmeier_TuringInstability_Soresina.m`
Analyzes the Turing instability and bifurcation structure of the modified Klausmeier model, including the computation of critical values and eigenmodes. Credits goes to Dr Cinzia Soresina, professor of the 2025 course of mathematical biology at the university of Trento.

## How to Run

Open the desired `.m` file in MATLAB and run it. Each script is self-contained and will produce figures illustrating the dynamics, pattern formation, and bifurcation diagrams.

## Notes

- The code uses finite difference discretization with Neumann boundary conditions.
- The noise is implemented as spatially uncorrelated Gaussian white noise.
- The main parameters to vary are `a0`, `a1`, `delta`, `m`, `eps`, and `sigma0`.
- The analytical and numerical results for the delayed pattern onset are compared in `chen_et_al_m1.m`.

## References

- Chen, Yuxin & Kolokolnikov, Theodore & Tzou, Justin & Gai, Chunyi. (2015). Patterned vegetation, tipping points, and the rate of climate change. European Journal of Applied Mathematics. -1. 1-14. 10.1017/S0956792515000261. 
- Klausmeier, C. A. (1999) Regular and irregular patterns in semiarid vegetation. Science 284(5421), 1826–1828

