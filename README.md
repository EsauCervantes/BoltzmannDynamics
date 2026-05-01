# Boltzmann Equation Solver for Cannibal Dark Sectors

This repository contains a **numerical solver for coupled Boltzmann equations** arising in self-interacting dark matter models produced via the freeze-in mechanism, as studied in

> **E. Cervantes**, *Freezing-in cannibal dark sectors* (2024).  
> [arXiv:2407.12104](https://arxiv.org/abs/2407.12104)

The code is written in **Wolfram Mathematica** and is designed to be fully reproducible.

---

In our work we consider two scenarios: a dark sector composed only by a singlet unstable dark matter candidate $\phi$ (cBE.wl) and a dark sector composed of:

- A **dark matter particle** $S$ with self number changing reactions of the form $3\leftrightarrow 2$,
- A **mediator** $\phi$,

found in the cBE with mediator.wl file. The code contains production from Higgs annihilation/decay as well as production from electroweak gauge bosons. For details see the paper.
The code includes:

- **Freeze-in production** of $S$ and/or $\phi$ from the SM bath (Higgs Portal)
- **Self-interacting processes** in the dark sector ($\phi \phi \phi \leftrightarrow \phi \phi $)
- **Hidden sector temperature evolution and co-moving number density**

---

## Features

### Coupled Boltzmann equations

The Boltzmann equation is an integro differential equation, whose solution is the phase-space distribution function of an ensemble of particles. Solving it in full generality is computationally expensive and often not really necessary. For instance, when the system is in thermal equilibrium, tracking the temperature and comoving number of particles is sufficient, as the system can be described by a Bose-Einstein or Fermi-Dirac distribution. This code does exactly this, and it is optimized to deal with stiffness during freeze-out. The solver also includes the $3 \leftrightarrow 2$ collision integral tabulated as a function of $m/T$.

### Hidden sector temperature tracking

The solver tracks the evolution of the **dark-sector temperature** $T'$, allowing for:

- Heating due to $3\to 2$
- Entropy exchange between the dark matter and the mediator

### Stiff ODE handling

The coupled Boltzmann equations are often **stiff** during cannibal phases.  
We use Mathematica’s ODE solvers with controlled precision and step sizes to obtain stable solutions. From the late-time asymptotic value of Y, the code computes the dark matter relic abundance

---

## Repository Structure

```text
Boltzmann-equations/
│
├── cBE.nb and cBE with mediator.wl      # Main Mathematica notebook with the coupled ODE solver
├── parameters/             # (Optional) Benchmark parameter files
├── plots/                  # (Optional) Generated figures: Y(x), T'/T, relic density, etc.
└── README.md               # This file
