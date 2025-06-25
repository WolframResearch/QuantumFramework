# Quantum Stochastic Master Equation Simulation in Mathematica

## Overview

This Mathematica notebook implements an efficient numerical scheme for simulating the **Quantum Stochastic Master Equation (SME)** based on the method proposed in [Macieszczak *et al.*, arXiv:1410.5345](https://arxiv.org/pdf/1410.5345). The SME describes the evolution of an open quantum system undergoing continuous monitoring and environmental decoherence, extending the Lindblad master equation with stochastic noise terms.

Instead of relying on Mathematica's built-in stochastic solvers (e.g., `NDSolve` or `ItoProcess`), this notebook uses a custom, manual evolution to preserve Bloch sphere dynamics and ensure positivity of the density matrix at each time step.

---

## References

* Macieszczak *et al.*, "Efficient numerical simulation of the stochastic master equation," arXiv:1410.5345.
* Prof. Gabriel T. Landi's melt notebook examples.

---

## License

This work is released under the MIT License.")}
