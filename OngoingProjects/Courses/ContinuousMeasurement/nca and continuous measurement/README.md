# NCA and Continuous Measurement

A framing note for this folder. The organizing question is one sentence: given the stochastic differential equation of a continuously monitored quantum system, what can noncommutative operator algebra tell you *exactly*, before any grid or trajectory, and how? The files here work that question in general; the two parent essays work it in full for one system.

## What is here

- [exact-relations-from-monitored-sdes.md](exact-relations-from-monitored-sdes.md): the guide. Starting from an SDE, the exact relations you can derive (closed forms, the ODE system a hierarchy closes into, the deterministic covariance Riccati), the questions they answer, and how `NonCommutativeAlgebra` computes and automates them from the equation alone. Every cell is kernel-verified.
- [nca-exact-relations-audit-plan.md](nca-exact-relations-audit-plan.md): the working plan and its adversarial audit. Literature grounding, the honest positioning (this is validation, not a new method; the conditioning term matters more than the drift), and an append-only revision log. A reference, not exposition.

The worked instance, the free particle under position monitoring, lives one level up:

- [../simulating-continuous-position-measurement.md](../simulating-continuous-position-measurement.md): its Part II is the canonical demonstration, one commutator handed to the algebra, out come the heating law, the Gaussian closure, the deterministic Riccati, and the minimum-uncertainty steady state.
- [../simulating-the-measurement-record.md](../simulating-the-measurement-record.md): the numerics-and-experiment companion, which uses those exact relations as benchmarks and identifies the deterministic Riccati with the Kalman filter that levitated-optomechanics experiments run in real time.

## The fundamental questions

Ordered weakest premise first. Each is paired with what answers it and why it matters.

1. **Is the exact content of a monitored SDE algebraic at all?** The premise everything rests on: every observable evolves by $\frac{d}{dt}\langle\hat O\rangle=\langle\mathcal L^\dagger\hat O\rangle$ in the mean, and by a bracket-built stochastic equation when conditioned, so the content is commutators and anticommutators, which `NonCommutativeAlgebra` carries from a single relation such as $[\hat x,\hat p]=i\hbar$. *Why it matters:* if this fails, nothing else stands; that it holds is exactly what makes the whole approach pre-numerical.

2. **When does the moment hierarchy close?** The pivotal question. Closure turns an infinite hierarchy into a finite linear system $\dot{\vec m}=M\vec m+\vec b$, from which the exact transient $e^{Mt}$, the exact rates $\operatorname{spec}M$, and the steady state $-M^{-1}\vec b$ all follow in closed form, and lets you skip Monte Carlo entirely. It closes for quadratic $\hat H$ and linear $\hat A$ (the Gaussian case) and breaks for nonlinear $\hat H$. *Why it matters:* it is the difference between a closed-form answer and a simulation, and deciding it from the equation is the single most useful thing the algebra does.

3. **Which relations test the measurement, and which only the drift?** The averaged relations (the heating slope $\lambda\hbar^2$, the closed moments) hold with the detector unplugged, so a simulator can reproduce every one of them with the conditioning term entirely wrong. Only the conditional-covariance Riccati, its efficiency dependence ($\sigma_x\propto\eta^{-3/8}$, steady purity $\sqrt\eta$), and the martingale structure probe the $dW$ term. *Why it matters:* a benchmark set that does not separate these certifies the easy half and misses the errors that actually occur in a continuous-measurement code.

4. **Given only the SDE, how do noncommutative-algebra features obtain and automate these?** Declare the algebra the monitored $\hat A$ generates (canonical, $\mathfrak{su}(2)$, Clifford), build $\mathcal L^\dagger$ from $\hat H$ and $\hat A$, normal-order with `NonCommutativeExpand`, apply it to a candidate operator basis, and read closure off the result, assembling $M$ or the Riccati. The mechanization is a single `applyGenerator` map over the basis. *Why it matters:* the pipeline runs from the equation, not from a hand derivation; its honest limit is that closure is searched over monomial spans and noncommutative bases need not terminate.

5. **Why should anyone care?** These are pre-simulation benchmarks that catch the errors a plausible-looking trajectory hides: a wrong sign or Itô-versus-Stratonovich convention on the conditioning term, a drifting norm, a grid too coarse to hold $\sqrt{\lambda\hbar^2 t}$ from above or resolve $\sigma_x$ from below. And the payoff is not academic: the deterministic-Riccati Kalman filter these relations describe is the exact algorithm that levitated-nanoparticle and membrane optomechanics run on hardware in real time.

## The open frontier

One question the corpus does not yet answer: does any of this survive a second system? Every exact relation demonstrated so far is the free particle under position monitoring. It becomes a method, rather than one worked example, only when it is carried to a monitored spin or a quadrature with a nontrivial quadratic $\hat H$, and to at least one deliberately non-closing case where the algebra has to report that no finite benchmark exists.
