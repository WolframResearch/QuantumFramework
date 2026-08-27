---
Template: Default
Title: "Backbone: A Cavity Cat, and Why It Decoheres Almost at Once"
Author: Mads Bahrami
---

# Backbone: A Cavity Cat, and Why It Decoheres Almost at Once

This page is the statement layer of the essay: a graph of assumptions, choices, and derived statements, each given in words and formulas only. Every node names its parents, so the dependency graph is recoverable from the page, and the layers order the reading from the setup to the headline and then to what holds when an assumption is weakened. Two conventions hold throughout. The equation of motion is never a root: it enters as the first derived node, translated from the physical setup by standard theory and cited, because the physically meaningful weakenings are acts on the setup (warm the bath, let it remember), and only the setup layer makes their effect on the equation canonical. And every assumption must be consumed by at least one statement that is not its own weakening; a root failing that test is not an assumption of the theory but only a direction along which the model can be varied. Node labels: R is a root, a setup assumption; C is a choice, a meaning we fixed rather than a fact; S is a derived statement, with S0 reserved for the equation of motion; H is the headline, the terminal statement; N is a recorded non-assumption; a prime marks the same statement under a weakened root. Nothing here executes; statements marked *checked* are verified computationally in the companion file *two-examples-from-minimal-assumptions.md*, downstream of this page, by routes that share no code with the derivations here.

## The setup, and the choice

**R1 (assumption).** The environment couples to the mode linearly: it attaches to the amplitude $\hat a$ itself, not to any higher power of it. In words: whatever the environment does, it does through the amplitude, one quantum at a time, which is what will make coherent states special and Gaussian overlaps the currency of the whole story. R1 alone implies no loss at all; which way the traffic flows is R2's claim, and how fast the environment forgets is R3's, so the three are independent axes (operator, direction, clock), and the corner models that separate them appear below.

**R2 (assumption).** The environment is at zero temperature: it absorbs quanta and never returns any.

**R3 (assumption).** The environment is memoryless: its correlation time lies far below every other timescale, so the loss proceeds at one constant rate $\gamma$ and the evolution from any instant depends only on the state at that instant.

**R4 (assumption).** The initial state is the even superposition of two coherent states, $|\alpha\rangle + |{-\alpha}\rangle$ normalized, with the lobes well separated: at separation $d = 2\alpha$ the overlap $\left|\langle\alpha|{-\alpha}\rangle\right| = e^{-d^2/2}$ is negligible.

**C1 (choice).** "Coherence" means the weight of the cross term when the state is expanded on its two shrinking lobes. Other reasonable numbers exist, the trace distance to the even mixture of the lobes among them; inside R4's scope they agree with this one up to corrections of the size of the overlap, and outside it they part company, so this choice and R4 travel together. *Checked.*

## The equation, derived from the setup

**S0 (from R1, R2, R3; cited).** Standard open-system theory turns the setup into the reduced dynamics
$$\dot\rho = \gamma\,\mathcal D[\hat a]\rho, \qquad \mathcal D[\hat c]\rho = \hat c\rho\hat c^\dagger - \tfrac12\left(\hat c^\dagger\hat c\,\rho + \rho\,\hat c^\dagger\hat c\right),$$
with each root visible in the equation's shape: the linear coupling R1 makes the jump operator $\hat a$ itself; zero temperature R2 is the absence of any upward term $\mathcal D[\hat a^\dagger]$; memorylessness R3 is the single constant rate. *Cited: Breuer and Petruccione, The Theory of Open Quantum Systems (Oxford 2002), chapter 3; the companion file computes with exactly this equation.*

## Consequences

**S1 (from S0).** On any coherent dyad the loss acts through a single exponent,
$$\mathcal E_t\big(|\alpha\rangle\langle\beta|\big) = \langle\beta|\alpha\rangle^{\,1 - e^{-\gamma t}}\;\big|\alpha\,e^{-\gamma t/2}\big\rangle\big\langle\beta\,e^{-\gamma t/2}\big|.$$
Its diagonal case $\beta = \alpha$ says the loss preserves coherent states and merely shrinks them, $|\alpha\rangle \to |\alpha\,e^{-\gamma t/2}\rangle$, with mean occupation decaying as $e^{-\gamma t}$: the energy clock. Off the diagonal, the dyad is additionally multiplied by the overlap of its two labels raised to the fraction of amplitude lost so far. The entire essay lives in that exponent. *Exact; cited: Haroche and Raimond, Exploring the Quantum (Oxford 2006), and Walls and Milburn, Quantum Optics; checked.*

## The headline

**H (from S1, R4, C1).** Setting $\beta = -\alpha$ in S1,
$$\left|C(t)\right| = \exp\left(-\tfrac{d^2}{2}\left(1 - e^{-\gamma t}\right)\right), \qquad \Gamma_{\mathrm{dec}} = \tfrac12\,\gamma\,d^2 \ \ \text{at short times}, \qquad \frac{\Gamma_{\mathrm{dec}}}{\gamma} = \frac{d^2}{2} = 2\alpha^2.$$
Mechanism, in words: the environment identifies which lobe it is absorbing from, and the leaked fields of the two lobes are distinguishable in proportion to the separation squared, so information drains faster than energy by exactly $d^2/2$. Doubling the cat quadruples its decoherence relative to its dimming; extrapolated to macroscopic separation, this one ratio is why large superpositions are never seen while the objects themselves are still bright. Audit: H rests on R1, R2, R3, R4, and C1, the full root set, with no conjecture in its closure. *Derived in place from S1; checked.*

## When an assumption is weakened

**S1$'$ (from R3 weakened, with R1 and R2 kept; S1).** Let the environment remember: the smallest such environment is one auxiliary mode, coupled at $g$, itself leaking at $\kappa$. The lost fraction $1 - e^{-\gamma t}$ becomes $1 - \left|u(t)\right|^2$ with
$$u'' + \tfrac{\kappa}{2}\,u' + g^2 u = 0, \qquad u(0) = 1, \quad u'(0) = 0,$$
so the decay starts flat ($1 - u^2 \approx g^2 t^2$), partially revives when $\kappa < 4g$ (the mode hands amplitude back before losing it), and recovers R3 in the corner $\kappa \gg g$ with $\gamma_{\mathrm{eff}} = 4g^2/\kappa$. The restructured roots sharpen the attribution: the $d^2$ amplification persists at every instant because it descends from R1, the linearity that keeps the lobes coherent and their overlap Gaussian, while R3 owned only the clock, the exponential profile and the monotonicity of decoherence. *Derived in place from the linear two-mode dynamics with vacuum inputs (the pseudomode analysis; cited: Garraway, Phys. Rev. A 55, 2290 (1997)); in the companion the u equation is decided symbolically on exact input and the channel reduction is checked against the full two-mode simulation.*

**S1$''$ (from R2 weakened, with R1 and R3 kept).** Warm the environment to occupation $\bar n$: a downward channel at rate $\gamma(\bar n + 1)$ and an upward one at $\gamma\bar n$. The short-time decoherence rate is multiplied by $(2\bar n + 1)$, and the diagonal case of S1 fails: a heated lobe is no longer a coherent state, so C1's expansion must be replaced by its mixture-distance rival. Temperature owns the prefactor, R3 the profile, R4 the amplification, R1 the currency. *Cited: Walls and Milburn, Quantum Optics (decoherence in a thermal reservoir); checked.*

**H$'$ (from R4 weakened).** Shrink the pair and the formula in H survives verbatim while the phenomenon vanishes: once $2\alpha^2 \lesssim 1$ the coherence outlives a fair share of the energy. R4 is load-bearing for the headline, never for the law. *Checked.*

**N (a recorded non-assumption).** No "no drive" assumption appears above, on purpose: a free Hamiltonian $\omega\,\hat a^\dagger\hat a$ commutes with the loss, rotating the whole picture rigidly and changing no statement on this page. The boundary of the deletion: any Hamiltonian that keeps coherent states coherent is harmless; a Kerr term is a different problem. *Derived in place (the rotation commutes with the dissipator); checked.*

## The abstract, read off the graph

A single cavity mode surrenders quanta, one at a time through its amplitude, to a memoryless zero-temperature environment. That loss preserves coherent states and drains energy at rate $\gamma$, but on the interference term between two coherent states it acts through one exponent: the overlap of the lobes raised to the fraction of amplitude lost. Since the overlap is $e^{-d^2/2}$ at separation $d$, a superposition of two lobes loses its coherence at rate $\tfrac12\gamma d^2$ while losing energy only at $\gamma$: the interference dies before the light dims, by a factor of the separation squared. Temperature multiplies the rate by $2\bar n + 1$; a bath with memory changes only the clock, a flat onset and partial revivals, never the $d^2$ amplification, which descends from the linearity of the loss rather than from its clock; a free drive changes nothing at all. Only the size assumption is essential to the phenomenon: shrink the pair and the same law describes coherence that outlives the energy.
