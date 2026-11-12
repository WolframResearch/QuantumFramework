# Quantum calculations QF supports beyond Qiskit's supplied capabilities

Research date: 3 October 2026. QF source: paclet 2.1.1, canonical checkout at `bf5a1ee0e251b4097d5a8b2c81e63b6e99928ce3`.

**Symbolic computation throughout circuits, states, channels and dynamics is the starting point of this comparison. The research question is which scientifically useful results that support makes possible.** The examples below were selected from implemented capabilities and fresh literature searches, independently of the earlier question draft.

The comparison concerns operations supplied by QF and Qiskit. An external Python package or a user-written symbolic implementation is recorded as an alternative; its existence does not turn it into a Qiskit capability. Actual Qiskit extensions do count, with their scope identified. In particular, [qiskit-symb](https://github.com/SimoneGasperini/qiskit-symb) supplies symbolic circuit operators, statevectors and density matrices. Thus “Qiskit has no symbolic capabilities” would be an incorrect ecosystem-wide claim.

**1. Derive a logical error channel as a function of coherent physical errors.**

Take a fixed small quantum code, its encoding isometry \(V\), syndrome projectors \(P_s\), recovery operators \(R_s\), and physical noise Kraus operators \(E_\mu(\boldsymbol\theta)\). Calculate the logical Kraus operators

\[
K_{s\mu}(\boldsymbol\theta)
=V^\dagger R_sP_sE_\mu(\boldsymbol\theta)V,
\qquad
\mathcal L_{\boldsymbol\theta}(\rho)
=\sum_{s,\mu}K_{s\mu}\rho K_{s\mu}^\dagger.
\]

For a recovery whose output lies in the code space, this gives the recovered logical channel. The desired output is its exact parameter-dependent process matrix, followed by the region where encoded storage has higher entanglement fidelity than one unencoded qubit exposed to the same single-qubit noise for the same interval. Assume ideal encoding, syndrome measurement and recovery. A second calculation compares the logical channel with that obtained after Pauli twirling the physical noise.

QF provides the symbolic operators and channels needed for this composition; the code and decoder must still be specified. Its [Kraus constructor](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumChannel/QuantumChannel.m:38) accepts operator lists, and its [channel transformations](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumChannel/Properties.m:83) preserve symbolic simplification. Wolfram Language then supplies the algebraic simplification and parameter-region solving. This is a source-supported workflow, not a claim that QF contains a named five-qubit-code performance solver or that this full application was executed here.

Qiskit can simulate the corresponding error-correction experiment and compose channels at numerical parameter values. Its [Kraus implementation](https://github.com/Qiskit/qiskit/blob/main/qiskit/quantum_info/operators/channel/kraus.py) uses numerical complex arrays. qiskit-symb supplies symbolic circuit operators; combined with user-written SymPy syndrome-branch contractions, it is a real alternative for coherent-error examples. The supplied-operation gap is native symbolic channel algebra. The comparison does not establish that the extended ecosystem cannot derive this result with that additional construction.

This is an established scientific calculation: [Rahn, Doherty and Mabuchi](https://doi.org/10.1103/PhysRevA.66.032304) derive coding maps and concatenation thresholds; [Huang, Doherty and Flammia](https://doi.org/10.1103/PhysRevA.99.022313) analyze coherent logical errors.

A concrete target is the Steane code under independent identical mixed dephasing and coherent Z errors, with the paper's Z-only minimum-distance decoder. Write the physical channel as \(\mathcal N(\rho)=(1-x)\rho+xZ\rho Z-iy[Z,\rho]\). For pure rotations \(e^{-i\theta Z}\), \(x=\sin^2\theta\) and \(y=\sin\theta\cos\theta\). The exact logical polynomial map includes the coherence-dependent contribution \(42(1-2x)^3y^4\) to \(x'\). Consequently, at small angle, \(x'\sim63\theta^4\) and \(y'\sim14\theta^3\), whereas twirling the physical noise first predicts \(x'\sim21\theta^4\), \(y'=0\). Thus it misses a factor of three in leading entanglement infidelity and changes the predicted diamond-distance scaling from cubic to quartic. These are a literature result and its immediate twirling comparison, not results produced by a QF run here. [Exact map, Eqs. (67)–(68)](https://arxiv.org/html/1805.08227#S6.E67).

**2. Find the interference condition for photon antibunching in coupled nonlinear cavities.**

Use two coupled, driven, lossy cavities and retain symbolic detuning, coupling, Kerr nonlinearity and loss. In a controlled weak-drive expansion, derive the driven cavity's two-photon amplitude and solve for its cancellation. Then use numerical dynamics to determine the drive range with \(g_a^{(2)}(0)=\langle a^{\dagger2}a^2\rangle/\langle a^\dagger a\rangle^2<1\), and the accompanying mean photon number. This makes the antibunching and brightness tradeoff explicit.

The useful output is a relation among physical parameters explaining *why* two excitation paths cancel. [Bamba and colleagues](https://doi.org/10.1103/PhysRevA.83.021802) establish this interference mechanism for weakly nonlinear photonic molecules. The cancellation belongs to the weak-drive approximation; it must not be described as an exact solution of the unrestricted driven oscillator model.

QF supplies [symbolic bosonic expressions and their finite representations](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/SecondQuantization/Utils.m:125), together with [symbolic and numerical evolution](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumEvolution.m:153). Its equation-returning path also allows symbolic stationary equations to be extracted and solved with Wolfram Language. [Qiskit Dynamics](https://qiskit-community.github.io/qiskit-dynamics/stubs/qiskit_dynamics.solvers.Solver.html) supplies numerical Hamiltonian and Lindblad evolution, including differentiable numerical workflows. The analytic cancellation condition requires a separate symbolic calculation. This is a promising concrete application of QF's symbolic dynamics; the paper-specific derivation remains to be implemented.

**3. Identify redundant directions and singularities in a variational circuit.**

For a normalized pure variational state, compute the Fubini–Study metric as a symbolic matrix, then determine its generic rank, local null directions and parameter values where its rank changes. This distinguishes an exactly null tangent direction from a small but nonzero numerical eigenvalue. On a constant-rank region, the analysis can guide a local reduced parameterization for optimization or multiparameter sensing; it does not automatically establish global redundancy.

QF's [`FubiniStudyMetricTensor`](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumOptimization.m:363) explicitly differentiates the symbolic state and constructs the parameter-dependent matrix. Null-space and rank-condition analysis use Wolfram Language, with real-parameter and normalization assumptions stated. QF's natural-gradient routine does not automatically remove singular directions.

Qiskit Algorithms' [`ReverseQGT`](https://qiskit-community.github.io/qiskit-algorithms/stubs/qiskit_algorithms.gradients.ReverseQGT.html) evaluates the tensor at supplied parameter values. qiskit-symb can provide a symbolic state from which a user could construct the metric and analyze it with SymPy; that is the concrete alternative, not an already supplied symbolic QGT operation. [Quantum circuit geometry research](https://doi.org/10.1103/PRXQuantum.2.040309) connects metric rank with effective circuit dimension.

**4. Calculate oscillator transition amplitudes without exponentiating a truncated ladder operator.**

For displaced harmonic potential surfaces of equal curvature, compute Franck–Condon transition probabilities as functions of displacement. QF directly supplies

\[
\langle m|D(\alpha)|n\rangle
=e^{-|\alpha|^2/2}\sqrt{\frac{n!}{m!}}\,
\alpha^{m-n}L_n^{(m-n)}(|\alpha|^2),
\qquad m\ge n,
\]

through [`BosonicMatrixElement`](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/SecondQuantization/BosonicAlgebra.m:232). The same implementation supplies squeezing matrix elements, relevant to oscillator frequency quenches. These expressions can remain symbolic or be evaluated numerically.

This differs mathematically from the inspected [C2QA bosonic Qiskit extension](https://raw.githubusercontent.com/C2QA/bosonic-qiskit/main/src/bosonic_qiskit/operators.py), which builds displacement by exponentiating truncated creation and annihilation matrices. In general, projecting an infinite-dimensional exponential is not the same as exponentiating the projected generator. Qiskit Nature's polynomial normal ordering does not supply this displacement kernel.

A thermal spectrum assembled from these exact transition elements still needs control of the tail of the state sum. The claim is exact individual matrix elements, not an automatically cutoff-free spectrum or a general multimode vibronic solver.

**5. Infer posterior uncertainty about entanglement from measurement counts.**

Given finite two-qubit tomography counts, sample physical density matrices conditioned on those counts. Calculate a posterior distribution for Bell-state overlap or negativity, and a posterior probability of entanglement. For example, the posterior probability that \(\langle\Phi^+|\rho|\Phi^+\rangle>1/2\) quantifies support for a sufficient entanglement witness.

QF's [`QuantumStateEstimate`](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumStateEstimation.m:285) returns a `"BayesianSampler"`, implemented with [Bures-based proposals and likelihood weighting](/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumStateEstimation.m:118). It supplies the sampling machinery; meaningful credible intervals still require assessing mixing and prior dependence. Samples also need the correct two-qubit subsystem interpretation before evaluating bipartite quantities.

[Qiskit Experiments 0.14.2 tomography](https://qiskit-community.github.io/qiskit-experiments/stubs/qiskit_experiments.library.tomography.TomographyAnalysis.html) supplies linear-inversion and least-squares fitters, plus bootstrap uncertainty for target fidelity. Those operations do not produce a posterior distribution over physical density matrices. A custom-fitter hook and random Bures states do not supply the missing posterior sampler. External packages such as [QInfer](https://pythonhosted.org/QInfer/guide/tomography.html) provide a separate Bayesian route. [Lougovski and van Enk](https://doi.org/10.1103/PhysRevA.80.052324) demonstrate the scientific use of posterior entanglement probabilities.

The fifth example is an additional statistical capability; it does not replace the symbolic examples.

**Research scope and verification.** Two new searches addressed [scientific results requiring symbolic and numerical operations together](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Quantum%20results%20requiring%20symbolic%20and%20numerical%20operations%20together) and [physical conclusions from Bayesian state posteriors](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Physical%20conclusions%20from%20Bayesian%20quantum%20state%20posteriors). They returned 297 and 37 records respectively; the search summaries and first 20 relevance entries of each were inspected, not all papers in full. Software conclusions rest on targeted live QF source and primary Qiskit documentation/source. Extension versions are inspected independently, not asserted to form one compatible installation.

Targeted primary-text review covered Bamba's model, weak-drive equations and limitations, and Huang's channel convention and Steane polynomial map. Neither paper was read cover to cover. The metric, oscillator matrix-element and posterior-sampling operations are directly present in QF source. The QEC and cavity applications combine existing operations with model-specific algebra; they are research directions supported by the implementation, not completed reproductions. Bounded kernel probes encountered dependency-loading failures and produced no usable application result. No end-to-end numerical comparison is claimed.

## Undermind research references

Saved searches for this comparison, including the earlier exploratory route. Their literature results do not establish software support by themselves.

- [Quantum results requiring symbolic and numerical operations together](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Quantum%20results%20requiring%20symbolic%20and%20numerical%20operations%20together) — Direct source search for the revised comparison: 297 records; the original review inspected the summary and first 20 relevance entries.
- [Physical conclusions from Bayesian quantum state posteriors](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Physical%20conclusions%20from%20Bayesian%20quantum%20state%20posteriors) — Direct source search for posterior-based physical conclusions: 37 records; the original review inspected the summary and first 20 relevance entries.
- [Quantum numerical applications with essential nonstandard representations](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Quantum%20numerical%20applications%20with%20essential%20nonstandard%20representations) — Earlier exploratory search, retained for reuse: 144 records on representation-sensitive numerical applications. It does not establish a Qiskit capability gap.

Workspace: [Roadmap for Wolfram quantum](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3); workspace ID: `bcedc042-7a0d-4f3d-b4ca-4e1c138571a3`. Links checked on 2026-10-04. For reuse, inspect the saved search goal and results before launching a new search; follow the primary-paper links for evidence.
