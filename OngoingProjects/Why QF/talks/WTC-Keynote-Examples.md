---
Template: Default
Title: Bosonic Algebra, Stabilizers, and Interferometers: Keynote Examples
Author: Mads Bahrami
---

# Bosonic Algebra, Stabilizers, and Interferometers: Keynote Examples

**Bosonic algebra.** Wolfram QuantumFramework now computes exactly in the Weyl algebra of light and oscillators: it disentangles squeezers, beam splitters and even cubic exponentials into normal-ordered closed forms, identities that physicists otherwise derive by hand from commutators, valid at every photon number and never truncated to a finite matrix.

**Stabilizer formalism.** Wolfram QuantumFramework runs Clifford circuits on the stabilizer tableau, $n$ Pauli generators instead of $2^n$ amplitudes, so a thousand-qubit state takes milliseconds and its entanglement comes out as an exact integer from linear algebra over $\mathbb{F}_2$. When a non-Clifford T gate arrives, the same object switches to an exact superposition of stabilizer states instead of failing.

**QuantumInterferometer.** In Wolfram QuantumFramework a beam splitter is an SU(2) rotation of the photon modes, so its photon-counting amplitudes are Wigner D-functions returned symbolically with the photon numbers themselves left as variables. One exact formula covers Hong-Ou-Mandel interference for every $n$, where other photonics software computes one photon count at a time, numerically.

Eighteen short examples from QuantumFramework. Each is one call that returns an exact answer, usually a formula, where a numerical tool would need a truncated Fock space, a $2^n$ state vector, or a table of numbers.

```wl
PacletInstall["https://www.wolfr.am/DevWQCF", ForceVersionInstall -> True]
Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]
```

## Bosonic algebra

Operators of light and oscillators, ordered and disentangled exactly. The field variables are the formal symbols $$ and $$; their creation operators are $^\dagger$ and $^\dagger$. Every result is an operator identity, valid at every photon number, with no Fock-space cutoff.

The ladder operators live in the Weyl algebra, the built-in `WeylAlgebra` with the creators as variables and the annihilators as their derivatives, so $\hat{a}\hat{a}^\dagger=\hat{a}^\dagger\hat{a}+1$, that is $[\hat{a},\hat{a}^\dagger]=1$. One mode and two modes:

```wl
weyl = WeylAlgebra[{{SuperDagger[\[FormalA]]}, {\[FormalA]}}]
weyl2 = WeylAlgebra[{{SuperDagger[\[FormalA]], SuperDagger[\[FormalB]]}, {\[FormalA], \[FormalB]}}];
```

The canonical form of the Weyl algebra, creators left of annihilators, is normal order. Normal-ordering $\hat{a}^2\hat{a}^{\dagger 2}=\hat{a}^{\dagger 2}\hat{a}^2+4\hat{a}^\dagger\hat{a}+2$ in QuantumFramework:

```wl
BosonicNormalOrder[\[FormalA]**\[FormalA]**SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]], weyl]
```

gives exactly the canonical form Wolfram Language computes in the algebra itself:

```wl
NonCommutativeExpand[\[FormalA]**\[FormalA]**SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]], weyl]
```

Every call below takes the algebra it works in as its second argument.

### 1. Squeezing, factored

The squeeze operator $\hat{S}(\xi)=\exp\!\big[\tfrac12(\xi^*\hat{a}^2-\xi\,\hat{a}^{\dagger 2})\big]$, $\xi=re^{i\theta}$. One call returns its exact factorization into creation, number and annihilation parts, the identity squeezed-light calculations rest on, for every $r$ and $\theta$:

$$\hat{S}(\xi)=\Big(e^{-\frac12 e^{i\theta}\tanh(r)\,\hat{a}^{\dagger 2}}\Big)\,\Big((\operatorname{sech}(r))^{\hat{a}^\dagger \hat{a}+\frac12}\Big)\,\Big(e^{\frac12 e^{-i\theta}\tanh(r)\,\hat{a}^2}\Big)$$

```wl
BosonicExpOrder[Exp[1/2 r E^(-I θ) \[FormalA]**\[FormalA] - 1/2 r E^(I θ) SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]], weyl, Assumptions -> r > 0]
```

### 2. Two-mode squeezing

The two-mode squeezer $\hat{S}_2(\xi)=\exp\!\big(\xi^*\hat{a}\hat{b}-\xi\,\hat{a}^\dagger \hat{b}^\dagger\big)$, the source of entangled photon pairs, disentangled just as exactly:

$$\hat{S}_2(\xi)=\operatorname{sech}(r)\,\Big(e^{-e^{i\theta}\tanh(r)\,\hat{a}^\dagger \hat{b}^\dagger}\Big)\,\Big((\operatorname{sech}(r))^{\hat{a}^\dagger \hat{a}+\hat{b}^\dagger \hat{b}}\Big)\,\Big(e^{e^{-i\theta}\tanh(r)\,\hat{a}\hat{b}}\Big)$$

```wl
BosonicExpOrder[Exp[r E^(-I θ) \[FormalA]**\[FormalB] - r E^(I θ) SuperDagger[\[FormalA]]**SuperDagger[\[FormalB]]], weyl2, Assumptions -> r > 0]
```

### 3. A beam splitter between two modes

Mode mixing $\exp(\alpha\,\hat{a}\hat{b}^\dagger+\beta\,\hat{a}^\dagger \hat{b})$ with arbitrary complex $\alpha,\beta$ and $\Omega=\sqrt{\alpha\beta}$: one exact operator identity covers every photon number:

$$e^{\alpha \hat{a}\hat{b}^\dagger+\beta \hat{a}^\dagger \hat{b}}=\Big(e^{\frac{\beta\tanh\Omega}{\Omega}\hat{a}^\dagger \hat{b}}\Big)\,\Big((\cosh\Omega)^{\hat{b}^\dagger \hat{b}-\hat{a}^\dagger \hat{a}}\Big)\,\Big(e^{\frac{\alpha\tanh\Omega}{\Omega}\hat{a}\hat{b}^\dagger}\Big)$$

```wl
Simplify[BosonicExpOrder[Exp[α \[FormalA]**SuperDagger[\[FormalB]] + β SuperDagger[\[FormalA]]**\[FormalB]], weyl2]]
```

### 4. A nonlinear exponential

An exponential cubic in the ladder operators, $e^{k\,\hat{a}^\dagger \hat{a} \hat{a}^\dagger}$, outside the quadratic cases that standard disentangling identities cover, still has an exact normal-ordered form:

$$e^{k\,\hat{a}^\dagger \hat{a} \hat{a}^\dagger}=\;:\!\Big(\frac{1}{1-k \hat{a}^\dagger}\Big)\,\Big(\exp\!\Big[\Big(\frac{1}{1-k \hat{a}^\dagger}-1\Big)\hat{a}^\dagger \hat{a}\Big]\Big)\!:$$

```wl
BosonicExpOrder[Exp[k SuperDagger[\[FormalA]]**\[FormalA]**SuperDagger[\[FormalA]]], weyl]
```

Its normal-ordered series in $k$, term by term:

```wl
Collect[BosonicExpOrder[Exp[k SuperDagger[\[FormalA]]**\[FormalA]**SuperDagger[\[FormalA]]], weyl]["Series", k, 3], k]
```

### 5. Reordering a product

Reordering $e^{k \hat{a}^2}e^{m \hat{a}^{\dagger 2}}$ into normal order, a page of commutator algebra by hand, in one call:

$$\Big(e^{k \hat{a}^2}\Big)\,\Big(e^{m \hat{a}^{\dagger 2}}\Big)=\frac{1}{\sqrt{1-4km}}\,\Big(e^{\frac{m}{1-4km}\hat{a}^{\dagger 2}}\Big)\,\Big((1-4km)^{-\hat{a}^\dagger \hat{a}}\Big)\,\Big(e^{\frac{k}{1-4km}\hat{a}^2}\Big)$$

```wl
BosonicExpOrder[Exp[k \[FormalA]**\[FormalA]]**Exp[m SuperDagger[\[FormalA]]**SuperDagger[\[FormalA]]], weyl]
```

and a displacement into anti-normal order, $e^{\lambda \hat{a}+\beta \hat{a}^\dagger}=e^{-\beta\lambda/2}\,\big(e^{\lambda \hat{a}}\big)\,\big(e^{\beta \hat{a}^\dagger}\big)$, the ordering behind the Husimi $Q$ function:

```wl
BosonicExpOrder[Exp[λ \[FormalA] + β SuperDagger[\[FormalA]]], weyl, "Ordering" -> "Antinormal"]
```

## Stabilizer formalism

A Clifford state stored as the group of Pauli operators that fix it, $\hat{P}|\psi\rangle=|\psi\rangle$: $n$ generators instead of $2^n$ amplitudes, so the cost grows polynomially in the number of qubits, not exponentially.

### 6. Same circuit, different engine

A five-qubit cluster-state circuit: Hadamards, then controlled-Z gates down the chain:

```wl
cluster5 = QuantumCircuitOperator[Join[Thread["H" -> Range[5]], {"CZ" -> {1, 2}, "CZ" -> {2, 3}, "CZ" -> {3, 4}, "CZ" -> {4, 5}}]];
cluster5["Diagram"]
```

The same circuit on the stabilizer engine, one option away, returns a tableau of generators instead of $2^5$ amplitudes:

```wl
stab5 = cluster5[Method -> "Stabilizer"]
```

Its generators are $\hat{K}_v=\hat{X}_v\prod_{w\sim v}\hat{Z}_w$, one per site:

```wl
stab5["PauliForm"]
```

Expectation values come straight from the generators, exactly and with no state vector: $\langle\hat{K}_1\rangle=1$, while $\hat{X}_1\hat{X}_2$, which anticommutes with $\hat{K}_1$, averages to zero:

```wl
<|"XZIII" -> stab5["Expectation", "XZIII"], "XXIII" -> stab5["Expectation", "XXIII"]|>
```

### 7. A thousand qubits, exact entanglement

Twenty thousand random Clifford gates on 1000 qubits, a state with $2^{1000}$ amplitudes, then its entanglement entropy as an exact integer, read off the generators restricted to half the chain $A$:

$$S(A)=\operatorname{rank}_{\mathbb{F}_2}\big(G_A\big)-|A|$$

A seeded random stream of Hadamard, CNOT and phase gates:

```wl
nq = 1000;
SeedRandom[1];
gates = Flatten@Table[{"H" -> RandomInteger[{1, nq}], "CNOT" -> RandomSample[Range[nq], 2], "S" -> RandomInteger[{1, nq}]}, {6666}];
Length@gates
```

Applied to the tableau in one compiled pass, timed:

```wl
(big = PauliStabilizer[nq]["ApplyCircuit", gates]) // AbsoluteTiming
```

The half-chain entropy, an exact integer rather than a sampled estimate:

```wl
big["Entropy", Range[nq/2]]
```

### 8. An 800-qubit GHZ state

The state $|\mathrm{GHZ}\rangle=\big(|0\rangle^{\otimes 800}+|1\rangle^{\otimes 800}\big)/\sqrt{2}$, built straight from its circuit name; a dense simulator would need $2^{800}$ amplitudes:

```wl
g = PauliStabilizer["GHZ"[800]]
```

It is fixed by $\hat{X}^{\otimes 800}$:

```wl
g["Expectation", StringJoin@ConstantArray["X", 800]]
```

and by $\hat{Z}_1\hat{Z}_{800}$, which links the two ends of the chain:

```wl
g["Expectation", "Z" <> StringJoin@ConstantArray["I", 798] <> "Z"]
```

### 9. Circuit, tableau, state, circuit

A three-qubit cluster circuit:

```wl
cluster3 = QuantumCircuitOperator[{"H" -> Range[3], "CZ", "CZ" -> {2, 3}}];
cluster3["Diagram"]
```

The stabilizer engine turns it into a tableau:

```wl
stab3 = cluster3[Method -> "Stabilizer"]
```

Expanded, the tableau gives exactly, not approximately, the state the default engine computes, $|\psi\rangle=\frac{1}{2\sqrt{2}}\sum_{x\in\{0,1\}^3}(-1)^{x_1x_2+x_2x_3}\,|x_1x_2x_3\rangle$:

```wl
QuantumState[stab3]["StateVector"] === cluster3[]["StateVector"]
```

And the tableau synthesizes a circuit that prepares it, closing the loop from circuit to tableau to state to circuit:

```wl
stab3["Circuit"]["Diagram"]
```

### 10. Beyond Clifford

A T gate leaves the stabilizer group. Instead of failing or approximating, the engine switches to a superposition of stabilizer states:

```wl
fr = QuantumCircuitOperator[{"H", "T"}][Method -> "Stabilizer"]
```

Still exact: $\hat{T}\hat{H}|0\rangle=\big(|0\rangle+e^{i\pi/4}|1\rangle\big)/\sqrt{2}$, and $e^{i\pi/4}/\sqrt{2}=\tfrac12+\tfrac{i}{2}$:

```wl
fr["State"] // FullSimplify // TraditionalForm
```

A ring cluster state, built directly from a graph:

```wl
qc = QuantumCircuitOperator["Graph"[CycleGraph[5]]];
qc["Diagram"]
```

Its tableau: the generators are $\hat{K}_v=\hat{Z}_{v-1}\hat{X}_v\hat{Z}_{v+1}$ around the ring, so the $Z$ block is the adjacency matrix of the cycle:

```wl
PauliStabilizer[qc]["TableauForm"]
```

## QuantumInterferometer

Linear optics with exact photon-counting amplitudes: permanents, hafnians and Wigner D-functions in closed form, with angles and photon numbers that can stay symbolic. Below, $D^{j}_{m m'}(\alpha,\beta,\gamma)$ is the Wigner D-function, `WignerD[{j, m, m'}, α, β, γ]`.

### 11. A beam splitter, any photon numbers

A beam splitter with symbolic angles, as one object:

```wl
bs = QuantumInterferometer[{"BS"[θ, φ]}]
```

```wl
bs["Diagram"]
```

The amplitude for $n_1,n_2$ photons in and $t_1,t_2$ out, every number symbolic, so one result covers every photon-number case, is a conservation delta times one Wigner D-function, $A=\delta_{n_1+n_2,\,t_1+t_2}\,D^{\frac{n_1+n_2}{2}}_{\frac{n_1-n_2}{2},\;t_1-\frac{n_1+n_2}{2}}(\varphi,2\theta,-\varphi)$:

```wl
bs["Amplitude", {n1, n2} -> {t1, t2}]
```

The same result in traditional notation:

```wl
% // TraditionalForm
```

Photon numbers can be expressions; $2n$ photons in one port, $n$ out of each: $A(2n,0\to n,n)=D^{\,n}_{n,\,0}(\varphi,2\theta,-\varphi)$:

```wl
bs["Amplitude", {2 n, 0} -> {n, n}]
```

```wl
% // TraditionalForm
```

Laser light stays laser light: a coherent state goes to a coherent state, $|\alpha,\beta\rangle\to|\alpha\cos\theta-e^{-i\varphi}\beta\sin\theta,\;\beta\cos\theta+e^{i\varphi}\alpha\sin\theta\rangle$:

```wl
bs["Coherent"[{α, β}]]
```

Its overlap with any coherent state $|\gamma,\delta\rangle$ is a closed-form Gaussian in their distance, not a sum over truncated Fock states, $P=\exp\!\big(-|\gamma-\alpha\cos\theta+e^{-i\varphi}\beta\sin\theta|^2-|\delta-\beta\cos\theta-e^{i\varphi}\alpha\sin\theta|^2\big)$:

```wl
bs["Probability", "Coherent"[{α, β}] -> "Coherent"[{γ, δ}]]
```

Loss only shrinks the amplitude, $|\alpha\rangle\to|\sqrt{\eta}\,\alpha\rangle$:

```wl
QuantumInterferometer[{"Loss"[η]}]["Coherent"[{α}]]
```

### 12. Hong-Ou-Mandel for any n

A balanced (50:50) beam splitter:

```wl
bs2 = QuantumInterferometer[{"BS"[Pi/4, 0]}];
bs2["Diagram"]
```

$n$ photons enter each port. For every $n$ at once, the amplitude for $k$ out of the first is $A(n,n\to k,2n-k)=D^{\,n}_{0,\;k-n}(0,\pi/2,0)$, which vanishes for every odd $k$: the Hong-Ou-Mandel effect, generalized to any photon number:

```wl
bs2["Amplitude", {n, n} -> {k, 2 n - k}]
```

The same result in traditional notation:

```wl
% // TraditionalForm
```

### 13. Photon number carried through a network

Two beam splitters in series on the same pair of modes:

```wl
chain = QuantumInterferometer[{"BS"[θ1, φ1], "BS"[θ2, φ2]}]
```

```wl
chain["Diagram"]
```

With $n$ photons in, for every $n$, the amplitude is an exact sum over the unseen photon number $K$ between the splitters, $A(n,0\to t,n-t)=\sum_{K=0}^{n}D^{\frac n2}_{\frac n2,\,K-\frac n2}(\varphi_1,2\theta_1,-\varphi_1)\,D^{\frac n2}_{K-\frac n2,\,t-\frac n2}(\varphi_2,2\theta_2,-\varphi_2)$:

```wl
chain["Amplitude", {n, 0} -> {t, n - t}]
```

Substituting $n=2$, $t=1$ and evaluating the sum reproduces the concrete two-photon amplitude exactly; the difference is zero:

```wl
FullSimplify[Activate[chain["Amplitude", {n, 0} -> {t, n - t}] /. {n -> 2, t -> 1}] - chain["Amplitude", {2, 0} -> {1, 1}]]
```

### 14. Three modes, any photon numbers

Three beam splitters across three modes, every angle symbolic:

```wl
bs3 = QuantumInterferometer[{"BS"[a1, b1] -> {1, 2}, "BS"[a2, b2] -> {2, 3}, "BS"[a3, b3] -> {1, 2}}, 3]
```

```wl
bs3["Diagram"]
```

With $N=n_1+n_2$ photons in and $T=t_1+t_2$, the amplitude for arbitrary photon numbers is one conservation delta and a single sum whose upper limit, $\min(N,T)$, follows where photons can flow:

$$A=\delta_{N,\,T+t_3}\sum_{K=0}^{\min(N,T)}D^{\frac N2}_{\frac{n_1-n_2}{2},\;K-\frac N2}(b_1,2a_1,-b_1)\;D^{\frac{N-K}{2}}_{\frac{N-K}{2},\;T-\frac{N+K}{2}}(b_2,2a_2,-b_2)\;D^{\frac T2}_{K-\frac T2,\;\frac{t_1-t_2}{2}}(b_3,2a_3,-b_3)$$

```wl
bs3["Amplitude", {n1, n2, 0} -> {t1, t2, t3}]
```

### 15. Outputs forbidden by symmetry

A four-mode Fourier network, decomposed exactly into beam splitters:

```wl
bsFT = QuantumInterferometer[FourierMatrix[4]]
```

```wl
bsFT["Diagram"]
```

One photon in each mode; every output probability, as an exact rational:

```wl
p = bsFT["Probabilities", {1, 1, 1, 1}]
```

The number of output patterns:

```wl
Length[p]
```

Most never occur, and they come out as exact zeros, not small numbers that floating point would blur into noise. With output modes labelled $0,\dots,3$ and $t_j$ photons in mode $j$, a pattern is forbidden unless $\sum_j j\,t_j\equiv 0 \pmod 4$, the suppression law of Fourier multiports:

```wl
Count[Values[p], 0]
```

### 16. Squeezed light

The chance a squeezer emits one photon pair, exact for complex squeezing $\xi$, $P(0\to2)=\tfrac12\operatorname{sech}|\xi|\,\tanh^2|\xi|$:

```wl
QuantumInterferometer[{"S"[ξ]}]["Probability", {0} -> {2}]
```

Two squeezers meeting at a balanced beam splitter, the building block of Gaussian boson sampling:

```wl
bs4 = QuantumInterferometer[{"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/4, 0]}]
```

```wl
bs4["Diagram"]
```

The coincidence amplitude, $A(0,0\to1,1)=\tfrac12\sqrt{\operatorname{sech}|x_1|\operatorname{sech}|x_2|}\,\big(e^{i\arg x_2}\tanh|x_2|-e^{i\arg x_1}\tanh|x_1|\big)$, vanishes for equal squeezing:

```wl
bs4["Amplitude", {0, 0} -> {1, 1}]
```

### 17. A lossy Mach-Zehnder, exactly

Two beam splitters with loss $\eta$ on the first mode between them, every parameter symbolic:

```wl
l2 = QuantumInterferometer[{"BS"[θ, φ], "Loss"[η] -> 1, "BS"[θ2, φ2]}]
```

```wl
l2["Diagram"]
```

One photon in. Loss enters the exact probabilities: the photon is lost with probability $(1-\eta)\cos^2\theta$, and the fringe term scales with $\sqrt{\eta}$:

$$P(1,0)=\eta\cos^2\theta\cos^2\theta_2+\sin^2\theta\sin^2\theta_2-2\sqrt{\eta}\,\cos\theta\sin\theta\cos\theta_2\sin\theta_2\cos(\varphi-\varphi_2)$$

```wl
l2["Probabilities", {1, 0}]
```

Probability is conserved exactly for any angles and loss, proved symbolically rather than checked at sample points:

```wl
FullSimplify[Total[l2["Probabilities", {1, 0}]], Element[{θ, φ, θ2, φ2}, Reals] && 0 <= η <= 1]
```

### 18. Loss and sampling

Of $n$ photons sent through a lossy channel, $k$ survive with the binomial law $P(n\to k)=\binom{n}{k}\,\eta^{k}(1-\eta)^{n-k}$, derived with $n$, $k$ and $\eta$ all symbolic:

```wl
QuantumInterferometer[{"Loss"[η]}]["Probability", {n} -> {k}]
```

Boson sampling, where each exact probability is a permanent and the full table grows explosively: a random $12\times12$ unitary,

```wl
SeedRandom[7];
mat = RandomVariate[CircularUnitaryMatrixDistribution[12]] // UnitaryMatrix
```

realized as a beam-splitter mesh:

```wl
sample = QuantumInterferometer[mat]
```

```wl
sample["Diagram"]
```

Six photons in, 1000 detector patterns drawn directly with the Clifford-Clifford algorithm, never building the probability table; the twenty most frequent, with their observed frequencies:

```wl
ReverseSort[Normalize[Counts@sample["Samples", {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}, 1000], Total]][[;; 20]]
```
