# Stabilizer Formalism — Test Formula Reference

A consolidated catalogue of formulas, identities, and test states extracted from
the 44 papers under `OngoingProjects/Stabilizer/tex/`. Each entry is a clean,
testable mathematical claim that any stabilizer simulator (CHP-style tableau,
graph-state, Choi-tableau, Wigner, GF(4), check-matrix, qudit, etc.) should
reproduce. References use the directory short-name; e.g. **Got97** is
`Got97_Thesis_quant-ph_9705052/Thesis.tex`.

The reference papers used here are:

**Foundational stabilizer / Clifford simulation:**
- **AarGot04** — Aaronson, Gottesman, *Improved Simulation of Stabilizer Circuits*
- **AndBri05** — Anders, Briegel, *Fast simulation of stabilizer circuits using a graph state representation*
- **Beaudrap11** — de Beaudrap, *A linearized stabilizer formalism for systems of finite dimension*
- **Biswas24** — Biswas, *Step-by-step Gottesman–Knill simulation*
- **DehMoo03** — Dehaene, De Moor, *Clifford group, stabilizer states, and linear and quadratic operations*
- **Gid21** — Gidney, *Stim: a fast stabilizer circuit simulator*
- **Got97 / Got00 / Got09** — Gottesman thesis (1997), pedagogical intros (2000, 2009)
- **Got98** — Gottesman, *The Heisenberg representation of quantum computers*
- **KocHuaLov17** — Kocia, Huang, Love, *Discrete Wigner derivation of the AG tableau algorithm*
- **KoeSmo14** — Koenig, Smolin, *How to efficiently select an arbitrary Clifford group element*
- **Mueller26** — Müller, *PauliEngine: fast arithmetic for Pauli strings*
- **Paler14** — Paler et al., *Tracking quantum errors in stabilizer codes*

**Inner-product / canonical-form / round-trip:**
- **GarMar15 / GarMarCro12** — Garcia, Markov (& Cross), *Stabilizer-state inner products & frames*
- **deSilSalYin23** — de Silva, Salfi, Yin, *Algorithms for stabiliser-state and Clifford-gate conversions*

**Graph states (qubit & qudit):**
- **EllEasCav07** — Elliott, Eastin, Caves, *Graphical description of Clifford operators on stabilizer states* (2007)
- **EllEasCav08** — Elliott, Eastin, Caves, *Graphical description of Pauli measurements on stabilizer states*
- **HeiEisBri03** — Hein, Eisert, Briegel, *Multiparty entanglement in graph states*
- **HeiDurEisRauNesBri06** — Hein, Dür, Eisert, Raussendorf, Van den Nest, Briegel, *Entanglement in Graph States and its Applications*
- **NesDehMoo03** — Van den Nest, Dehaene, De Moor, *Graphical description of local Clifford on graph states*
- **BahBei06** — Bahramgiri, Beigi, *Graph States Under the Action of Local Clifford in Non-Binary Case*

**Pedagogical / textbook QEC introductions:**
- **DevMunNem09** — Devitt, Munro, Nemoto, *Quantum error correction for beginners*
- **Girvin21** — Girvin, *Introduction to QEC and fault tolerance* (Les Houches 2019)
- **Fujii15** — Fujii, *Quantum Computation with Topological Codes* (SpringerBrief)
- **MonPar23** — Monsalve, Parameswaran, *Systematic ECC encoder design from stabilizers*
- **BraDalEva25** — Bradshaw, Dale, Evans, *Introduction to QEC with stabilizer codes* (homological algebra view)
- **BalCenHub20** — Ball, Centelles, Huber, *Quantum error-correcting codes and their geometries*

**Qudits / general Abelian groups:**
- **HosDehMoo04** — Hostens, Dehaene, De Moor — qudit Clifford / stabilizer formalism
- **Gheorghiu11** — Gheorghiu, *Standard form of qudit stabilizer groups* (composite D)
- **BraCheKofKue26** — Brandl, Cherniak, Kofler, Kueng, *QuickQudits — Extended stabilizer tableau for noisy qudit Clifford*
- **BerNes12** — Bermejo-Vega, Van den Nest, *Classical sims of Abelian-group normalizer circuits / Gottesman–Knill for finite Abelian groups*

**Measurement-based QC / fusions / one-way:**
- **BroBri06** — Browne, Briegel, *One-way Quantum Computation tutorial introduction*
- **PatGuh26** — Patil, Guha, *Graphical Rule Book for Clifford Manipulations of Stabilizer States*

**ZX-calculus stabilizer:**
- **Ranchin14** — Ranchin, *Depicting qudit QM and mutually unbiased qudit theories*
- **RanCoe13** — Ranchin, Coecke, *Complete set of circuit equations for stabilizer QM*

**Type theory / symbolic / channel formalism:**
- **FangYing23** — Fang & Ying, *SymPhase: phase symbolization for stabilizer simulation*
- **PayWin24a / WinPay24b** — Payne, Winderl, *λC and λPC: type theories for stabilizer/Clifford*
- **Reid24** — Reid, *Tableau manipulation for noisy circuits*
- **RuhDev25** — Ruh, Devitt, *Pauli tracking library for MBQC scheduling*
- **Winderl23** — Winderl, *Stabilizer-circuit synthesis via inverse symplectic*
- **Yashin25** — Yashin, *Choi-tableau formalism for arbitrary Clifford channels*

Notation conventions in this document: `I, X, Y, Z` are the Pauli matrices;
`H` is the Hadamard; `S = P = diag(1,i)` is the phase gate; `T = diag(1, e^{iπ/4})`
is the π/8 gate; `CNOT_{a→b}` and `CZ_{a,b}` are the standard two-qubit Cliffords.
"Sign-free Pauli string" means a tensor product of `{I,X,Y,Z}` without any phase
factor; "Pauli operator" includes a phase `±1, ±i`. Single-qubit `Y = iXZ`.
Vectors / addition over `Z_2` use `⊕`. The Pauli group on `n` qubits is `P_n`
with `|P_n| = 4^{n+1}`. The Clifford group is `C_n = {U ∈ U(2^n) : UP_nU† = P_n}`.

---

## 0. Quick smoke tests (dimensions, cardinalities, basic checks)

Use these as the simplest correctness gates before any deeper test.

| Quantity | Value | Source |
|---|---|---|
| `\|P_n\|` (full Pauli group, with phase) | `4^{n+1}` | Got97 §3, Got98, AarGot04 |
| `\|P_n / phase\|` (sign-free Pauli strings) | `4^n` | Got98 |
| Number of n-qubit pure stabilizer states | `2^n · ∏_{k=0}^{n-1}(2^{n-k}+1) = 2^{(1/2+o(1))n^2}` | AarGot04 Prop. 2 |
| `\|N(S)\|` for an `[[n,k,d]]` code | `4 · 2^{n+k}` | Got97 §3.2 |
| `\|N(S)/S\|` for an `[[n,k,d]]` code | `4^{k+1}` | Got09 §2.4 |
| Bits to specify a stabilizer state | `n(2n+1)` | Gid21, Got09 |
| Bits to specify the AG tableau (incl. destabilizers) | `2n(2n+1) ≈ 4n²` | AarGot04 |
| Bits per Pauli operator (xz encoding) | `2n + 1` | Gid21 §2.1 |
| Conjugate-tuple bits for a Clifford | `4n² + 2n` | deSilSalYin23 |

**Pauli matrix identities** (Got97, AarGot04, Got98, MonPar23, Biswas24):

```
X² = Y² = Z² = I
XY = iZ      YZ = iX     ZX = iY
YX = -iZ     ZY = -iX    XZ = -iY
Y = iXZ
{X,Y} = {Y,Z} = {Z,X} = 0
[X,Y] = 2iZ  [Y,Z] = 2iX  [Z,X] = 2iY
```

**Single-qubit eigenstate table** (AarGot04 §2):

```
+X : |0⟩+|1⟩       -X : |0⟩-|1⟩
+Y : |0⟩+i|1⟩      -Y : |0⟩-i|1⟩
+Z : |0⟩            -Z : |1⟩
```

`I` stabilizes every state; `−I` stabilizes none.

---

## 1. Pauli encoding round-trips

The core abstraction every package shares: an `(x,z)` bit pair encodes a Pauli
modulo phase.

### 1.1 xz encoding (Got97 §3.4, AarGot04, Gid21)

```
encode_xz(I) = (0,0)   encode_xz(X) = (1,0)
encode_xz(Y) = (1,1)   encode_xz(Z) = (0,1)
```

**Multiplication is XOR mod phase**:
```
encode_xz(P₁ P₂) = (x₁⊕x₂, z₁⊕z₂)        (Gid21 eq. 1)
```

**Round-trip (TEST)**: for every Pauli `P ∈ {I,X,Y,Z}`,
`decode_xz(encode_xz(P)) === P`.

### 1.2 GF(4) encoding (Got97 §3.4, Got00 §5)

`{0, 1, ω, ω²}` corresponds to `{I, Z, X, Y}` with `ω³=1`, `ω+ω²=1`.

```
trace ω = trace ω² = 1     trace 0 = trace 1 = 0
[P,Q] = 0  ⟺  Tr(u · v̄) = 0  (where v̄ swaps ω↔ω²)
```

Linear (vs. additive) GF(4) codes correspond exactly to stabilizer groups
closed under multiplication by `ω` (Got97, Got00).

### 1.3 Symplectic / commutation predicate

For Paulis encoded as `(x | z)` bit-vectors over `Z_2^{2n}` (Got97, AarGot04,
KocHuaLov17, RuhDev25):

```
[P, Q] = 0   ⟺   x_P · z_Q + x_Q · z_P = 0   (mod 2)         (Got97 eq. 25)
```

Equivalently, with the symplectic matrix `J = [[0, I],[I, 0]]`, two rows
`(u | c), (u' | c')` commute iff `[u, u'] := u^T J u' = 0` in `Z_2`
(Yashin25 §2.1, Beaudrap11).

**TEST**: brute-force matrix-level commutation must agree with the symplectic
inner product on every pair of `4^n × 4^n` Paulis.

### 1.4 Phase update under Pauli product (Yashin25, Beaudrap11, AarGot04)

For Pauli observables `P(u | c) = (-1)^c · i^{-⟨z,x⟩} · Z₁^{z₁}X₁^{x₁}…`,

```
P(u|c) · P(u'|c') = i^{β(u,u')} · P(u⊕u' | c⊕c')

β(u,u') = ⟨z⊕z', x⊕x'⟩ - ⟨z,x⟩ - ⟨z',x'⟩ - 2⟨z',x⟩    (Yashin25 §2.1)
```

When `[u,u']=0`, `β` is divisible by 2 and the sign update is `c⊕c'⊕β/2`.

---

## 2. Clifford action on Paulis (transversal "Heisenberg" rules)

These are the absolute foundation. Every simulator must reproduce them
character-for-character.

### 2.1 One-qubit gates (Got97 §5.3, Got98 Table 1, Got09 §2.6, Gid21)

```
H X H† = Z          H Z H† = X         H Y H† = -Y
S X S† = Y          S Z S† = Z         S Y S† = -X
T X T† = (X-Y)/√2   T Z T† = Z         (T ∉ Clifford — see §11)

X X X† = X    X Y X† = -Y    X Z X† = -Z
Y X Y† = -X   Y Y Y† = Y     Y Z Y† = -Z
Z X Z† = -X   Z Y Z† = -Y    Z Z Z† = Z
```

Higher-order single-qubit Cliffords (six total automorphisms of `P_1`):

```
R X R† = Z, R Z R† = X      (R ≡ H)
P X P† = Y, P Z P† = Z       (P ≡ S)
Q X Q† = X, Q Z Q† = -Y      (Q = P† R P, swaps Y↔Z)
T X T† = Y, T Z T† = X       (T = R P†, X→Y→Z cyclic)   (Got97 §5.3)
```

### 2.2 CNOT (control a, target b) — Got98 Table 1

```
X⊗I → X⊗X        I⊗X → I⊗X
Z⊗I → Z⊗I        I⊗Z → Z⊗Z
```

### 2.3 CZ — Got09, AndBri05

```
X⊗I → X⊗Z        I⊗X → Z⊗X
Z⊗I → Z⊗I        I⊗Z → I⊗Z
Y⊗I → Y⊗Z        I⊗Y → Z⊗Y
X⊗X → Y⊗Y · phase(-1)        (i.e. CZ X⊗X CZ = Y⊗Y up to sign)
```

CZ commutes with `Z` on either qubit; both qubits behave symmetrically.

### 2.4 CNOT/CZ identity (Got98, Beaudrap11)

```
CZ_{a,b} = H_b · CNOT_{a→b} · H_b              (basis for graph-state algos)
```

### 2.5 Three CNOTs swap (Got98 fig. 3)

`CNOT_{1→2} · CNOT_{2→1} · CNOT_{1→2}` = SWAP. **TEST**: applying this
sequence to two qubits maps `X⊗I ↔ I⊗X` and `Z⊗I ↔ I⊗Z`.

### 2.6 The "anything in N(P)" 4-qubit transversal gate (Got97 §5.4, eq. 16-17)

A particularly diagnostic 4-qubit Clifford that is a valid transversal
operation for *any* stabilizer code:

```
X⊗I⊗I⊗I → X⊗X⊗X⊗I
I⊗X⊗I⊗I → I⊗X⊗X⊗X
I⊗I⊗X⊗I → X⊗I⊗X⊗X
I⊗I⊗I⊗X → X⊗X⊗I⊗X
Z⊗I⊗I⊗I → Z⊗Z⊗Z⊗I
…   (analogous for Z)
```

If this gate is implemented correctly by a package, conjugating any
stabilizer through it lands inside `S × S × S × S`.

### 2.7 Tableau update primitives (AarGot04 §3, Gid21 §3, KocHuaLov17 §4)

For a stored row `(x_i, z_i, r_i)` and a single gate, the tableau updates
are exactly:

```
Hadamard(a):           r_i ← r_i ⊕ x_{ia} z_{ia};                    swap(x_{ia}, z_{ia})
Phase S(a):            r_i ← r_i ⊕ x_{ia} z_{ia};                    z_{ia} ← z_{ia} ⊕ x_{ia}
CNOT(a→b):             r_i ← r_i ⊕ x_{ia} z_{ib} (x_{ib} ⊕ z_{ia} ⊕ 1)
                       x_{ib} ← x_{ib} ⊕ x_{ia}
                       z_{ia} ← z_{ia} ⊕ z_{ib}
```

**TEST**: construct a random tableau, apply gate `U` then `U†`; tableau must
return to the starting one bit-for-bit.

### 2.8 Rowsum / Pauli string multiplication (AarGot04, KocHuaLov17, Mueller26)

For two rows that anti/commute and store xz bits, the phase function

```
g(0,0,*,*) = 0
g(1,1,x₂,z₂) = z₂ - x₂
g(1,0,x₂,z₂) = z₂(2x₂ - 1)
g(0,1,x₂,z₂) = x₂(1 - 2z₂)
```

drives the rule

```
2 r_h + 2 r_i + Σⱼ g(x_{ij}, z_{ij}, x_{hj}, z_{hj}) ≡ 0 or 2 (mod 4)
                                  → r_h := 0 or 1 respectively
```

with `x_{hj} ← x_{ij} ⊕ x_{hj}`, `z_{hj} ← z_{ij} ⊕ z_{hj}`.

---

## 3. Standard test states & their stabilizers (golden states)

### 3.1 Computational basis — Got98, AarGot04

```
|0…0⟩  : ⟨ +Z₁, +Z₂, …, +Z_n ⟩
|x⟩    : ⟨ (-1)^{x_i} Z_i ⟩_{i=1..n}
```

### 3.2 Plus-state / equal superposition — Got98

```
|+⟩^⊗n : ⟨ +X₁, +X₂, …, +X_n ⟩
|−⟩^⊗n : ⟨ -X₁, -X₂, …, -X_n ⟩
```

### 3.3 Bell states — Got98 §5, Biswas24, PatGuh26

```
|Φ⁺⟩ = (|00⟩ + |11⟩)/√2 : ⟨ +XX, +ZZ ⟩       (also: -YY)
|Φ⁻⟩ = (|00⟩ - |11⟩)/√2 : ⟨ -XX, +ZZ ⟩       (also: +YY)
|Ψ⁺⟩ = (|01⟩ + |10⟩)/√2 : ⟨ +XX, -ZZ ⟩       (singlet has 3 minus signs)
|Ψ⁻⟩ = (|01⟩ - |10⟩)/√2 : ⟨ -XX, -ZZ ⟩       (singlet)
```

**Round-trip TEST**: build via `H_1 · CNOT_{1→2} |00⟩` and verify the
tableau, then convert state ↔ stabilizer in both directions.

### 3.4 GHZ / cat states — Got97 glossary, AarGot04, PatGuh26

```
|GHZ_n⟩ = (|0…0⟩ + |1…1⟩)/√2
       : ⟨ X⊗X⊗…⊗X, Z⊗Z⊗I⊗…, Z⊗I⊗Z⊗I⊗…, …, Z⊗I⊗…⊗I⊗Z ⟩
       (one global X, n-1 pairwise ZZ between qubit 1 and each of 2..n)
```

### 3.5 Cluster / graph states — AndBri05, EllEasCav08, PatGuh26

For graph `G = (V, E)` with adjacency matrix `Γ`,

```
|G⟩ = ∏_{{a,b}∈E} CZ_{a,b} · |+⟩^⊗n            (AndBri05 eq. 8)

K_G^{(a)} = X_a · ∏_{b ∈ N(a)} Z_b              stabilizer generator at vertex a
K_G^{(a)} |G⟩ = |G⟩
```

**TEST**: every graph state has stabilizer rank `n` and the generators
`{K_G^{(a)}}` pairwise commute exactly when `a ≠ b`.

### 3.6 Stabilizer of standard ECC codes (every textbook, MonPar23)

5-qubit `[[5,1,3]]` perfect code (Got97 Tab. 5, Got00 Tab. 5):
```
M₁ = X Z Z X I       Logical:
M₂ = I X Z Z X       X̄ = X X X X X
M₃ = X I X Z Z       Z̄ = Z Z Z Z Z
M₄ = Z X I X Z
```

Steane `[[7,1,3]]` CSS code (Got97 Tab. 7, Got00 Tab. 6):
```
M₁ = X X X X I I I       (rows of Hamming H = parity check)
M₂ = X X I I X X I       Z generators: same support but with Z literals
M₃ = X I X I X I X       (M₄ = Z Z Z Z I I I, M₅, M₆ analogously)
X̄ = I I I I X X X
Z̄ = I I I I Z Z Z
```

Shor `[[9,1,3]]` code (Got97 Tab. 1, Got00 Tab. 4): six `ZZ`s pairing
qubits within each block of three, plus two weight-6 `X` operators:
```
Z Z I I I I I I I            X X X X X X I I I
I Z Z I I I I I I            I I I X X X X X X
I I I Z Z I I I I
I I I I Z Z I I I
I I I I I I Z Z I
I I I I I I I Z Z
```

`[[8,3,3]]` code (Got97 Tab. 4) and `[[4,2,2]]` code (Got97 Tab. 6) are
also standard test fodder; the latter is just `M_1, M_3 M_4` of the
five-qubit code with the last qubit dropped.

**TEST**: package's "preset codes" should reproduce these exact tableaux
(modulo row reordering / row multiplications).

### 3.7 Five-qubit code basis codeword (Got97 §3.2)

```
|0̄⟩ = (Σ_{M ∈ S} M) |00000⟩

Expanded form (16 terms with signs ±1):
  +|00000⟩ +|10010⟩ +|01001⟩ +|10100⟩ +|01010⟩ -|11011⟩ -|00110⟩ -|11000⟩
  -|11101⟩ -|00011⟩ -|11110⟩ -|01111⟩ -|10001⟩ -|01100⟩ -|10111⟩ +|00101⟩
```

**TEST**: reconstructing the dense state vector from the stabilizer must
match this 16-term expansion.

---

## 4. Stabilizer state ↔ density matrix

### 4.1 Pure-state density matrix from stabilizer (AarGot04 §6.1)

For a pure stabilizer state with `n` independent generators `{M_i}`:

```
ρ = (1/2^n) ∏_{i=1}^n (I + M_i) = (1/2^n) Σ_{M ∈ S} M
```

**Mixed stabilizer state** (`r < n` generators):

```
ρ_mix = (1/2^r) ∏_{i=1}^r (I + M_i)        (rank 2^{n−r})
```

**TEST**: trace = 1, hermitian, positive, idempotent (pure case).

### 4.2 Projector onto +1-eigenspace of any Pauli `M ∈ P_n`

```
Π_+ = (I + M) / 2,    Π_- = (I - M) / 2
M = Π_+ - Π_-,    Π_+ + Π_- = I
```

**TEST**: `Π_±² = Π_±`, `Π_+ · Π_- = 0`. Used everywhere in measurement
analysis (Got97, PatGuh26, EllEasCav08).

---

## 5. Stabilizer measurement rules (the canonical update)

The single most-tested behavior in any stabilizer simulator.

### 5.1 Single-Pauli measurement of `M ∈ P_n` on stabilizer state `|ψ⟩` with generators `S = ⟨g_1, …, g_n⟩` (Got97, Got98, Biswas24, PatGuh26)

**Case 1 (deterministic).** If `M` commutes with every `g_i`, then
`±M ∈ S` and the outcome is fixed: `+1` if `M ∈ S`, `-1` if `−M ∈ S`. The
state and tableau are unchanged.

**Case 2 (random).** If `M` anti-commutes with `g_{i_1}, …, g_{i_k}`:

1. Replace `g_{i_2}, …, g_{i_k}` by `g_{i_1} · g_{i_j}` (commutes with `M`).
2. Replace `g_{i_1}` by `(-1)^a · M`, where `a ∈ {0,1}` chosen uniformly.
3. The post-measurement state is the +1 eigenstate of the new generator set.

**TEST**: marginal probability of `a=0` is exactly 1/2; the post-measurement
state is the projection `(I±M) |ψ⟩ / √(⟨ψ|(I±M)/2|ψ⟩)`.

### 5.2 AarGot04 tableau-based measurement (§3, Cases I/II)

Measurement of qubit `a` in the `Z` basis examines column `a` of the
*stabilizer* half of the tableau (rows `n+1 … 2n`):

- **Case I**: ∃ `p ∈ {n+1,…,2n}` with `x_{pa} = 1` ⇒ outcome random.
  Apply `rowsum(i, p)` to all `i ≠ p` with `x_{ia}=1`; copy row `p` into
  row `p−n`; zero row `p` except `r_p ∈ {0,1}` random and `z_{pa}=1`.
- **Case II**: no such `p` ⇒ outcome deterministic. Use scratch row
  `2n+1`, `rowsum(2n+1, i+n)` for every `i ∈ {1..n}` with `x_{ia}=1`;
  return `r_{2n+1}`.

**TEST**: for every initial state and every measurement, marginalizing
over many runs of Case I must recover the exact density-matrix probability.

### 5.3 Graph-state Pauli measurement (AndBri05 §4)

For the graph-state representation `|G; C⟩`, measuring qubit `a`:

- Apply local Clifford that turns `M_a` into `Z`.
- `Z` measurement of vertex `a`: sample outcome; remove all edges
  incident to `a`; right-multiply VOPs of `a` and its neighbors by
  `(X_a ∏_{b∈N(a)} Z_b)^{ζ̃} H_a`.
- `Y` measurement: complement edges among `N(a)`; multiply VOPs of
  `N(a) ∪ {a}` by `√(±iZ)`.
- `X` measurement (with chosen neighbor `b ∈ N(a)`): see AndBri05 eqs.
  for triple symmetric-difference of edge-sets; VOPs receive `√iY` etc.

**TEST**: the post-measurement graph state, converted back to a
stabilizer tableau, agrees with §5.1 / §5.2 outputs.

### 5.4 EllEasCav08 graphical "Z-type" measurement rule

After three simplifications (turn all measured `M_j` into `Z`, reduce
graph, disconnect hollow measured nodes), define

```
M = set of measured nodes
M_S = M \ H            (measured solid)
M_H = M ∩ H            (measured hollow)
M_SE = { j ∈ M_S : |N(j) ∩ M_H| ≡ 0 (mod 2) }
b = |M_H ∩ Z|
```

Then **outcome is deterministic ⟺ M_SE = ∅**, and equals `(-1)^b`. If
`M_SE ≠ ∅`, pick a chosen node `p ∈ M_SE`, apply the four-step graph
rewrite (toggle edges between neighbors of `p` and unchosen `M_SE`,
sign-flip rules, disconnect `p`, make hollow). Outcome is `(-1)^a` with
`a` uniform.

### 5.5 Multi-Pauli (joint) measurement (PatGuh26 §2.5)

For `l` simultaneous Pauli measurements `M = ⟨M_1,…,M_l⟩` on `S = ⟨g_1,…,g_n⟩`:

- If every `M_i` commutes with every `g_k`: state unchanged.
- Else, partition generators so each `M_i` has exactly one anti-commuting
  partner `g_{i,1}` (using product-replacements). Replace `g_{i,1}` with
  `m_i M_i`, leave the rest untouched.

This generalizes single-Pauli measurement and is the engine for fusions
(see §13).

---

## 6. Stabilizer state inner products — **Aaronson–Gottesman theorem**

### 6.1 Orthogonality criterion (GarMarCro12 Th. 1)

Two stabilizer states `|ψ⟩, |φ⟩` are orthogonal iff there exist `P ∈ S(|ψ⟩)`
and `Q ∈ S(|φ⟩)` with `P = -Q`.

### 6.2 Inner-product magnitude (AarGot04 §3 last ¶, GarMarCro12 Th. 2)

If `|ψ⟩` and `|φ⟩` are non-orthogonal stabilizer states, choose generator
sets `{P_i}, {Q_i}` minimizing `s = #{i : P_i ≠ Q_i}`. Then

```
|⟨ψ|φ⟩| = 2^{-s/2}
```

**TEST**: `⟨XX,ZZ | ZI,IZ⟩ = 1/√2` because `⟨ZI,IZ⟩ = ⟨ZI,ZZ⟩` differs
from `⟨XX,ZZ⟩` in only one generator.

**TEST (canonical-form)**: `|⟨0…0 | ψ⟩|² = 2^{-k}` where `k` is the number
of `X/Y` literals in the row-reduced echelon form of `S(|ψ⟩)`.

### 6.3 Brute-force verification

For small `n` (say `n ≤ 5`), every package should round-trip through dense
amplitudes: enumerate all `~2^{n²/2}` stabilizer states, compute the
amplitude inner products, and confirm the formula above.

---

## 7. Stabilizer projection onto codeword (the "sum trick")

Got97 §3.2, Got00, Biswas24, EllEasCav08:

```
|ψ⟩ = (Σ_{M ∈ S} M) |φ⟩
```

is in the codespace for any `|φ⟩`, since multiplying by `M ∈ S` permutes
the sum. **TEST**: pick any starting product state and verify the result is
stabilized by every generator of `S`.

For CSS codes (Got00 §5):

```
|ū⟩ = Σ_{w ∈ C₂^⊥} |u + w⟩,    u ∈ C₁
```

A Hadamard transform on each qubit swaps `C₁ ↔ C₂` and gives the dual form.

---

## 8. CSS codes — the special structure

Got97 §3, Got00 §5, MonPar23, AndBri05.

**Definition**. Stabilizer splits into a pure-`X` block and a pure-`Z`
block. Equivalently, exists `(P, Q)` classical parity-check matrices with
`P · Q^T = 0` over `Z_2`. Distance `d ≥ min(d_P, d_Q)`.

**Transversal CNOT** (Got97 §5.1): bitwise CNOT between two CSS-code
blocks implements the logical CNOT — *iff* the code is CSS. The converse
also holds: bitwise CNOT being a valid encoded CNOT forces CSS structure.

**TEST**: build a 7-qubit Steane block, perform `CNOT_{i,i+7}` for
`i=1..7`, and verify the logical action is exactly logical-CNOT on the
two encoded qubits.

---

## 9. The MacWilliams identity (Got97 §6.2, eq. 24)

For a stabilizer of an `[[n,k,d]]` code, the weight enumerators
`A(z) = Σ A_d z^d` (counts elements of `S` of weight `d`) and
`B(z) = Σ B_d z^d` (elements of `N(S)`) satisfy:

```
B(z) = (1/2^{n-k}) (1+3z)^n A((1-z)/(1+3z))
```

**Companion shadow enumerator** (Got97 §6.2, eq. 25):

```
S(z) = (1/2^{n-k}) (1+3z)^n A((z-1)/(1+3z))
```

Constraints: `A_0 = B_0 = 1`, `B_d ≥ A_d ≥ 0`, `S_d ≥ 0`.

**TEST for `[[5,1,3]]`**: linear programming gives the unique solution
`A = (1, 0, 0, 0, 15, 0)`, `B = (1, 0, 0, 30, 15, 18)` — every five-qubit
distance-3 code must reproduce these counts.

---

## 10. Bounds on quantum codes (Got97, Got00 §4)

| Bound | Statement |
|---|---|
| Quantum Hamming | `(Σ_{j=0}^t 3^j C(n,j)) · 2^k ≤ 2^n` (nondegenerate) |
| Quantum Gilbert–Varshamov | `(Σ_{j=0}^{d-1} 3^j C(n,j)) · 2^k ≤ 2^n` ⇒ `[[n,k,d]]` exists |
| Quantum Singleton (Knill–Laflamme) | `n - k ≥ 2(d-1)` (any code) |
| Asymptotic capacity | `1 - 2p log₂ 3 - H(2p) ≤ R ≤ 1 - p log₂ 3 - H(p)` |

Where `R = k/n`, `p = d/(2n)`, `H(x) = -x log₂ x - (1-x) log₂(1-x)`.

---

## 11. Magic states and the Clifford hierarchy (Got09, AarGot04 §6.3)

**Clifford hierarchy** `C_k`:

```
C_1 = P_n  (Pauli)
C_2 = Cl_n  (Clifford)
C_k = { U : U Q U† ∈ C_{k-1}  ∀Q ∈ C_1 }
```

T-gate is in `C_3`:

```
T X T† = (X - Y)/√2 = e^{iπ/4} X S†
T Z T† = Z
```

**Test (golden magic states)**:

```
|T⟩ ≡ T |+⟩ = (|0⟩ + e^{iπ/4} |1⟩) / √2
|H⟩ ≡ (|0⟩ + cos(π/8)|+⟩ + sin(π/8)|−⟩) / √2     (per Bravyi–Kitaev)
```

For diagonal gate teleportation (Got09 fig.):

```
T |α⟩ = correction · ⟨gates⟩ · ( |+⟩ ⊗ T|+⟩ )
```

If the package supports T-gate "stabilizer-extension" simulation, **TEST**:
each T-gate doubles the size of the stabilizer-decomposition state count.

---

## 12. Symplectic / Clifford ↔ Pauli quotient (Got09 §2.6, Beaudrap11, KocHuaLov17, RuhDev25)

```
Cl_n / (P_n · phases) ≅ Sp(2n, Z_2)
```

A Clifford operation is determined (up to phase and a Pauli) by its
`2n × 2n` symplectic matrix `S` over `Z_2`, satisfying

```
S J S^T = J,    J = [[0, I_n], [I_n, 0]]
```

with inverse `S^{-1} = J S^T J`.

**TEST**: `H, S, CNOT, CZ` produce the standard symplectic generators:

```
H_i  → swap rows i, n+i and negate (over Z_2 trivial)
S_i  → add row i to row n+i
CNOT(a→b) → row b += row a;  row n+a += row n+b
```

KoeSmo14 §2 gives an `O(n³)` constructive algorithm for any symplectic
matrix as a product of those four; **TEST**: synthesize a random
symplectic, decode back via Heisenberg propagation, recover identity.

---

## 13. Bell-state measurement, gate teleportation, fusions

### 13.1 BSM as a Clifford-rotated Z⊗Z measurement (PatGuh26 §2.4)

Circuit `CNOT_{1→2} → H_1 → measure Z_1, Z_2` is exactly the four-outcome
projective measurement onto `{|Φ⁺⟩,|Φ⁻⟩,|Ψ⁺⟩,|Ψ⁻⟩}` — equivalently the
joint stabilizer measurement of `⟨X⊗X, Z⊗Z⟩`.

### 13.2 Quantum teleportation (Got98 §6, Got97, AarGot04 §4)

Initial state: `|ψ⟩_A ⊗ (|00⟩ + |11⟩)_{BC} / √2`.

```
Stabilizers   start: ⟨ I⊗X⊗X, I⊗Z⊗Z ⟩         logicals: X̄ = X⊗I⊗I, Z̄ = Z⊗I⊗I
After CNOT(A→B):     ⟨ I⊗X⊗X, Z⊗Z⊗Z ⟩
After H(A):          ⟨ X⊗X⊗X, Z⊗I⊗I (transformed) ⟩
After Z-measure A:                            ← random outcome m_A
After Z-measure B:                            ← random outcome m_B
End:           qubit C carries  X^{m_B} Z^{m_A} |ψ⟩    (Pauli correction)
```

**TEST**: 1000 sampled runs reconstruct `|ψ⟩` on Bob's qubit after
`X^{m_B} Z^{m_A}` correction, with marginal `(m_A, m_B) ∼ uniform`.

### 13.3 Two-qubit fusion (PatGuh26 §2.6)

For two Bell pairs `|Ψ⁺⟩_{1,2} |Ψ⁺⟩_{3,4}` with stabilizers
`⟨X₁X₂, Z₁Z₂, X₃X₄, Z₃Z₄⟩`, fusion on qubits 2 and 3 with
identities `R_c=R_t=I` and outcomes `m_1, m_3`:

```
Result:  ⟨ m_1 X_2 X_4, m_3 Z_2 Z_4 ⟩  ⇒  Bell pair on (2,4) after correction
```

### 13.4 n-qubit GHZ fusion (PatGuh26 §2.7)

Generalizes to `n` Bell pairs and a star pattern of CNOTs:

```
Stabilizers of n-GHZ on outputs:
  ⟨ X⊗…⊗X,  Z₁Z₂, Z₁Z₃, …, Z₁Z_n ⟩
```

### 13.5 Remote XOR (Got98 §7)

One Bell pair shared between Alice and Bob lets Alice apply CNOT to Bob's
qubit using two CNOTs locally and one bit of classical communication each
way. Final tableau: `Z̄_A → Z̄_A · Z̄_B`, `X̄_B → X̄_A · X̄_B`.

---

## 14. Local-complementation rules for graph states (AndBri05 §3)

Local complementation `L_a` on graph `G`:

```
L_a(V, E) = (V, E ⊕ { {b,c} : b,c ∈ N(a) })
|L_a G⟩ = √(-iX_a) · ∏_{b∈N(a)} √(iZ_b) · |G⟩  ∝ √K_G^{(a)} · |G⟩
```

Together with single-qubit Clifford Vertex Operators (VOPs), local
complementation generates all stabilizer states. There are exactly
**24 single-qubit Cliffords** (the local Clifford group `C_1` mod global
phase) — every package should hard-code this 24-element multiplication
table.

---

## 15. Canonical forms

### 15.1 AarGot04 11-block canonical form (§4 Theorem 8)

Every Clifford circuit `U` has a representative

```
U  ≡  H · C · P · C · P · C · H · P · C · P · C
```

(Hadamard, CNOT, Phase blocks). Corollary: `O(n²/log n)` gate count
suffices via Patel–Markov–Hayes CNOT minimization.

### 15.2 GarMarCro12 5-block canonical form

```
U  ≡  H · C · CZ · P · H        (basis-normalization)
```

Every stabilizer state has a unique mapping to a basis state via this
canonical form.

### 15.3 Dehaene–De Moor quadratic-form triple (DehMoo03, deSilSalYin23)

Every stabilizer state can be written as

```
|s⟩ ∝ Σ_{z ∈ V} (-1)^{Q(z)} i^{ℓ(z)} |z + z₀⟩
```

where `V ⊆ Z_2^n` is a `k`-dim subspace, `z₀ ∈ Z_2^n`, `Q : V → Z_2`
quadratic, `ℓ : V → Z_2` linear. **TEST**: every stabilizer simulator
should be able to convert `(Q, ℓ, z₀, V)` ↔ amplitude vector ↔ check matrix
in poly-time, all three round-tripping exactly.

---

## 16. Choi-tableau for Clifford channels (Yashin25)

Every Clifford channel `Φ : T(H_A) → T(H_B)` has a stabilizer Choi state,
hence is determined by a `k × (2|A| + 2|B| + 1)` Boolean tableau
`[U_A | U_B | c]`. Standard channel-level identities:

```
preparation |0⟩:     [⋅ | 1 0 | 0]                preparation |+⟩: [⋅ | 0 1 | 0]
preparation |1⟩:     [⋅ | 1 0 | 1]                preparation |−⟩: [⋅ | 0 1 | 1]
discard / Tr:        [⋅ ]
identity I:          [[1 0 | 1 0 | 0], [0 1 | 0 1 | 0]]
Z-dephasing:         [1 0 | 1 0 | 0]
X-dephasing:         [0 1 | 0 1 | 0]
Pauli Z:             [[1 0 | 1 0 | 0], [0 1 | 0 1 | 1]]
Pauli X:             [[1 0 | 1 0 | 1], [0 1 | 0 1 | 0]]
Hadamard:            [[1 0 | 0 1 | 0], [0 1 | 1 0 | 0]]
Phase S:             [[1 0 | 1 0 | 0], [0 1 | 1 1 | 0]]
CNOT (2-qubit):      4 rows; column structure mirrors the Heisenberg table
```

**TEST**: composing the tableaus of `H · S · H · S · H · S` should yield
the identity tableau (since `(HS)^3 = ωI`).

---

## 17. Qudit (d-dimensional) extensions

### 17.1 Got97 §3.7 (`C_d D_ω` basis, prime `d`)

```
D_ω |i⟩ = ω^i |i⟩,           C_d |i⟩ = |i+1 mod d⟩
ω = e^{2πi/d} primitive
C_d D_ω = ω D_ω C_d
(C_d^a D_ω^b)(C_d^c D_ω^d) = ω^{ad - bc} (C_d^c D_ω^d)(C_d^a D_ω^b)
```

### 17.2 HosDehMoo04 — generalized Pauli & Clifford for arbitrary `d`

Let `ζ² = ω`, `XZ(a) := X^{v₁}Z^{w₁} ⊗ … ⊗ X^{v_n}Z^{w_n}`. Then

```
ζ^δ · XZ(a) · ζ^ε · XZ(b) = ζ^{δ+ε+2 a^T U b} · XZ(a+b)         (eq. 4)
XZ(a) · XZ(b) = ω^{a^T P b} · XZ(b) · XZ(a),   P = U - U^T

C^T P C = P (mod d)            symplectic constraint
(d-1) Vdiag(C^T U C) + h ≡ 0 (mod 2)        (phase-vector constraint)
```

A Clifford `(C, h)` is determined uniquely up to global phase.

### 17.3 Beaudrap11 — Weyl-operator linearization (any `d`)

```
W_{a,b} = τ^{-ab} Z^a X^b,    τ² = e^{2πi/d},  D = order(τ) ∈ {d, 2d}
W_v · W_w = τ^{[v,w]} W_{v+w}
[v,w] = v^T σ_{2n} w  (symplectic form in Z_d^{2n})
```

The Weyl basis removes quadratic phase dependencies that plague the bare
`Z^a X^b` formalism.

**TEST (qudit smoke)**: for `d = 3, 5, 7`, a single-qudit Fourier `F`
implements `Z ↦ X†`, `X ↦ Z`. Two-qudit SUM gate generalizes CNOT:
`|x⟩|y⟩ → |x⟩|x+y⟩`.

---

## 18. Pauli tracking (RuhDev25, Gid21)

Tracking a single Pauli frame through a Clifford circuit:

```
F → C^{-1} · F · C       (each gate g_i ∈ Cl_n / P_n acts as conjugation)
F → P · F                (each Pauli gate acts by left-multiplication)
```

Memory cost: `O(n)` per frame; time per gate `O(1)` for one-/two-qubit
gates, `O(n²)` worst-case for global Cliffords.

**TEST**: feed a noise-free reference circuit and a noisy circuit through
the simulator with the same RNG seed; bit-by-bit XOR of the two output
streams equals the Pauli-frame contents at each measurement.

---

## 19. Universal identities every test suite should include

A high-coverage minimal set of bool-level checks:

1. `H² = I`, `S² = Z`, `S⁴ = I`, `T⁸ = I`, `(HT)^∞ ≠ identity for any
   finite power` — a useful marker that a putative T-gate is *outside*
   the Clifford group.
2. `(HS)^3 = ωI`,  `H S H = S† X S†`,  `H = S · (S+iX)/√2` — exhaustive
   equivalences among `{H, S, X, Y, Z}`.
3. `(CNOT_{a→b})² = I`,  `(CZ)² = I`.
4. `H_b · CNOT_{a→b} · H_b = CZ_{a,b}`.
5. `CNOT_{a→b} · CNOT_{b→a} · CNOT_{a→b} = SWAP_{a,b}`.
6. `H_a H_b · CZ_{a,b} · H_a H_b = CNOT_{a→b}`.
7. For any `n`-qubit stabilizer state, `Σ_{P ∈ S} P` projects to the
   one-dimensional codespace, `(Σ P)² = 2^n · Σ P`, and tracelessness
   of any `M ∈ P_n \ {±I}` gives `Tr(ρ M) = ±1` if `±M ∈ S` else 0.
8. For Bell state `|Φ⁺⟩`: `⟨Φ⁺ | XX | Φ⁺⟩ = +1`,  `⟨Φ⁺ | YY | Φ⁺⟩ = -1`,
   `⟨Φ⁺ | ZZ | Φ⁺⟩ = +1`. Mismatched eigenvalues are immediate red flag
   for a sign-tracking bug.
9. **Toffoli ∉ Clifford**: `CCX |+⟩|+⟩|0⟩` produces a non-stabilizer
   state (the "magic-augmented" hyper-cube state). Any package claiming
   to handle CCX as a stabilizer gate is wrong.

---

## 20. Round-trip tests across representations

Following deSilSalYin23, every implementation should round-trip:

| Rep. (S₁) Amplitudes ↔ (S₂) Quadratic-form triple ↔ (S₃) Check matrix |
|---|
| (S₁ → S₂): Algorithm 1, time `O(N)` |
| (S₂ → S₁): Algorithm 2 via Gray code, time `O(2^k · n)` |
| (S₃ → S₂): Algorithm 3, time `O(n³)` |
| (S₂ → S₃): Algorithm 4, time `O(n³)` |
| (C₁) Clifford matrix ↔ (C₂) conjugate tuple, time `O(N · n)` |

**TEST**: for `n ≤ 6`, a randomly-sampled Clifford `U` round-tripped
through (C₁ → C₂ → C₁) reproduces `U` exactly (up to global phase).

---

## 21. Stim-style bench identities (Gid21)

Specific assertions Stim's test suite uses:

- `tableau(C_Y) X₁ = X₁ Y₂`, `tableau(C_Y) Z_1 = Z_1 Z_2`,
  `tableau(C_Y) X₂ = Z_1 X_2`, `tableau(C_Y) Z₂ = Z_1 Z_2`.
- For random tableau `T`, both `T · T^{-1}` and `T^{-1} · T` equal
  identity tableau (sign bits all zero).
- Pauli-frame mode ⊕ noise-free reference ≡ noisy direct sample
  (statistically over many shots, exactly bit-by-bit per shot when
  RNG is seeded identically).

---

## 22. Cross-package canonical fixtures

Every test suite should include the following small, hand-verifiable
fixtures (states + tableaux, both `n` small enough for full enumeration):

1. `|0⟩, |1⟩, |+⟩, |−⟩, |+i⟩, |−i⟩` — six 1-qubit stabilizer states.
2. The 4 Bell states (with explicit `(±XX, ±ZZ)` tableaus).
3. 3-qubit GHZ, 3-qubit cluster (path graph), W-state (NOT a stabilizer
   state — should be flagged).
4. Five-qubit code, encoded `|0̄⟩` and `|1̄⟩`.
5. Steane code, encoded `|0̄⟩` and a logical Hadamard test:
   `H^⊗7 |0̄⟩ = |+̄⟩`.
6. Four-qubit error-detecting code (Got98 fig. 6) — single-qubit error
   detection only, with explicit syndrome lookup table.

---

## 23. Caveats / pitfalls (warnings to encode as failing tests)

Drawn from across Got97, AarGot04, Beaudrap11, deSilSalYin23, Yashin25:

- The "phase-bit" of a Pauli operator is *not* the same as its global
  phase; mixing them up breaks deterministic-measurement signs (esp. in
  AarGot04 case II).
- Two equivalent generator sets can have different `(x | z | r)` matrices
  — equivalence must be checked group-theoretically (does same group?),
  not by row-equality.
- The `Y = iXZ` convention differs in some older papers (Got97 used a
  variant with `Y = -iXZ`). Conversions to Heisenberg form must declare
  which convention is active.
- For `d` even, qudit Weyl operators `W_v` and `W_{v+d·e}` differ by a
  sign (Beaudrap11 eq. 18). A simulator that simply "reduces mod d" will
  silently invert phases.
- Three Bell-state stabilizer choices `⟨XX, ZZ⟩`, `⟨XX, -YY⟩`, `⟨ZZ, -YY⟩`
  describe the *same* state; orthogonality tests must be group-aware.
- The W state `(|001⟩ + |010⟩ + |100⟩)/√3` is NOT a stabilizer state. A
  package that "successfully" produces a tableau for the W state has a
  silent bug.
- `-I^⊗n` is excluded from any valid stabilizer group; if it ever appears
  during reduction, the stabilizer is degenerate (zero state).

---

## 24. Useful concrete numerical fixtures

For final amplitude-level testing, package these dense states (each is a
golden expected output):

```
|+⟩       = [1, 1] / √2
|GHZ_3⟩   = [1, 0, 0, 0, 0, 0, 0, 1] / √2
|cluster_3⟩  (linear chain) = [+, +, +, +, +, +, +, +] / 2√2     before signs
                 (signs follow Γ_{ij} adjacency; explicit table in AndBri05)
|0̄⟩_{[5,1,3]} = (1/4) · 16-term sum from §3.7 above
|0̄⟩_{[7,1,3]} = (1/√8) · [|0000000⟩ + |1111000⟩ + |1100110⟩ + |1010101⟩
                          + |0011110⟩ + |0101101⟩ + |0110011⟩ + |1001011⟩]
```

---

## 25. End-to-end "stabilizer Olympics" checklist

A package passing all of the following can be considered formula-conformant:

1. ✅ Encodes / decodes `(x, z, r)` for every Pauli matching §1.
2. ✅ Symplectic commutation predicate matches matrix commutator on all
   `4^n × 4^n` pairs for `n ≤ 4`.
3. ✅ Reproduces every entry of the Heisenberg propagation table (§2)
   for `H, S, S†, X, Y, Z, T-conjugation-by-Clifford, CNOT, CZ`.
4. ✅ All 8 stabilizers of the standard 5-qubit and 7-qubit codes are
   produced exactly (sign included).
5. ✅ Bell, GHZ, cluster, and `[[5,1,3]]` codeword states round-trip
   between tableau ↔ amplitudes ↔ quadratic-form triple.
6. ✅ MacWilliams identity holds for the explicit `(A_d, B_d)` of every
   built-in code.
7. ✅ Inner product `|⟨ψ|φ⟩|² = 2^{-s}` holds for 100 random pairs at
   `n ≤ 6`.
8. ✅ Single-qubit and joint Pauli measurements obey §5 and produce the
   correct deterministic vs. random classification, and the random
   marginal is ½.
9. ✅ Quantum teleportation circuit reproduces the input on Bob with a
   `X^{m_B} Z^{m_A}` correction; statistics over 1000 shots are uniform.
10. ✅ Toffoli on `|+⟩|+⟩|0⟩` flagged as "non-stabilizer" (or routed to
    a magic-state extension).
11. ✅ Local-complementation rules reproduce AndBri05 graph dynamics.
12. ✅ Choi-tableau composition (Yashin25 Tab. I) of `H S H S H S` =
    identity tableau (i.e. `(HS)^3 = ωI` mod global phase).
13. ✅ Qudit (`d=3, 5`) Weyl-operator commutation `W_a W_b = ω^{[a,b]} W_b W_a`
    holds for all `a, b ∈ Z_d^{2n}`.

A failure on any single item localizes the bug to one of the surveyed
papers' rules, providing a directly testable regression case.

---

## 26. Sign-convention pitfalls (testable failure modes)

Drawn from a 2026 audit of pedagogical sources; each item is a place where naive
implementations silently get a sign wrong.

### 26.1 `Y_L = -Y₁Y₂Y₃` for the 3-qubit repetition code (Girvin21 §QEC)

For `|0_L⟩ = |000⟩`, `|1_L⟩ = |111⟩`:
```
X_L = X₁X₂X₃,    Z_L = Z₁Z₂Z₃,    Y_L = i X_L Z_L = -Y₁Y₂Y₃
```
Three factors of `i` from `Y_i = iX_iZ_i` combine with the global `i` to give `−1`.

**TEST**: `[X_L, Y_L] = 2iZ_L` only when the minus sign is included.

### 26.2 Steane transversal `S^⊗7` enacts logical `S†`, NOT `S` (DevMunNem09 §IV-D)

Bitwise `P^⊗7` (where `P = S = diag(1,i)`) on a Steane-encoded state does NOT
enact logical `S`; it enacts logical `S†`. Codes with `Z_L = Z^⊗n` and
`n ≡ 3 (mod 4)` flip the sign.

**TEST**: Apply `S^⊗7` to `|+_L⟩` of Steane; the resulting stabilizer change
is `X^⊗7 → −i (XZ)^⊗7`, matching `S†` not `S`.

### 26.3 The 5-qubit `[[5,1,3]]` code is non-CSS — `H^⊗5` is NOT transversal (DevMunNem09 §IV)

Under bitwise Hadamard:
```
K¹ = XZZXI →  ZXXZI       K² = IXZZX →  IZXXZ
K³ = XIXZZ →  ZIZXX       K⁴ = ZXIXZ →  XZIZX
```
None of these images are in `⟨XZZXI, IXZZX, XIXZZ, ZXIXZ⟩` (the original group).

**TEST**: `H^⊗5` should NOT preserve the [[5,1,3]] stabilizer; any package that
silently accepts it as a logical Hadamard has a bug.

### 26.4 Reed–Muller `[[15,1,3]]` transversal `T` enacts logical `T†` (Fujii15 §6)

```
T^⊗15 |0_L⟩ = e^{+iπ/8} |0_L⟩,    T^⊗15 |1_L⟩ = e^{-iπ/8} |1_L⟩
```
i.e. the bitwise `T` realizes logical `T†` (not `T`). Same pattern as Steane-`S`.

### 26.5 `Y` convention (Got97 vs. modern)

`Y = iXZ` (modern, Gid21, AarGot04) but Got97 occasionally uses `Y = −iXZ`.
Heisenberg conjugation tables can disagree by an overall `−` on every `Y` row
unless the convention is pinned.

---

## 27. Additional codeword amplitude expansions

Concrete dense-state fixtures, useful as numerical golden expected outputs.

### 27.1 `[[4,2,2]]` error-detecting code (DevMunNem09 §IV-A)

```
|00⟩_L = (|0000⟩ + |1111⟩)/√2
|01⟩_L = (|1100⟩ + |0011⟩)/√2
|10⟩_L = (|1010⟩ + |0101⟩)/√2
|11⟩_L = (|0110⟩ + |1001⟩)/√2
```
All four are +1 eigenstates of `XXXX` and `ZZZZ`. Logicals (Got97 §3.6):
`X̄₁ = X₁X₂`, `X̄₂ = X₁X₃`, `Z̄₁ = Z₁Z₃`, `Z̄₂ = Z₁Z₂`.

**TEST**: `X̄_i |00⟩_L = |0..1..0⟩_L` (with i-th bit flipped), `Z̄_i |b₁b₂⟩_L = (-1)^{b_i}|b₁b₂⟩_L`.

### 27.2 Hexacode `[[6,0,4]]_2` weight enumerator (BalCenHub20 §8.2)

```
A(x,y) = x⁶ + 45 x²y⁴ + 18 y⁶
```
Self-dual: `B = A`. Built from the GF(4) hexacode `[6,3,4]_4` via the F4-trick.

**TEST**: stabilizer group of the hexacode contains exactly `1 + 45 + 18 = 64 = 2⁶`
elements; weight distribution matches `(1, 0, 0, 0, 45, 0, 18)`.

### 27.3 5-qubit code alternative form (DevMunNem09 §IV-B)

A different 16-term `|0⟩_L` for `[[5,1,3]]` (related to the Got97 §3.7 form by
local rotations):
```
|0⟩_L = (1/4) [+|00000⟩+|01010⟩+|10100⟩-|11110⟩
              +|01001⟩-|00011⟩-|11101⟩-|10111⟩
              +|10010⟩-|11000⟩-|00110⟩-|01100⟩
              -|11011⟩-|10001⟩-|01111⟩+|00101⟩]
```
constructed via `|0⟩_L = ∏_{i=1..4}(I + Kⁱ)/√16 |00000⟩` from stabilizers
`XZZXI, IXZZX, XIXZZ, ZXIXZ`.

### 27.4 Smallest planar code `m=1, n=2` — 5 qubits (BraDalEva25 §5.4)

Stabilizers (4 generators):
```
Z₀Z₂Z₃,   Z₁Z₂Z₄,   X₀X₁X₂,   X₂X₃X₄
```
Logicals: `X̄ = X₀X₃`, `Z̄ = Z₀Z₁`. Logical states:
```
|0⟩_L = ½(|10010⟩ + |01110⟩ + |10101⟩ + |01001⟩)
|1⟩_L = ½(|00000⟩ + |11101⟩ + |00111⟩ + |11011⟩)
```
A complete 16-entry syndrome table maps each of the 15 single-qubit Pauli
errors to a unique 4-bit syndrome (with one degeneracy: `X₀` and `X₃` have the
same syndrome `1000`).

### 27.5 Toric code logical-zero (BraDalEva25 §5.3)

On an `m × n` torus:
```
|00⟩_L ∝ Σ_{b ∈ ∂C₂(L*)} X_{d_v + b} |0…0⟩
|10⟩_L ∝ Σ_b X_b |0…0⟩
|01⟩_L ∝ Σ_b X_{d_h + d_v + b} |0…0⟩
|11⟩_L ∝ Σ_b X_{d_h + b} |0…0⟩
```
where `d_v, d_h` are non-trivial cycles on the dual lattice.

---

## 28. Graph-state algebraic identities (NesDehMoo03, EllEasCav07, HeiDurEisRauNesBri06)

### 28.1 LC at vertex `i` — explicit formula on adjacency matrix (NesDehMoo03 Thm. 2)

```
g_i(θ) = θ + θ Λ_i θ + Λ
```
where `Λ_i` is the matrix with 1 at `(i,i)` and zero elsewhere, and `Λ` is
diagonal chosen to zero the diagonal of the result. As a fractional-linear
map on symmetric binary matrices, with `Q_i = [[I, diag(θ_i)], [Λ_i, I]]`:
```
g_i(θ) = (Aθ + B)(Cθ + D)^{-1}
```

### 28.2 Edge-LC formula (EllEasCav07 eq. 8; NesDehMoo03 Lemma 2)

`g_{jk} = g_j g_k g_j` (3-fold product) acts on adjacency by:
```
Γ'_{lm} = Γ_{lm} + (Γ_{jl}+δ_{jl})(Γ_{km}+δ_{km})
                 + (Γ_{jm}+δ_{jm})(Γ_{kl}+δ_{kl})
```
This expression is symmetric in `j, k`.

### 28.3 GHZ orbit count (NesDehMoo03 Final note)

Under local Clifford, the n-qubit GHZ state's orbit contains exactly `n+1`
inequivalent graphs: the `n` stars (each centered at a different vertex) plus
the complete graph `K_n`.

A graph state with adjacency `θ` is LC-equivalent to GHZ iff
`rank θ_i = 1 = rank θ_{i,j}` for all single-vertex / vertex-pair submatrices.

### 28.4 Reduced-state rank ↔ submatrix rank (NesDehMoo03; HeiEisBri03 eq. 947)

For graph state `|G⟩` with adjacency `θ`, partition into `A ∪ B`:
```
log₂ rank tr_A[|G⟩⟨G|] = rank_{F₂} Γ_{AB}
```
where `Γ_{AB}` is the off-diagonal block. Polynomial-time exact Schmidt-rank
via GF(2) row reduction.

### 28.5 Reduced state via stabilizer subgroup (HeiDurEisRauNesBri06 eq. reduced_GS_1)

```
ρ_G^A = tr_B|G⟩⟨G| = (1/2^|A|) Σ_{σ ∈ S_A} σ
```
where `S_A = {σ ∈ S(|G⟩) : supp(σ) ⊆ A}` (stabilizer elements supported in `A`).
Also: `(ρ_G^A)² = (|S_A|/2^|A|) ρ_G^A` (proportional to a projector).

### 28.6 Pauli error → Pauli-Z error on graph state (HeiDurEisRauNesBri06 eq. 3737)

```
σ_x^a |G⟩ = σ_z^{N_a} |G⟩,    σ_y^a |G⟩ = ±i σ_z^{N_a + a} |G⟩
```
Any Pauli channel acting on `|G⟩⟨G|` produces a graph-diagonal mixed state.

### 28.7 Three-vertex GHZ-style identity (HeiDurEisRauNesBri06 eq. 3170)

For `a, b, c` pairwise adjacent in graph `G`:
```
K_a K_b K_c = − σ_x^a σ_x^b σ_x^c σ_z^{N_a + N_b + N_c}
```
The `−1` sign drives Bell-inequality violations for graph states.

### 28.8 Graph-state projector and basis (HeiDurEisRauNesBri06 eq. 734)

```
|G⟩⟨G| = (1/2^N) Σ_{σ ∈ S} σ
|W⟩ := σ_z^W |G⟩,    K_a |W⟩ = (-1)^{W_a} |W⟩
```
The set `{|W⟩ : W ⊆ V}` is an orthonormal basis of `(C²)^V` (graph-state basis).

### 28.9 Pauli measurement rules on a graph state (HeiEisBri03 eq. proposition `pauli`)

- `σ_z` measurement at `a`:  `G' = G − {a}`,   `U_{z,+} = I`, `U_{z,−} = ∏_{b∈N_a} σ_z^b`
- `σ_y` measurement at `a`:  `G' = G Δ E(N_a, N_a) − {a}`,   `U_{y,±} = ∏_{b∈N_a} (∓i σ_z^b)^{1/2}`
- `σ_x` measurement at `a` with chosen neighbor `b₀ ∈ N_a`:
  ```
  G' = G  Δ E(N_{b₀}, N_a)  Δ E(N_{b₀} ∩ N_a, N_{b₀} ∩ N_a)  Δ E({b₀}, N_a − {b₀})  − {a}
  U_{x,+} = (+iσ_y^{b₀})^{1/2} ∏_{b ∈ N_a − N_{b₀} − {b₀}} σ_z^b
  U_{x,−} = (−iσ_y^{b₀})^{1/2} ∏_{b ∈ N_{b₀} − N_a − {a}} σ_z^b
  ```

### 28.10 Stabilizer-state graph decoration scheme (EllEasCav07 §III.B)

A stabilizer state's graph is the simple graph from adjacency `B` (right half of
canonical `(I|B)` form), plus three decorations:
- `H` on qubit `j` → hollow node at vertex `j`
- `S` on qubit `j` → loop on vertex `j`
- `Z` on qubit `j` → minus sign on vertex `j`

Every stabilizer state (not just graph states) has such a representation.
Reduced stabilizer graphs (those with the minimum number of hollow nodes among
all equivalent graphs) are a canonical-form invariant.

### 28.11 Bipartite Schmidt-measure formulas (HeiEisBri03)

For graph `G` on `N` vertices with adjacency `θ`:
- `E_S(|G⟩) ≤ E_S(|G'⟩) + 1` and `E_S(|G'⟩) ≤ E_S(|G⟩)` under any single-Pauli measurement.
- Even-`N` ring: `E_S = N/2`.
- 3D cluster `(L_x, L_y, L_z)`: `E_S = ⌊L_x L_y L_z / 2⌋`.

---

## 29. Topological codes (chain-complex view, BraDalEva25 §5)

### 29.1 Stabilizer construction from a 2-cell complex

Given cell complex `X = (V, E, F)` with chain complex
`0 → C₂ →^{∂₂} C₁ →^{∂₁} C₀ → 0` over `Z₂`:
```
Z_f = ∏_{e ∈ ∂₂(f)} Z_e          (one Z-stab per 2-cell)
X_v = ∏_{e ∈ ∂*(v)} X_e          (one X-stab per 0-cell)
[Z_c, X_d] = 0  ⟺  c · d = Σ c_i d_i = 0 (mod 2)        (intersection number)
```

### 29.2 Logical operators ↔ first homology (BraDalEva25 Prop. 5.1)

```
N(S)/S ≅ H₁(L, Z₂) ⊕ H₁(L*, Z₂)
d = min{ ℓ(c) : 1 ≠ c ∈ H₁(L, Z₂) ∪ H₁(L*, Z₂) }
```

### 29.3 Toric code parameters

- `[[2mn, 2, min(m,n)]]` on `m × n` torus
- `2mn − 2` independent stabilizer generators (one redundancy per cell type)
- 4-fold ground-state degeneracy (logical qubit count `k = 2`)
- Constraints: `∏_s A_s = ∏_p B_p = I`

### 29.4 Planar code parameters

- `[[2mn + n − m, 1, min(m,n)]]` on `m × n` lattice (smooth top/bottom, rough sides)
- `mn` `Z_f` plaquettes + `(m+1)(n−1)` `X_v` stars = `2mn − m + n − 1` generators

### 29.5 Smallest planar code `m=1, n=2` (5 qubits)

(See §27.4 above for the full state-vector data.)

**TEST**: build the 4 stabilizers, verify `dim(codespace) = 2`; build the 16-entry
syndrome lookup table; verify each Pauli error maps to a unique 4-bit syndrome.

---

## 30. Quantum MacWilliams identity — Shor–Laflamme form (BalCenHub20 §8.2)

### 30.1 Shor–Laflamme weight enumerators

```
A_j = Σ_{wt(E) = j} Tr(E P) Tr(E† P)
B_j = Σ_{wt(E) = j} Tr(E P E† P)
```
where `P` is the codespace projector and `wt(E)` is Pauli weight.

### 30.2 The MacWilliams duality

```
q^n B(x, y) = A(x + (q² − 1) y, x − y)
q^n A(x, y) = B(x + (q² − 1) y, x − y)
```
For qubits (`q = 2`): `2^n B(x,y) = A(x + 3y, x − y)`.

### 30.3 Distance certification

A code has distance `d` iff `K · B_j = A_j` for `0 ≤ j < d`, where `K = dim P`.

### 30.4 Pure-code self-dual identity

For `[[n,0,d]]` codes (k = 0): `A(x, y) = B(x, y)` (self-dual weight enumerator).

### 30.5 TEST — hexacode `[[6,0,4]]_2`

Compute `(A_0, …, A_6) = (1, 0, 0, 0, 45, 0, 18)`; verify `B = A`; verify the
polynomial `x^6 + 45 x² y^4 + 18 y^6` is fixed by `2^6 B(x,y) = A(x + 3y, x − y)`.

---

## 31. Extended quantum bounds

### 31.1 Quantum Hamming (BraDalEva25 eq. 3.1)

For nondegenerate `[[n, k, d]]` codes correcting `t = ⌊(d−1)/2⌋` errors:
```
Σ_{j=0}^{t} C(n, j) · 3^j · 2^k ≤ 2^n
```
Verify saturation for `[[5,1,3]]` (perfect code): `2 + 5·3·2 = 32 = 2^5`.

### 31.2 Quantum Singleton — strong form (BalCenHub20 Thm. 7.3)

```
n ≥ k + 2(d − 1)        for any [[n,k,d]] code with k ≥ 1
```
Forbids small-`d` codes that classical Singleton would allow (e.g. `[[3,1,2]]`).

### 31.3 Quantum-MDS dimension bound (BalCenHub20 eq. 7.5)

If a QMDS `[[n, n−2d+2, d]]_q` exists, then
```
n ≤ q² + d − 2
```

### 31.4 QMDS partial-trace family (BalCenHub20 §7.4)

If `[[n, n−2d+2, d]]` QMDS exists, so does `[[n−s, n−2d+2+s, d−s]]` for
`0 ≤ s ≤ d`. For every `S ⊂ {1,…,n}` with `|S| ≤ (n+k)/2`:
```
tr_{S^c}(P) ∝ I
```

### 31.5 Threshold theorem clean form (DevMunNem09 §VII)

For `g`-level concatenation with a `[[n, k, 3]]` code:
```
p_L^g = (c p)^{2^g} / c,    p_th = 1/c
```

---

## 32. Diagonal-unitary canonical decomposition (BroBri06)

### 32.1 Pauli-Z exponential expansion (BroBri06 eq. 15)

Any n-qubit diagonal unitary can be written (up to global phase) as
```
D_n = ∏_{m ∈ {0,1}^n \ {0}} exp[ i (φ_m / 2) Z_1^{m_1} Z_2^{m_2} … Z_n^{m_n} ]
```
i.e. a product over the `2^n − 1` non-trivial Pauli-Z strings.

### 32.2 Toffoli-Z (CCZ) coefficients (BroBri06 Fig. 7)

```
CCZ = D_3   with   θ_{100} = θ_{010} = θ_{001} = θ_{111} = −π/4,
                   θ_{110} = θ_{101} = θ_{011} = +π/4
```
Plug into eq. (15) and exponentiate; result is `diag(1,1,1,1,1,1,1,−1)`.

### 32.3 Three-pattern arbitrary single-qubit unitary (BroBri06 §2.1)

On a 4-qubit linear graph, measurement angles `(φ_1, φ_2, φ_3)` with adaptive
corrections `φ_2 = (−1)^{m_1} β`, `φ_3 = (−1)^{m_2} γ` realize:
```
U_z(γ) U_x(β) U_z(α)        (modulo a Pauli by-product)
```
deterministically, regardless of measurement outcomes `m_1, m_2, m_3`.

### 32.4 ZZ-rotation pattern (3-qubit star)

Central ancilla `a` connected to qubits 1, 2; measuring `a` in basis
`(−φ, −π/2)` realizes `e^{−i (φ/2) Z_1 Z_2}` on the input qubits, with optional
`Z_1 Z_2` byproduct on outcome `1`.

### 32.5 2n-qubit Clifford pattern

Any n-qubit Clifford is implementable on a 2n-qubit one-way pattern using only
Pauli measurements; no adaptivity needed. Strengthens Gottesman–Knill.

---

## 33. Composite-D qudit caveats (Gheorghiu11, BraCheKofKue26)

### 33.1 Composite-D size theorem (Gheorghiu11 Thm. 1)

```
K · |S| = D^n
```
for any qudit stabilizer code with arbitrary integer `D ≥ 2`. `K = dim C` need
not be a power of `D`.

**Counterexample**: for `D = 4`, `S = ⟨X², Z²⟩` has `|S| = 4`, `K = 1` (it
stabilizes a single state, not 4); needs **2 generators**, not 1.

### 33.2 Composite-D needs up to 2n generators

For composite `D`, a stabilizer can require up to `2n` generators (not `n`).
The prime-D claim "always `n` generators" fails because non-prime `Z_D` is not
a field.

### 33.3 Standard form via Smith Normal Form (Gheorghiu11 Thm. 2)

Any qudit stabilizer over `Z_D` reduces (by Clifford + row ops) to:
```
S' = | M    0    | Z_1   Z_3 |    rank r ≤ n
     | 0    0    | Z_2   Z_4 |    last k − r rows
```
with `M = diag(m_1, …, m_r)`, `m_i | D`. For prime `D`: `M = I_r`, `Z_2 = 0`,
`Z_1` symmetric (recovers Got97/Nielsen–Chuang form).

### 33.4 Qudit SWAP factorization (Gheorghiu11 Tab. I)

```
SWAP_{ab} = CNOT_{ab} · CNOT_{ba}^† · CNOT_{ab} · (F_a² ⊗ I_b)
```
The extra `F²` is necessary for `D ≥ 3`.

### 33.5 QuickQudits τ-phase convention (BraCheKofKue26 §2)

```
ω = e^{2πi/d},   τ = exp(πi(d² + 1)/d),   τ² = ω
Y := τ^{-1} X Z
W(a, b) = τ^{-ab} X^a Z^b
W(a, b) · W(a', b') = τ^{a' b − a b'} W(a + a', b + b')
[W(a, b), W(a', b')] = 0  ⟺  a' b − a b' ≡ 0 (mod d)
```

### 33.6 Tableau update branching for `S` (BraCheKofKue26 Tab. 2)

```
S applied to qudit i:    z_i ← z_i + x_i;    r ← r + Δ
                         Δ = x_i²            (d even)
                         Δ = x_i (x_i − 1)   (d odd)
```
Branching by parity of `d` is **required** for correctness.

### 33.7 Eigenspace projector for prime p (BahBei06 Lemma 3)

For `g = ω^c X(a) Z(b) ∈ G_n` non-trivial:
```
P_j = (1/p) [I + ω^{-j} g + ω^{-2j} g² + … + ω^{-(p-1)j} g^{p-1}]
```
All `p` projectors have equal rank `(dim H)/p = p^{n-1}` per qudit.

### 33.8 Power formula for qudit Pauli (BahBei06 Lemma 2)

For odd prime `p` and `g = ω^c X(a) Z(b)`:
```
g^k = ω^{(k c + C(k,2) a · b)} · X(k a) Z(k b)
```
In particular `g^p = I` (using `p | C(p, 2)` for odd `p`).

---

## 34. Stabilizer formalism over arbitrary finite Abelian groups (BerNes12)

For `G = Z_{d_1} × … × Z_{d_m}`:

### 34.1 Pauli operators

```
X(g) |h⟩ = |h + g⟩
Z(g) |h⟩ = χ_g(h) |h⟩,    χ_g(h) = exp(2πi Σ_i g(i) h(i) / d_i)
σ(a, g, h) = γ^a Z(g) X(h),    γ = e^{iπ / |G|},    a ∈ Z_{2|G|}
```

### 34.2 Commutation predicate

```
σ(a₁, g₁, h₁) σ(a₂, g₂, h₂) commute  ⟺  χ_{g₁}(h₂) = χ_{g₂}(h₁)
```

### 34.3 Stabilizer-code dimension theorem (BerNes12 Cor. 6.3)

For any stabilizer group `S` over Abelian `G` with codespace `V`:
```
dim V = |G| / |S|
```
generalizing `K · |S| = D^n`. `V` is one-dimensional iff `H = D^⊥` (where `H`
is the X-label subgroup, `D` the Z-only subgroup).

### 34.4 Stabilizer-state normal form (BerNes12 Thm. 8)

Every stabilizer state has the form
```
|φ⟩ = α (1/√|H|) Σ_{h ∈ H} ξ(h) |s + h⟩
```
where `α` is a global phase, `s ∈ G`, `H` is the X-label subgroup, and
`ξ : H → C*` is a quadratic function:
```
ξ(g + h) = ξ(g) ξ(h) B(g, h)        for some bilinear B
ξ(g)^{2|G|} = 1
ξ(0) = 1
```
For `G = Z₂^m` recovers the Dehaene–De Moor form
`(−1)^{q(x)} i^{ℓ(x)} |x + s⟩`.

### 34.5 Pauli measurement update (BerNes12 Thm. 7)

After measuring `σ` on a stabilizer state with eigenvalue `λ`:
```
S_m = ⟨ λ̄ σ, C_S(σ) ⟩
```
where `C_S(σ)` is the centralizer of `σ` inside `S`.

### 34.6 Adaptive > unitary normalizer circuits (BerNes12 §7)

The state `|ψ⟩ = (|0⟩ + |2⟩)/√2 ∈ C^4` (a coset state of `H = ⟨2⟩ < Z_4`)
**cannot** be obtained by any unitary normalizer circuit applied to `|0⟩`
(its stabilizer needs ≥ 2 generators, while a unitary preserves generator
count). It IS reachable with measurement.

This is a canonical witness that adaptive normalizer circuits are strictly
more powerful than unitary ones.

### 34.7 Subgroup-orthogonality lemma (BerNes12 Lem. 4.1)

For `H, K ≤ G`:
```
(H^⊥)^⊥ = H,    |H^⊥| = |G|/|H|,    H ⊆ K ⟺ K^⊥ ⊆ H^⊥,
(H ∩ K)^⊥ = ⟨H^⊥, K^⊥⟩
```

---

## 35. ZX-calculus stabilizer circuit identities (RanCoe13, Ranchin14)

A complete set of stabilizer-circuit equations follows from ZX-completeness;
the most directly testable are unitary-level identities.

### 35.1 π-commutation

```
X · R_Z(α) = R_Z(−α) · X
Z · R_X(α) = R_X(−α) · Z
```

### 35.2 Hadamard Euler decomposition

```
H = R_Z(π/2) · R_X(π/2) · R_Z(π/2)        (up to global phase)
```

### 35.3 Spider fusion (S1, S2)

Two same-color spiders connected by a wire fuse: their phases add modulo 2π.
For stabilizer-restricted angles, the phase group is `Z_4 = {0, π/2, π, 3π/2}`.

### 35.4 CNOT-elimination on |0⟩-control or |+⟩-target (B1)

```
CNOT (|0⟩_c ⊗ |ψ⟩_t) = |0⟩ ⊗ |ψ⟩
CNOT (|ψ⟩_c ⊗ |+⟩_t) = |ψ⟩ ⊗ |+⟩
```

### 35.5 CNOT classical-bit identity (K1)

```
CNOT (|1⟩_c ⊗ |ψ⟩_t) = |1⟩ ⊗ X|ψ⟩
CNOT (|ψ⟩_c ⊗ |−⟩_t) = (Z|ψ⟩) ⊗ |−⟩
```

### 35.6 Phase composition rules

```
R_Z(α) · R_Z(β) = R_Z(α + β)
R_X(α) · R_X(β) = R_X(α + β)
```

### 35.7 D-torus picture for qudit stabilizer states (Ranchin14 §6)

For prime `D`, every single-qudit stabilizer state is a vector in one of
`D + 1` mutually unbiased bases (axes `Z, X, XZ, XZ², …, XZ^{D−1}`).
The stabilizer states unbiased to `Z` form a `D-torus` parametrized by
`(α₁, …, α_{D−1}) ∈ (Z_D)^{D−1}`, with `D²` points per torus.

For `D = 3` (qutrit): exactly **12 single-qutrit stabilizer states**, falling on
4 MUBs of size 3.

---

## 36. Updated "Stabilizer Olympics" — extended conformance checklist

The following items extend §25's 13-point checklist with formulas added in 2026.
A package passing **all 23 items** is "thoroughly conformant":

14. ✅ `Y_L = −Y₁Y₂Y₃` for the 3-qubit repetition code (sign convention).
15. ✅ `S^⊗7` on Steane enacts `S†` (logical phase-conjugation).
16. ✅ `H^⊗5` does NOT preserve the [[5,1,3]] stabilizer (non-CSS check).
17. ✅ `[[4,2,2]]` codeword amplitudes match §27.1 explicitly.
18. ✅ NesDehMoo03 LC formula `g_i(θ) = θ + θ Λ_i θ + Λ` reproduces the standard
    local-complementation graph operation on every `n ≤ 6` graph state.
19. ✅ Reduced-state rank `log₂ rank tr_A|G⟩⟨G| = rank Γ_{AB}` holds for
    every bipartition of every connected `n ≤ 5` graph state.
20. ✅ Quantum Hamming bound saturates exactly for `[[5,1,3]]`, fails to
    saturate for Steane and Shor (degenerate codes).
21. ✅ Hexacode `[[6,0,4]]_2` weight enumerator matches MacWilliams self-dual
    fixed point `B = A = x^6 + 45 x²y^4 + 18 y^6`.
22. ✅ Toric code on `m × n` torus reports `[[2mn, 2, min(m,n)]]` parameters
    with 4-dim ground-space.
23. ✅ Smallest planar code `m=1, n=2` reproduces the §27.4 amplitudes and
    16-entry syndrome lookup.
