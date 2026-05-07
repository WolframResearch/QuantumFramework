# Stabilizer Formalism — Test Formula Reference

A consolidated catalogue of formulas, identities, and test states extracted from
the 28 papers under `OngoingProjects/Stabilizer/tex/`. Each entry is a clean,
testable mathematical claim that any stabilizer simulator (CHP-style tableau,
graph-state, Choi-tableau, Wigner, GF(4), check-matrix, qudit, etc.) should
reproduce. References use the directory short-name; e.g. **Got97** is
`Got97_Thesis_quant-ph_9705052/Thesis.tex`.

The reference papers used here are:

- **AarGot04** — Aaronson, Gottesman, *Improved Simulation of Stabilizer Circuits*
- **AndBri05** — Anders, Briegel, *Fast simulation of stabilizer circuits using a graph state representation*
- **Beaudrap11** — de Beaudrap, *A linearized stabilizer formalism for systems of finite dimension*
- **Biswas24** — Biswas, *Step-by-step Gottesman–Knill simulation*
- **DehMoo03** — Dehaene, De Moor, *Clifford group, stabilizer states, and linear and quadratic operations*
- **EllEasCav08** — Elliott, Eastin, Caves, *Graphical description of Pauli measurements on stabilizer states*
- **FangYing23** — Fang & Ying, *SymPhase: phase symbolization for stabilizer simulation*
- **GarMar15 / GarMarCro12** — Garcia, Markov (& Cross), *Stabilizer-state inner products & frames*
- **Gid21** — Gidney, *Stim: a fast stabilizer circuit simulator*
- **Got97 / Got00 / Got09** — Gottesman thesis (1997), pedagogical intros (2000, 2009)
- **Got98** — Gottesman, *The Heisenberg representation of quantum computers*
- **HosDehMoo04** — Hostens, Dehaene, De Moor — qudit Clifford / stabilizer formalism
- **KocHuaLov17** — Kocia, Huang, Love, *Discrete Wigner derivation of the AG tableau algorithm*
- **KoeSmo14** — Koenig, Smolin, *How to efficiently select an arbitrary Clifford group element*
- **MonPar23** — Monsalve, Parameswaran, *Systematic ECC encoder design from stabilizers*
- **Mueller26** — Müller, *PauliEngine: fast arithmetic for Pauli strings*
- **PatGuh26** — Patil, Guha, *Stabilizer formalism for Clifford+fusion linear-optical quantum computing*
- **Paler14** — Paler et al., *Tracking quantum errors in stabilizer codes*
- **PayWin24a / WinPay24b** — Payne, Winderl, *λC and λPC: type theories for stabilizer/Clifford*
- **Reid24** — Reid, *Tableau manipulation for noisy circuits*
- **RuhDev25** — Ruh, Devitt, *Pauli tracking library for MBQC scheduling*
- **Winderl23** — Winderl, *Stabilizer-circuit synthesis via inverse symplectic*
- **Yashin25** — Yashin, *Choi-tableau formalism for arbitrary Clifford channels*
- **deSilSalYin23** — de Silva, Salfi, Yin, *Algorithms for stabiliser-state and Clifford-gate conversions*

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
