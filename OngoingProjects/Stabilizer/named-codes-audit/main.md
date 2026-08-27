# Named Stabilizer Codes for the `PauliStabilizer` Catalog

*A naming survey, a proposed catalog, and a contributor checklist.*

Wolfram `QuantumFramework` maintenance notes, August 2026.

> **Abstract.** The shipped `PauliStabilizer` engine of the Wolfram `QuantumFramework` paclet (version 2.1.1) carries a small set of named quantum error-correcting codes: the five-qubit perfect code [[5,1,3]], the Steane code [[7,1,3]], and the Shor code [[9,1,3]], each stored as a codeword state. This document surveys how stabilizer codes are named and referred to in the literature, records canonical stabilizer generators for the common small codes with a traceable citation for each, and proposes a prioritized list of codes worth adding to the catalog. It closes with a step-by-step checklist a contributor can follow to add a named code, matching the conventions already used in `NamedCodes.m`. No generator set is quoted from memory: every one has been reproduced against a primary reference and checked, as a stabilizer state, against the live kernel.

## 1. How stabilizer codes are named in the literature

### 1.1 The [[n,k,d]] convention

A quantum stabilizer code is specified by the triple [[n,k,d]]: $n$ is the number of physical qubits, $k$ the number of logical (encoded) qubits, and $d$ the code distance, the minimum weight of a Pauli operator that maps one codeword to another while commuting with every stabilizer. The double-bracket notation distinguishes a quantum code from a classical [n,k,d] code and was made standard by Calderbank, Rains, Shor and Sloane [CRSS 1998] and by Gottesman's thesis [Gottesman 1997]. A distance $d$ code detects up to $d-1$ errors and corrects up to $\lfloor (d-1)/2 \rfloor$; thus $d=1$ codes detect nothing (they are pure states or classical repetition memories), $d=2$ codes detect a single error but correct none, and $d=3$ is the smallest distance that corrects an arbitrary single-qubit error.

A stabilizer code on $n$ qubits with $k$ logical qubits is defined by $n-k$ independent, mutually commuting Pauli operators (the stabilizer generators, or parity checks) that generate an abelian group $S$ not containing $-I$. The code space is the simultaneous $+1$ eigenspace of $S$, of dimension $2^k$. Adding $k$ further commuting, independent Pauli operators (a choice of logical $\overline{Z}_1,\dots,\overline{Z}_k$) pins a single state within that space; this is the representation the shipped `PauliStabilizer` object uses (Section 1.4).

### 1.2 Folklore names and their aliases

The same code is very often referred to by several names: an author's name, a qubit count, a geometric picture, or a classical-coding pedigree. The most frequent are collected in Table 1. The aliases matter for a catalog because a user is as likely to reach for `"ShorCode"` as for `"9QubitCode"`. Two aggregators are used throughout as cross-checks: the [Error Correction Zoo](https://errorcorrectionzoo.org) and [Grassl's code tables](http://www.codetables.de).

**Table 1. Common stabilizer codes and the names they go by.**

| Parameters | Common names / aliases | Primary references |
|---|---|---|
| [[3,1,1]] | bit-flip code; (majority) repetition code; three-qubit code | [Shor 1995; Nielsen & Chuang 2010] |
| [[3,1,1]] | phase-flip code; dual repetition code (Hadamard of the bit-flip code) | [Shor 1995; Nielsen & Chuang 2010] |
| [[4,2,2]] | $C_4$ code; four-qubit code; "iceberg" (m=2); little Shor code | [Vaidman et al. 1996; Grassl et al. 1997; Rains 1999] |
| [[2m,2m-2,2]] | iceberg code; [[n,n-2,2]] error-detecting code | [Steane 1996b; Gottesman 1997] |
| [[5,1,3]] | five-qubit code; perfect code; Laflamme code; DiVincenzo-Shor code | [Laflamme et al. 1996; Bennett et al. 1996; DiVincenzo & Shor 1996] |
| [[7,1,3]] | Steane code; seven-qubit code; CSS Hamming code (r=3) | [Steane 1996a; Steane 1996b] |
| [[8,3,2]] | smallest interesting color code; 3D color code; hypercube code | [Kubica et al. 2015; Campbell 2016; Vasmer & Kubica 2022] |
| [[9,1,3]] | Shor code; nine-qubit code; concatenated [[3,1,1]] code | [Shor 1995] |
| [[9,1,3]] | rotated surface code (d=3); Surface-17 | [Tomita & Svore 2014; Bombin & Martin-Delgado 2007; Horsman et al. 2012] |
| [[9,1,3]] | Bacon-Shor code (subsystem code, gauge qubits) | [Bacon 2006] |
| [[13,1,3]] | unrotated (standard) surface code (d=3); planar code | [Bravyi & Kitaev 1998; Dennis et al. 2002; Fowler et al. 2012] |
| [[15,1,3]] | quantum Reed-Muller code; tetrahedral color code; punctured RM code | [Steane 1999; Bravyi & Kitaev 2005] |
| [[15,7,3]] | quantum Hamming code (r=4); CSS Hamming code | [Steane 1996b; Calderbank & Shor 1996] |
| [[23,1,7]] | Golay code; quantum Golay code | [Golay 1949; Grassl 2007] |
| family | toric code; surface code (planar) | [Kitaev 2003; Bravyi & Kitaev 1998; Fowler et al. 2012] |
| family | color code (triangular, 6.6.6 / 4.8.8) | [Bombin & Martin-Delgado 2006] |

Naming pitfalls worth stating once:

- "Shor code" and "9-qubit code" are the same object, but "9-qubit code" is badly ambiguous: at least three distinct codes share [[9,1,3]], namely Shor's concatenated code, the $d=3$ rotated surface code, and the Bacon-Shor code. The first two are stabilizer codes; the Bacon-Shor code is a *subsystem* code with gauge qubits [Bacon 2006], so it is not a single stabilizer state and falls outside the `PauliStabilizer` state form.
- "Quantum Hamming code" is used for two unrelated families: the CSS codes [[2^r-1, 2^r-1-2r, 3]] built from the classical Hamming code (Steane's family, of which [[7,1,3]] is $r=3$ and [[15,7,3]] is $r=4$) [Steane 1996b], and the non-CSS family [[2^r, 2^r-r-2, 3]] that saturates the quantum Hamming bound (Gottesman's family, of which $r=4$ is [[16,11,3]]) [Gottesman 1996]. These must be kept apart.
- "Surface code" without qualification usually means the rotated variant in current hardware papers, but the unrotated [[13,1,3]] is the original planar layout. The rotated and unrotated distance-3 codes differ in $n$ (9 vs 13).
- The five-qubit code is unique up to local Clifford equivalence, so its many names all denote one code; the presentation as a cyclic code (Table 2) is one representative.

### 1.3 Codes versus states: what `PauliStabilizer` stores

A code with $k \geq 1$ is a subspace, not a single state. The shipped `PauliStabilizer` object, however, is a *state*: it stores $n$ generators, so it represents one codeword rather than the whole code space. The shipped named codes finesse this by fixing the logical value: the $n-k$ parity checks are followed by $k$ logical $\overline{Z}$ operators, which selects the logical all-zeros codeword $|0_L\rangle$ (or $|0\cdots0_L\rangle$). For $k=1$ this is natural and the $|1_L\rangle$ partner is obtained by flipping the sign of the single logical $\overline{Z}$. For $k>1$ the state form must pin every logical qubit, and a single sign-flipped variant reaches only one of the $2^k$ logical basis states; this is a real limitation of the state representation, flagged again in Section 2 for each $k>1$ entry. (The separate, not-yet-shipped `QEC` package stores a code as its parity checks only, i.e. $n-k$ generators, sidestepping the choice of logical state.)

### 1.4 The `PauliStabilizer` generator convention

The conventions below were read from the `Constructors.m` and `NamedCodes.m` sources under `Kernel/Stabilizer/` and confirmed against the live 2.1.1 kernel.

- **Signed Pauli strings.** A generator is a string over the alphabet `{I,X,Y,Z}` with an optional leading sign `"+"` or `"-"`; the constructor maps `"-"` to phase $-1$ and `"+"` (or no prefix) to $+1$. Example: `"-ZZZZZ"`.
- **Qubit 1 is leftmost.** Character position $i$ (counting from the left, after any sign) acts on qubit $i$. A single-qubit $Z$ on qubit 1 of a three-qubit register is `"ZII"`.
- **Order: parity checks, then logical $\overline{Z}$.** The first $n-k$ generators are the code's parity checks; the last $k$ are the logical $\overline{Z}_1,\dots,\overline{Z}_k$. This ordering is what makes the base name the $|0_L\rangle$ state.
- **Base state and its sign-flip partner.** The base name builds $|0_L\rangle$. The `"...1"` variant flips the sign of the last generator (the logical $\overline{Z}$), giving $|1_L\rangle$. For the five-qubit code, `"5QubitCode"` ends in `"ZZZZZ"` and `"5QubitCode1"` ends in `"-ZZZZZ"`.
- **Logical $\overline{Z}$ need not be $Z$-type.** For the standard CSS-oriented codes the last generator is a $Z$-string, but the label "logical $\overline{Z}$" means "the operator whose $\pm 1$ eigenvalue labels the logical basis", which for the phase-flip code is an $X$-string (Section 2.2).

The three shipped named codes are reproduced in Table 2 for reference. A useful internal consistency check: the Steane $X$-checks `IIIXXXX`, `XIXIXIX`, `IXXIIXX` are exactly the three parity-check hyperplanes of the classical [7,4] Hamming code (qubit $i$ sits at the binary expansion of $i$, and check $j$ collects the qubits whose bit $j$ is set). The same indexing generates the proposed [[15,7,3]] and [[15,1,3]] entries in Section 2, so those follow the shipped convention exactly.

**Table 2. Currently shipped named codes (`NamedCodes.m`), as $|0_L\rangle$ states.** The `"...1"` variant of each flips the sign of the last (logical $\overline{Z}$) generator.

| Name | Parameters | Generators (qubit 1 leftmost) |
|---|---|---|
| `5QubitCode` | [[5,1,3]] | `XZZXI` `IXZZX` `XIXZZ` `ZXIXZ` `ZZZZZ` |
| `SteaneCode` (= `7QubitCode`) | [[7,1,3]] | `IIIXXXX` `XIXIXIX` `IXXIIXX`<br>`IIIZZZZ` `ZIZIZIZ` `IZZIIZZ` `ZZZZZZZ` |
| `9QubitCode` (Shor) | [[9,1,3]] | `ZZIIIIIII` `IZZIIIIII` `IIIZZIIII`<br>`IIIIZZIII` `IIIIIIZZI` `IIIIIIIZZ`<br>`XXXXXXIII` `IIIXXXXXX` `ZZZZZZZZZ` |

## 2. Proposed catalog

Codes are prioritized by three criteria: (a) how often they are referenced by name, (b) whether they are small enough to be useful in an exact small-$n$ setting (roughly $n \leq 20$ to 30; the stabilizer tableau itself is only $O(n^2)$, but the named states are meant to be read, written, and turned into dense vectors, which is $2^n$), and (c) whether a canonical generator set exists. Every generator set below was checked against the live kernel: it builds a valid `PauliStabilizer` state on the stated number of qubits (`StabilizerStateQ` returns `True`). Code distances are taken from the cited literature, not recomputed, because distance is a property of the code space, not of the single stored state.

### 2.1 Priority 0: aliases for codes already shipped

The cheapest and highest-value change adds *names*, not codes. The ecosystem's own `QEC` package already calls the Shor code `"ShorCode"`, and the five-qubit code is at least as often called the perfect or Laflamme code. Proposed aliases (pointing at the existing rules):

| New alias | Resolves to |
|---|---|
| `"ShorCode"` | `"9QubitCode"` |
| `"PerfectCode"`, `"LaflammeCode"` | `"5QubitCode"` |
| `"HammingCode"` (with a doc note on the ambiguity) | `"7QubitCode"` |

### 2.2 Priority 1: the repetition codes (d=1)

The three-qubit bit-flip and phase-flip codes are the standard first example in every textbook [Shor 1995; Nielsen & Chuang 2010] and are currently unnamed; users must type the strings by hand. Both are $k=1$, so they fit the shipped $|0_L\rangle$/$|1_L\rangle$ convention exactly.

**Table 3. Repetition codes, as $|0_L\rangle$ states.** For the phase-flip code the logical $\overline{Z}$ is the $X$-string `XXX`: flipping its sign gives $|1_L\rangle = |---\rangle$.

| Proposed name | Parameters | Parity checks | Logical $\overline{Z}$ |
|---|---|---|---|
| `"BitFlipCode"` (`"RepetitionCode"`) | [[3,1,1]] | `ZZI` `IZZ` | `ZZZ` |
| `"PhaseFlipCode"` | [[3,1,1]] | `XXI` `IXX` | `XXX` |

The full $|0_L\rangle$ generator lists are therefore `{ZZI, IZZ, ZZZ}` and `{XXI, IXX, XXX}`; the $|1_L\rangle$ variants replace the last entry with `"-ZZZ"` and `"-XXX"`.

### 2.3 Priority 2: distance-3 codes commonly named (k=1)

These correct one error, are $k=1$ (so they fit the convention cleanly), and are named constantly in current work.

#### Rotated surface code, Surface-17, [[9,1,3]]

Distinct from the Shor code of the same parameters. Data qubits sit on a $3\times 3$ rotated lattice; the stabilizers are weight-2 on the boundary and weight-4 in the bulk. The explicit generators below are the Tomita-Svore Surface-17 layout [Tomita & Svore 2014] as tabulated by the [Error Correction Zoo](https://errorcorrectionzoo.org); the geometric construction is due to Bombin and Martin-Delgado [Bombin & Martin-Delgado 2007]. *Flag: the exact strings depend on the lattice-labelling convention (which boundaries carry $X$ vs $Z$ checks); the set below is one valid, verified representative, not the only presentation.* A logical $\overline{Z} = Z_3 Z_6 Z_9$ (a $Z$-string joining the two $Z$-boundaries) completes the $|0_L\rangle$ state.

**Table 4. Surface-17 [[9,1,3]] $|0_L\rangle$, data qubits 1 to 9.** The last generator is the logical $\overline{Z}$.

| Type | Support | String |
|---|---|---|
| $X$ check | $X_2 X_8$ | `IXIIIIIXI` |
| $X$ check | $X_4 X_5$ | `IIIXXIIII` |
| $X$ check | $X_3 X_6 X_7 X_8$ | `IIXIIXXXI` |
| $X$ check | $X_1 X_4 X_6 X_9$ | `XIIXIXIIX` |
| $Z$ check | $Z_1 Z_9$ | `ZIIIIIIIZ` |
| $Z$ check | $Z_3 Z_7$ | `IIZIIIZII` |
| $Z$ check | $Z_4 Z_5 Z_6 Z_7$ | `IIIZZZZII` |
| $Z$ check | $Z_2 Z_6 Z_8 Z_9$ | `IZIIIZIZZ` |
| $\overline{Z}$ | $Z_3 Z_6 Z_9$ | `IIZIIZIIZ` |

#### Quantum Reed-Muller code, [[15,1,3]]

The tetrahedral color code, famous for a transversal $T$ gate [Steane 1999; Bravyi & Kitaev 2005]. Index qubit $i$ ($1 \le i \le 15$) by the 4-bit binary expansion of $i$. The four $X$-checks are the coordinate hyperplanes (weight 8); the ten $Z$-checks are the same four hyperplanes (weight 8) plus six weight-4 faces. A logical $\overline{Z}$ is any weight-3 $Z$-string along a tetrahedron edge; $Z_1 Z_2 Z_3$ works and completes $|0_L\rangle$. Supports are listed in Table 5; the same hyperplane indexing generates the shipped Steane checks, so the convention is identical.

**Table 5. Quantum Reed-Muller [[15,1,3]] $|0_L\rangle$: 14 parity checks and one logical $\overline{Z}$, by qubit support.**

| Generator | Support (qubit indices) |
|---|---|
| $Z$ check (bit 1) | 1, 3, 5, 7, 9, 11, 13, 15 |
| $Z$ check (bit 2) | 2, 3, 6, 7, 10, 11, 14, 15 |
| $Z$ check (bit 3) | 4, 5, 6, 7, 12, 13, 14, 15 |
| $Z$ check (bit 4) | 8, 9, 10, 11, 12, 13, 14, 15 |
| $Z$ face | 3, 7, 11, 15 |
| $Z$ face | 5, 7, 13, 15 |
| $Z$ face | 6, 7, 14, 15 |
| $Z$ face | 10, 11, 14, 15 |
| $Z$ face | 12, 13, 14, 15 |
| $Z$ face | 9, 11, 13, 15 |
| $X$ check (bit 1) | 1, 3, 5, 7, 9, 11, 13, 15 |
| $X$ check (bit 2) | 2, 3, 6, 7, 10, 11, 14, 15 |
| $X$ check (bit 3) | 4, 5, 6, 7, 12, 13, 14, 15 |
| $X$ check (bit 4) | 8, 9, 10, 11, 12, 13, 14, 15 |
| $\overline{Z}$ (edge) | 1, 2, 3 |

### 2.4 Priority 3: detection codes and higher-rate codes (k>1)

These are common and small, but $k>1$: as `PauliStabilizer` *states* they must pin all logical qubits to $|0\rangle$, and the `"...1"` sign-flip reaches only one of the $2^k$ logical basis states. *Flag for every entry in this section: the logical-operator basis is a convention, so the $|0\cdots0_L\rangle$ state is one verified representative; a doc note should say which logical qubit the sign-flip toggles.* An alternative is to expose these as their parity checks only (as the `QEC` package does) and leave the logical completion to the user.

#### C4 / iceberg code, [[4,2,2]]

The smallest error-*detecting* code [Vaidman et al. 1996; Grassl et al. 1997; Rains 1999]. Parity checks `XXXX`, `ZZZZ`. A verified $|00_L\rangle$ completion adds logical $\overline{Z}_1 = Z_3 Z_4$ (`IIZZ`) and $\overline{Z}_2 = Z_2 Z_4$ (`IZIZ`), giving the state $(|0000\rangle + |1111\rangle)/\sqrt{2}$.

#### Iceberg family, [[2m,2m-2,2]]

Parity checks are the single all-$X$ and all-$Z$ strings on $2m$ qubits [Steane 1996b; Gottesman 1997]. A verified $|0\cdots0_L\rangle$ completion uses $\overline{Z}_i = Z_1 Z_{i+1}$ for $i = 1,\dots,2m-2$. For [[6,4,2]] the full list is `{XXXXXX, ZZZZZZ, ZZIIII, ZIZIII, ZIIZII, ZIIIZI}`; for [[8,6,2]] extend the all-$X$/all-$Z$ pair with $Z_1 Z_{i+1}$ up to $i=6$.

#### Smallest interesting color code, [[8,3,2]]

Eight qubits on the vertices of a cube; admits a transversal (non-Clifford) $CCZ$ [Kubica et al. 2015; Campbell 2016; Vasmer & Kubica 2022]. The five parity checks ([Error Correction Zoo](https://errorcorrectionzoo.org) tableau) are `ZZIIIIZZ`, `ZIZZIIZI`, `IIZIZIZZ`, `ZIZIIZIZ`, `XXXXXXXX`. A verified $|000_L\rangle$ completion adds three cube-edge logical operators $\overline{Z}_1 = Z_1 Z_2$ (`ZZIIIIII`), $\overline{Z}_2 = Z_1 Z_3$ (`ZIZIIIII`), $\overline{Z}_3 = Z_1 Z_5$ (`ZIIIZIII`).

#### Quantum Hamming code, [[15,7,3]]

The $r=4$ member of Steane's CSS Hamming family (the shipped Steane code is $r=3$) [Steane 1996b; Calderbank & Shor 1996]. The $X$- and $Z$-checks are both the four parity-check hyperplanes of the classical [15,11] Hamming code, i.e. the four "bit $j$" supports listed at the top of Table 5 (as $X$ then as $Z$). It has $k=7$, so a state completion needs seven logical $\overline{Z}$; this makes it a poor fit for the state form and it is better offered as a code (parity checks only) or as a family generator `HammingCSSCode[r]` mirroring the `QEC` package.

### 2.5 Optional / large: reference entries

#### Unrotated surface code, [[13,1,3]]

The original planar layout [Bravyi & Kitaev 1998; Dennis et al. 2002; Fowler et al. 2012]: 13 data qubits, weight-3 and weight-4 star/plaquette checks. Worth documenting as the "other" distance-3 surface code, but a specific string set is layout-dependent and should be drawn from a single cited figure before being frozen into a rule.

#### Golay code, [[23,1,7]]

$k=1$ and distance 7, the smallest single-logical code that corrects three errors [Golay 1949; Grassl 2007]. It is a CSS code built from the self-dual-containing cyclic [23,12,7] binary Golay code. *Flag: the 22 explicit Pauli strings depend on the chosen generator matrix of the classical code*, so they are not tabulated here; a contributor should lift a fixed parity-check matrix from Grassl's tables [Grassl 2007] and verify the resulting state before adding a rule. At $n=23$ it is near the upper edge of the "small" range and is best marked as a stress-test / reference entry.

#### Families (out of scope for fixed strings)

The toric and color codes [Kitaev 2003; Bombin & Martin-Delgado 2006] and the hypergraph-product / quantum-LDPC families are parametrized by a lattice or a pair of classical codes; they belong behind a generator function (`SurfaceCode[d]`, `ToricCode[L]`, `ColorCode[d]`) rather than as fixed named strings, and are noted here only so the naming survey is complete.

### 2.6 Summary of the proposal

**Table 6. Proposed additions, in priority order.** "State fit" notes how well the code matches the shipped $|0_L\rangle$/$|1_L\rangle$ state convention.

| Priority | Proposed name(s) | Parameters | State fit |
|---|---|---|---|
| 0 | `ShorCode`, `PerfectCode`, `LaflammeCode` | (aliases) | exact |
| 1 | `BitFlipCode` / `RepetitionCode` | [[3,1,1]] | exact ($k=1$) |
| 1 | `PhaseFlipCode` | [[3,1,1]] | exact ($k=1$, $\overline{Z}$ is $X$-type) |
| 2 | `RotatedSurfaceCode` / `Surface17` | [[9,1,3]] | exact ($k=1$) |
| 2 | `ReedMullerCode` / `TetrahedralCode` | [[15,1,3]] | exact ($k=1$) |
| 3 | `C4Code` / `IcebergCode` | [[4,2,2]] | pin \|00_L⟩ ($k=2$) |
| 3 | iceberg family | [[2m,2m-2,2]] | pin \|0…0_L⟩ |
| 3 | `ColorCode832` | [[8,3,2]] | pin \|000_L⟩ ($k=3$) |
| 3 | `HammingCode154` | [[15,7,3]] | poor ($k=7$); prefer code/family |
| opt | `SurfaceCode13` | [[13,1,3]] | exact ($k=1$), layout-dependent |
| opt | `GolayCode` | [[23,1,7]] | exact ($k=1$), basis-dependent |

## 3. How to add a named code: a contributor checklist

Adding a named code touches exactly two places in the kernel plus one test file. The description below is prose only; it edits no code in this document.

1. **Fix the generator set from a citation.** Choose the $n-k$ parity checks and $k$ logical $\overline{Z}$ operators from a primary reference. Write them as signed Pauli strings, qubit 1 leftmost, in the order "parity checks, then logical $\overline{Z}_1,\dots,\overline{Z}_k$" (Section 1.4). Keep the reference; it will go in the code comment.

2. **Verify the set is a valid stabilizer state before writing the rule.** In a scratch kernel (not in the paclet sources), load the working tree and confirm the strings build a state:

   ```wl
   PacletDirectoryLoad["/path/to/QuantumFramework"];
   Needs["Wolfram`QuantumFramework`"];
   ps = PauliStabilizer[{"...", "...", ...}];
   ps["Qubits"]           (* should equal n *)
   StabilizerStateQ[ps]   (* True iff generators commute and are independent *)
   ```

   If `StabilizerStateQ` is `False`, the generators either fail to commute or are not independent; fix the set before proceeding.

3. **Register the name(s) in `$PauliStabilizerNames`.** Add the new name string, and any aliases and the `"...1"` variant, to the `$PauliStabilizerNames` list at the top of `NamedCodes.m`. This is the catalog of recognized names; keep related names grouped on one line as the existing entries are.

4. **Add the resolution rule(s).** Below the catalog, add a rule of the form

   ```wl
   PauliStabilizer["MyCode"] :=
       PauliStabilizer[{"g1", "g2", ..., "Zbar"}]
   ```

   Group aliases with an alternative pattern, matching the shipped style, for example

   ```wl
   PauliStabilizer["MyCode" | "MyAlias"] :=
       PauliStabilizer[{ ... }]
   ```

   For a $k=1$ code, add the $|1_L\rangle$ partner by flipping the sign of the last generator:

   ```wl
   PauliStabilizer["MyCode1"] :=
       PauliStabilizer[{"g1", ..., "-Zbar"}]
   ```

5. **Comment with the physics, not the change.** Precede the rule with a short comment naming the code, its [[n,k,d]], and the reference section (as the existing `Got97`, `Got00`, Nielsen-Chuang comments do). Describe what the rule builds, not that it was added.

6. **For $k>1$, document the state form.** State explicitly, in the comment and the doc page, that the rule builds $|0\cdots0_L\rangle$, that the logical-operator basis is a convention, and which logical qubit the `"...1"` sign-flip toggles. Consider whether the code is better exposed as parity checks only.

7. **Add tests in `Tests/Stabilizer/`.** Put correctness tests next to the existing ones. `PauliStabilizer.wlt` holds the constructor and property tests; `Correctness_TextbookResults.wlt` holds against-the-literature checks. At minimum assert:

   ```wl
   PauliStabilizer["MyCode"]["Qubits"] === n
   StabilizerStateQ[PauliStabilizer["MyCode"]] === True
   PauliStabilizer["MyCode"]["Stabilizers"]     (* matches intended strings *)
   ```

   and, for $k=1$, that the $|0_L\rangle$ and $|1_L\rangle$ variants are distinct and orthogonal. A shape-only test (right number of qubits) is not enough: it passes on a wrong-matrix-right-shape set, so assert the actual generator strings.

8. **Update the doc page and the name catalog audit.** If the symbol's documentation lists the named codes, add the new name and its aliases there, and record the count change so the audit's "`$PauliStabilizerNames` unchanged" note stays accurate.

## 4. References

Every stabilizer generator set quoted in this document is traceable to one of the entries below. arXiv identifiers are given where the paper predates or accompanies a journal version and the identifier has been verified.

- **[Bacon 2006]** D. Bacon, "Operator quantum error-correcting subsystems for self-correcting quantum memories," Physical Review A 73, 012340 (2006). [arXiv:quant-ph/0506023](https://arxiv.org/abs/quant-ph/0506023).
- **[Bennett et al. 1996]** C. H. Bennett, D. P. DiVincenzo, J. A. Smolin, and W. K. Wootters, "Mixed-state entanglement and quantum error correction," Physical Review A 54, 3824 (1996). [arXiv:quant-ph/9604024](https://arxiv.org/abs/quant-ph/9604024).
- **[Bombin & Martin-Delgado 2006]** H. Bombin and M. A. Martin-Delgado, "Topological quantum distillation," Physical Review Letters 97, 180501 (2006). [arXiv:quant-ph/0605138](https://arxiv.org/abs/quant-ph/0605138).
- **[Bombin & Martin-Delgado 2007]** H. Bombin and M. A. Martin-Delgado, "Optimal resources for topological two-dimensional stabilizer codes: Comparative study," Physical Review A 76, 012305 (2007). [arXiv:quant-ph/0703272](https://arxiv.org/abs/quant-ph/0703272).
- **[Bravyi & Kitaev 1998]** S. B. Bravyi and A. Yu. Kitaev, "Quantum codes on a lattice with boundary" (1998). [arXiv:quant-ph/9811052](https://arxiv.org/abs/quant-ph/9811052).
- **[Bravyi & Kitaev 2005]** S. Bravyi and A. Kitaev, "Universal quantum computation with ideal Clifford gates and noisy ancillas," Physical Review A 71, 022316 (2005). [arXiv:quant-ph/0403025](https://arxiv.org/abs/quant-ph/0403025).
- **[Calderbank & Shor 1996]** A. R. Calderbank and P. W. Shor, "Good quantum error-correcting codes exist," Physical Review A 54, 1098 (1996). [arXiv:quant-ph/9512032](https://arxiv.org/abs/quant-ph/9512032).
- **[Campbell 2016]** E. T. Campbell, "The smallest interesting colour code," blog post (2016), https://earltcampbell.com/2016/09/26/the-smallest-interesting-colour-code/.
- **[CRSS 1998]** A. R. Calderbank, E. M. Rains, P. W. Shor, and N. J. A. Sloane, "Quantum error correction via codes over GF(4)," IEEE Transactions on Information Theory 44, 1369 (1998). [arXiv:quant-ph/9608006](https://arxiv.org/abs/quant-ph/9608006).
- **[Dennis et al. 2002]** E. Dennis, A. Kitaev, A. Landahl, and J. Preskill, "Topological quantum memory," Journal of Mathematical Physics 43, 4452 (2002). [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143).
- **[DiVincenzo & Shor 1996]** D. P. DiVincenzo and P. W. Shor, "Fault-tolerant error correction with efficient quantum codes," Physical Review Letters 77, 3260 (1996). [arXiv:quant-ph/9605031](https://arxiv.org/abs/quant-ph/9605031).
- **[Error Correction Zoo]** V. V. Albert, P. Faist, and others (eds.), *The Error Correction Zoo*, https://errorcorrectionzoo.org (2024).
- **[Fowler et al. 2012]** A. G. Fowler, M. Mariantoni, J. M. Martinis, and A. N. Cleland, "Surface codes: Towards practical large-scale quantum computation," Physical Review A 86, 032324 (2012). [arXiv:1208.0928](https://arxiv.org/abs/1208.0928).
- **[Golay 1949]** M. J. E. Golay, "Notes on digital coding," Proceedings of the IRE 37, 657 (1949).
- **[Gottesman 1996]** D. Gottesman, "Class of quantum error-correcting codes saturating the quantum Hamming bound," Physical Review A 54, 1862 (1996). [arXiv:quant-ph/9604038](https://arxiv.org/abs/quant-ph/9604038).
- **[Gottesman 1997]** D. Gottesman, "Stabilizer Codes and Quantum Error Correction," PhD thesis, California Institute of Technology (1997). [arXiv:quant-ph/9705052](https://arxiv.org/abs/quant-ph/9705052).
- **[Grassl 2007]** M. Grassl, "Bounds on the minimum distance of linear codes and quantum codes," online at http://www.codetables.de (2007; accessed 2026).
- **[Grassl et al. 1997]** M. Grassl, T. Beth, and T. Pellizzari, "Codes for the quantum erasure channel," Physical Review A 56, 33 (1997). [arXiv:quant-ph/9610042](https://arxiv.org/abs/quant-ph/9610042).
- **[Horsman et al. 2012]** C. Horsman, A. G. Fowler, S. Devitt, and R. Van Meter, "Surface code quantum computing by lattice surgery," New Journal of Physics 14, 123011 (2012). [arXiv:1111.4022](https://arxiv.org/abs/1111.4022).
- **[Kitaev 2003]** A. Yu. Kitaev, "Fault-tolerant quantum computation by anyons," Annals of Physics 303, 2 (2003). [arXiv:quant-ph/9707021](https://arxiv.org/abs/quant-ph/9707021).
- **[Kubica et al. 2015]** A. Kubica, B. Yoshida, and F. Pastawski, "Unfolding the color code," New Journal of Physics 17, 083026 (2015). [arXiv:1503.02065](https://arxiv.org/abs/1503.02065).
- **[Laflamme et al. 1996]** R. Laflamme, C. Miquel, J. P. Paz, and W. H. Zurek, "Perfect quantum error correcting code," Physical Review Letters 77, 198 (1996). [arXiv:quant-ph/9602019](https://arxiv.org/abs/quant-ph/9602019).
- **[Nielsen & Chuang 2010]** M. A. Nielsen and I. L. Chuang, *Quantum Computation and Quantum Information*, 10th Anniversary Edition, Cambridge University Press (2010). Chapter 10, quantum error-correction.
- **[Rains 1999]** E. M. Rains, "Quantum codes of minimum distance two," IEEE Transactions on Information Theory 45, 266 (1999). [arXiv:quant-ph/9704043](https://arxiv.org/abs/quant-ph/9704043).
- **[Shor 1995]** P. W. Shor, "Scheme for reducing decoherence in quantum computer memory," Physical Review A 52, R2493 (1995).
- **[Steane 1996a]** A. M. Steane, "Multiple-particle interference and quantum error correction," Proceedings of the Royal Society A 452, 2551 (1996). [arXiv:quant-ph/9601029](https://arxiv.org/abs/quant-ph/9601029).
- **[Steane 1996b]** A. M. Steane, "Simple quantum error-correcting codes," Physical Review A 54, 4741 (1996). [arXiv:quant-ph/9605021](https://arxiv.org/abs/quant-ph/9605021).
- **[Steane 1999]** A. M. Steane, "Quantum Reed-Muller codes," IEEE Transactions on Information Theory 45, 1701 (1999). [arXiv:quant-ph/9608026](https://arxiv.org/abs/quant-ph/9608026).
- **[Tomita & Svore 2014]** Y. Tomita and K. M. Svore, "Low-distance surface codes under realistic quantum noise," Physical Review A 90, 062320 (2014). [arXiv:1404.3747](https://arxiv.org/abs/1404.3747).
- **[Vaidman et al. 1996]** L. Vaidman, L. Goldenberg, and S. Wiesner, "Error prevention scheme with four particles," Physical Review A 54, R1745 (1996). [arXiv:quant-ph/9603031](https://arxiv.org/abs/quant-ph/9603031).
- **[Vasmer & Kubica 2022]** M. Vasmer and A. Kubica, "Morphing quantum codes," PRX Quantum 3, 030319 (2022). [arXiv:2112.01446](https://arxiv.org/abs/2112.01446).
