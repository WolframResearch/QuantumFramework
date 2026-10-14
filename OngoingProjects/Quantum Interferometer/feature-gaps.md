# QuantumInterferometer: feature gaps vs. other photonic toolkits

Compared against the established photonic software (Xanadu's Strawberry Fields and The
Walrus, Quandela's Perceval, Piquasso, SOQCS), `QuantumInterferometer` is an exact,
symbolic, single-object calculator and diagram tool. Below are the capabilities those
toolkits have and it does not, in order of how much they limit real use. Timings were
measured on an 8-core Windows laptop with Wolfram Language 15.0 and The Walrus 0.22.

## 1. It cannot draw random outcomes the way an experiment does

A real photonic chip fires photons in and a detector clicks, giving one random pattern
per run; repeat it and a distribution builds up. `QuantumInterferometer` can score one
output pattern (`"Probability"`) or build the whole probability table (`"Probabilities"`),
but it cannot generate samples. The table grows explosively with the number of photons
and modes, so it is already the bottleneck at six photons:

| photons in 2n modes | output patterns | full table | one sample (prototype) |
|---|---|---|---|
| 4 | 330 | 0.036 s | 1.5 ms |
| 6 | 12,376 | 1.26 s | 3.2 ms |
| 10 | 2 × 10⁷ | out of reach | 13 ms |
| 14 | 3.5 × 10¹⁰ | out of reach | 0.14 s |

The prototype in the last column is the Clifford & Clifford algorithm in about ten lines of
Wolfram Language on top of the built-in `Permanent`. On 3 photons in 6 modes, 10⁵ samples
land 0.008 in total-variation distance from the exact distribution, which is the expected
sampling noise. It also skips the algorithm's Laplace-expansion trick, which would save
another factor of about n. This gap is the largest, and also the cheapest to close.

*Refs:* P. Clifford and R. Clifford, "The Classical Complexity of Boson Sampling,"
*SODA 2018*, arXiv:1706.01260 (the standard efficient sampling algorithm the Python
toolkits use); Perceval and Strawberry Fields both expose samplers.

## 2. It assumes a perfect, lossless device

The constructor only accepts a unitary, and the only elements are beam splitters, phase
shifters and squeezers. Real hardware loses photons, the photons are never perfectly
identical, and detectors are imperfect. These are three separate problems:

- **Loss** can already be modeled by hand: a beam splitter into an extra vacuum mode,
  traced out afterwards. A built-in loss element would make it routine.
- **Partial distinguishability** needs different mathematics (sums of permanents weighted
  by the photons' overlaps), not a new element.
- **Detector imperfections** (efficiency, dark counts) act on the output distribution.

*Refs:* J. Osca and J. Vala, "SOQCS: A Stochastic Optical Quantum Circuit Simulator,"
*SoftwareX* (2024), arXiv:2307.06965 (built specifically for loss, distinguishability, and
detector imperfections); Piquasso also models photon loss.

## 3. Gaussian boson sampling stops short of real experiments

The squeezer support computes photon-number probabilities with hafnians, but two
ingredients of experimental Gaussian boson sampling are missing:

- **Threshold detectors.** Large experiments use click/no-click detectors, whose
  probabilities are torontonians, not hafnians.
- **Displacement.** With no displacement element there are no displaced Gaussian states,
  and so no loop hafnian.

*Refs:* N. Quesada, J. M. Arrazola, N. Killoran, "Gaussian boson sampling using threshold
detectors," *Phys. Rev. A* **98**, 062322 (2018); A. Björklund, B. Gupt, N. Quesada,
arXiv:1805.12498 (hafnian and loop-hafnian algorithms). The Walrus provides both.

## 4. The hafnian and permanent are 20–30× slower than The Walrus

On the same random complex matrices:

| size | hafnian, The Walrus | hafnian, ours | permanent, The Walrus | `Permanent` |
|---|---|---|---|---|
| 16 | 0.91 ms | 15 ms | 2.3 ms | 27 ms |
| 24 | 18 ms | 0.47 s | 0.85 s | 20 s |
| 28 | 0.11 s | 2.4 s | 16 s | 419 s |
| 32 | 0.76 s | 13.7 s | | |

Both libraries grow at the same exponential rate, so this is a constant factor, not a
ceiling, and the values agree to 3 × 10⁻¹¹ or better. The permanent is the built-in
`Permanent`; the hafnian is our own Wolfram Language implementation of the
Björklund–Gupt–Quesada trace formula.

Compiling the hafnian's loop over subsets with `Compile`, run in parallel across 8 threads,
gives 2–9×, with the gain shrinking as the matrices grow. At large sizes the time goes to
arithmetic, not interpretation: each subset needs m matrix powers, O(m·n³) work. The next
step is algorithmic: get the power traces from the characteristic polynomial in O(n³) per
subset, as The Walrus does. For photon-number patterns with repeated modes, a
repeated-row hafnian avoids expanding the matrix at all.

**Status (2026-10-14):**

- ✅ **Completed: compiled hafnian (step 1).** Machine-precision matrices above 8×8 now go
  through a parallel `Compile` of the trace formula, compiled on first use. Exact, symbolic
  and arbitrary-precision input keep the uncompiled code. Measured against the uncompiled
  version in the same kernel: 14× at 20×20, 5× at 24×24, 3× at 28×28, 2× at 32×32, with
  agreement to 10⁻¹⁰ or better. Gotcha: `Tr` inside `Compile` falls back to the evaluator and
  serializes the parallel loop, so the diagonal is summed by hand. Pinning MKL to one thread
  made no difference.
- ~~**Next: characteristic polynomial (step 2).**~~ **Tried, not adopted.** Prototyped twice:
  Householder reduction to Hessenberg form, the La Budde recurrence for det(I − ηC), then
  the η^m coefficient of det(I − ηC)^(−1/2).
  - As a `Compile` (WVM) function, it was 8× *slower* at 16×16, broke even at 24×24, and
    was 1.9× faster at 32×32. The WVM interprets the scalar loops, while step 1's `Dot`
    calls go to BLAS.
  - As a `FunctionCompile` (native) function on one thread, it was 3× slower at 16–20,
    broke even at 24, and was 2.1× faster at 32. Compiling it takes 120 s, which rules
    out compiling on first use.
  - Both versions lose 1–2 digits against step 1: 2 × 10⁻¹¹ relative error versus a 40-digit
    reference at 20×20, where step 1 gets 10⁻¹³.
  - Worth revisiting only with a precompiled native library run on several threads. That
    could plausibly match The Walrus at 32×32. At those sizes, sampling (gap 1) is the
    better tool anyway.

*Ref:* B. Gupt, J. Izaac, N. Quesada, "The Walrus," *J. Open Source Softw.* **4**(44),
1705 (2019), doi:10.21105/joss.01705.

---

*None of these are bugs. They are scope decisions: the tool is positioned as exact and
symbolic-first, not as a high-throughput sampler or a hardware-noise simulator. What it does
that the Python toolkits do not is exact and symbolic Fock amplitudes and probabilities in
one object, from the decomposition to the photon statistics; Perceval also accepts
symbolic circuit parameters, but stops at the unitary. It also offers a cosine-sine mesh
decomposition, which is numeric only.*
