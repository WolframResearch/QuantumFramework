# QuantumInterferometer: three feature gaps vs. other photonic toolkits

Compared against the established photonic software (Xanadu's Strawberry Fields and The
Walrus, Quandela's Perceval, Piquasso, SOQCS), `QuantumInterferometer` is an exact,
symbolic, single-object calculator and diagram tool. Three capabilities those toolkits
have and it does not, in order of how much they limit real use.

## 1. It cannot draw random outcomes the way an experiment does

A real photonic chip fires photons in and a detector clicks, giving one random pattern
per run; repeat it and a distribution builds up. `QuantumInterferometer` instead computes
the *entire* probability table at once. That table grows explosively with the number of
photons and modes, so the exact approach only works for small setups. Every other toolkit
ships a dedicated routine that produces sample outcomes efficiently, without ever building
the full table, which is how they reach interesting sizes. This is the single biggest
practical gap: as written, "boson sampling" here means scoring small cases exactly, never
generating a stream of samples.

*Refs:* P. Clifford and R. Clifford, "The Classical Complexity of Boson Sampling,"
*SODA 2018*, arXiv:1706.01260 (the standard efficient sampling algorithm the Python
toolkits use); Perceval and Strawberry Fields both expose samplers.

## 2. The core counting is pure Wolfram Language, so it is slow

The heavy arithmetic behind the photon statistics (the permanent and hafnian of a matrix)
is done directly in Wolfram Language. It is correct and, unusually, it works on symbolic
inputs, but it is far slower than the specialized compiled C++ libraries the other toolkits
call. Fine for teaching and symbolic work; a hard ceiling for anything numeric at scale.

*Ref:* B. Gupt, J. Izaac, N. Quesada, "The Walrus," *J. Open Source Softw.* **4**(44),
1705 (2019), doi:10.21105/joss.01705: the compiled hafnian/permanent library the Xanadu
stack relies on.

## 3. It assumes a perfect, lossless device

The constructor only accepts an ideal lossless network and refuses anything else. Real
hardware loses photons, the photons are never perfectly identical, and detectors are
imperfect. Modeling those effects is the whole reason some of the other toolkits exist, and
`QuantumInterferometer` models none of them.

*Refs:* J. Osca and J. Vala, "SOQCS: A Stochastic Optical Quantum Circuit Simulator,"
*SoftwareX* **25**, 101603 (2024), arXiv:2307.06965 (built specifically for loss,
distinguishability, and detector imperfections); Piquasso also models photon loss.

---

*None of these are bugs. They are scope decisions: the tool is positioned as exact and
symbolic-first, not as a high-throughput sampler or a hardware-noise simulator. The two
things it does that the Python toolkits do not are the cosine-sine mesh decomposition and
full symbolic evaluation.*
