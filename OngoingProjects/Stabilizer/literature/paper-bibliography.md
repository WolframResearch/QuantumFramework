# Stabilizer formalism papers for SDK building

##### [**Undermind**](https://undermind.ai)

---

**Research Goal:** Find educational, pedagogical, or application-oriented academic papers on the stabilizer formalism in quantum information, with emphasis on newer arXiv-friendly or tutorial-style papers. Keep the scope broad across the stabilizer formalism, but prioritize papers that would be useful for building an SDK for the Pauli stabilizer formalism, especially papers that explain tableau representations, Clifford gate update rules, and Pauli measurement/update rules.

*Found 33 papers · April 30, 2026 · Estimated coverage of relevant papers: 100%*

## Summary of Results

The SDK-facing literature is organized around a durable core—Heisenberg/stabilizer propagation \[1\], tableau-based simulation with explicit Clifford and Pauli-measurement updates \[2\], and modern high-performance engineering of the same ideas \[3\].

#### Core representations

- **Tableaux/check matrices** remain the central implementation model for qubit SDKs: Clifford gates act by Pauli conjugation \[1\], and the Aaronson–Gottesman tableau makes this operational for simulation, mixed states, and measurement updates \[2\].
- **Binary symplectic / GF(2) formulations** expose the algebra under the tableau rules and are useful for serialization, verification, and synthesis \[4\], \[5\].
- **Graph-state representations** are an alternative when measurement structure is primary, especially for MBQC-style manipulations \[6\], \[7\].

#### SDK-relevant algorithmic patterns

- **Measurement is the key branching point**: \[2\] is the canonical source for deterministic vs random Pauli measurements and row/phase update rules; \[3\] accelerates deterministic measurement by tracking an inverse tableau.
- **Tableaux as executable specifications** appear in compilation/synthesis work, where matching instantaneous Pauli conjugation is the main invariant \[8\], \[9\].
- **Pauli-frame propagation** connects tableau methods to software architecture for fault-tolerant execution and scheduling \[10\], \[11\], \[12\].

#### Pedagogical entry points

- Broad stabilizer/QEC introductions: \[13\], \[14\], \[15\].
- Newer tutorial-style or accessible material: \[16\], \[17\], \[18\].

## Paper Catalog (33 papers)

|  | Year | Cit/yr | Title | Authors | Journal |
|---:|:--:|:--:|:---|:---|:---|
| 1 | 2021 | 106 | Stim: a fast stabilizer circuit simulator ([link](https://doi.org/10.22331/q-2021-07-06-497)) | C. Gidney | Quantum |
| 2 | 2004 | 69 | Improved Simulation of Stabilizer Circuits ([link](https://doi.org/10.1103/PhysRevA.70.052328)) | S. Aaronson and D. Gottesman | ArXiv |
| 3 | 1998 | 51 | The Heisenberg Representation of Quantum Computers ([link](https://www.semanticscholar.org/paper/4c7e6a179e0f096553c9388b372f021a5f457966)) | D. Gottesman |  |
| 4 | 2003 | 8.5 | Clifford group, stabilizer states, and linear and quadratic operations over GF(2) ([link](https://doi.org/10.1103/PhysRevA.68.042318)) | J. Dehaene and B. Moor | Physical Review A |
| 5 | 2005 | 7.9 | Fast simulation of stabilizer circuits using a graph-state representation ([link](https://doi.org/10.1103/PhysRevA.73.022334)) | S. Anders and H. Briegel | Physical Review A |
| 6 | 2008 | 0.8 | Graphical description of Pauli measurements on stabilizer states ([link](https://doi.org/10.1088/1751-8113/43/2/025301)) | M. Elliott, Bryan Eastin, and C. Caves | Journal of Physics A: Mathematical and Theoretical |
| 7 | 2017 | 1.5 | Discrete Wigner Function Derivation of the Aaronson-Gottesman Tableau Algorithm ([link](https://doi.org/10.3390/e19070353)) | L. Kocia, Yifei Huang, and P. Love | Entropy |
| 8 | 2024 |  | A Simple Method for Compiling Quantum Stabilizer Circuits ([link](https://doi.org/10.1109/QCE60285.2024.00135)) | B. Reid | 2024 IEEE International Conference on Quantum Computing and Engineering (QCE) |
| 9 | 2023 | 2.0 | Fast algorithms for classical specifications of stabiliser states and Clifford gates ([link](https://doi.org/10.22331/q-2025-01-08-1586)) | Nadish de Silva, Wilfred Salmon, and Ming Yin | ArXiv |
| 10 | 2025 |  | Quantum Circuit Optimization and MBQC Scheduling With a Pauli Tracking Library ([link](https://doi.org/10.1109/TQE.2025.3610112)) | Jannis Ruh and Simon Devitt | IEEE Transactions on Quantum Engineering |
| 11 | 2014 | 1.9 | Software-based Pauli tracking in fault-tolerant quantum circuits ([link](https://doi.org/10.7873/DATE.2014.137)) | A. Paler, S. Devitt, K. Nemoto, and I. Polian | 2014 Design, Automation & Test in Europe Conference & Exhibition (DATE) |
| 12 | 2017 | 5.2 | Pauli frames for quantum computer architectures ([link](https://doi.org/10.1145/3061639.3062300)) | L. Riesebos, Xiang Fu, Savvas Varsamopoulos, C. G. Almudéver, and K. Bertels | 2017 54th ACM/EDAC/IEEE Design Automation Conference (DAC) |
| 13 | 2026 |  | A Graphical Rule Book for Clifford Manipulations of Stabilizer States ([link](https://doi.org/10.1109/TQE.2026.3653200)) | A. Patil and Saikat Guha | IEEE Transactions on Quantum Engineering |
| 14 | 2023 | 3.3 | Clifford Manipulations of Stabilizer States: A graphical rule book for Clifford unitaries and measurements on cluster states, and application to photonic quantum computing ([link](https://www.semanticscholar.org/paper/7d5acb74121c2484692160e3dbaff93ee1b23231)) | A. Patil and Saikat Guha |  |
| 15 | 2025 | 1.0 | A Streamlined Demonstration That Stabilizer Circuits Simulation Reduces to Boolean Linear Algebra ([link](https://doi.org/10.1134/S1995080225607842)) | V. Yashin | Lobachevskii Journal of Mathematics |
| 16 | 2023 | 2.3 | Architecture-Aware Synthesis of Stabilizer Circuits from Clifford Tableaus ([link](https://www.semanticscholar.org/paper/5858a45a548ec4f5db68ea489062665104ae2bf7)) | David Winderl, Qunsheng Huang, A. M. D. Griend, and Richie Yeung |  |
| 17 | 2014 | 8.2 | How to efficiently select an arbitrary Clifford group element ([link](https://doi.org/10.1063/1.4903507)) | Robert Koenig and J. Smolin | Journal of Mathematical Physics |
| 18 | 2004 | 6.4 | Stabilizer states and Clifford operations for systems of arbitrary dimensions and modular arithmetic ([link](https://doi.org/10.1103/PhysRevA.71.042315)) | E. Hostens, J. Dehaene, and B. Moor | Physical Review A |
| 19 | 2011 | 2.9 | A linearized stabilizer formalism for systems of finite dimension ([link](https://doi.org/10.26421/QIC13.1-2-6)) | N. D. Beaudrap | Quantum Inf. Comput. |
| 20 | 1997 | 101 | Stabilizer Codes and Quantum Error Correction ([link](https://doi.org/10.7907/RZR7-DT72.)) | D. Gottesman | arXiv: Quantum Physics |
| 21 | 2000 | 5.5 | An Introduction to Quantum Error Correction ([link](https://doi.org/10.1090/psapm/058/1922900)) | Daniel Gottesman | arXiv: Quantum Physics |
| 22 | 2009 | 39 | An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation ([link](https://doi.org/10.1090/psapm/068/2762145)) | D. Gottesman | arXiv: Quantum Physics |
| 23 | 2023 | 3.5 | Quantum Circuits for Stabilizer Error Correcting Codes: A Tutorial ([link](https://doi.org/10.1109/MCAS.2024.3349668)) | Arijit Mondal and Keshab K. Parhi | IEEE Circuits and Systems Magazine |
| 24 | 2024 | 0.5 | On Classical Simulation of Quantum Circuits Composed of Clifford Gates ([link](https://doi.org/10.12743/quanta.v13i1.265)) | G. Biswas | Quanta |
| 25 | 2015 | 0.2 | Stabilizer Formalism and Its Applications ([link](https://doi.org/10.1007/978-981-287-996-7_2)) | K. Fujii |  |
| 26 | 2014 |  | 1 Quantum Error Correction and The Stabilizer Formalism ( Revision ) ([link](https://www.semanticscholar.org/paper/9ca582bd85bf395daadeebf658449d9e458a9c59)) | D. Browne |  |
| 27 | 2024 | 2.5 | Lecture notes on quantum entanglement: From stabilizer states to stabilizer channels ([link](https://doi.org/10.1007/s11467-024-1397-4)) | Amir R. Arab | Frontiers of Physics |
| 28 | 2024 | 0.6 | Qudit Quantum Programming with Projective Cliffords ([link](https://doi.org/10.1145/3776646)) | Jennifer Paykin and Sam Winnick | Proceedings of the ACM on Programming Languages |
| 29 | 2024 | 0.6 | Condensed encodings of projective Clifford operations in arbitrary dimension ([link](https://doi.org/10.1063/5.0239974)) | Sam Winnick and Jennifer Paykin | Journal of Mathematical Physics |
| 30 | 2026 |  | PauliEngine: High-Performant Symbolic Arithmetic for Quantum Operations ([link](https://doi.org/10.48550/arXiv.2601.02233)) | Leon Müller, Adelina Bärligea, Alexander Knapp, and Jakob S. Kottmann | ArXiv |
| 31 | 2023 | 2.0 | SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits ([link](https://doi.org/10.1145/3649329.3655902)) | Wang Fang and Mingsheng Ying | 2024 61st ACM/IEEE Design Automation Conference (DAC) |
| 32 | 2015 | 2.0 | Simulation of Quantum Circuits via Stabilizer Frames ([link](https://doi.org/10.1109/TC.2014.2360532)) | Héctor J. García and I. Markov | IEEE Transactions on Computers |
| 33 | 2012 | 2.1 | Efficient Inner-product Algorithm for Stabilizer States ([link](https://www.semanticscholar.org/paper/8bb6b59914487f634d3f840c66ab4dfa1882b657)) | Héctor J. García, I. Markov, and Andrew W. Cross | ArXiv |

### Paper Details

1\. · 100% match · 2021 · 106 cit/yr\
**Stim: a fast stabilizer circuit simulator** ([link](https://doi.org/10.22331/q-2021-07-06-497))\
C. Gidney\
*Quantum* · Mar 3, 2021 · 547 citations

> This paper presents “Stim”, a fast simulator for quantum stabilizer circuits. The paper explains how Stim works and compares it to existing tools. With no foreknowledge, Stim can analyze a distance 100 surface code circuit (20 thousand qubits, 8 million gates, 1 million measurements) in 15 seconds and then begin sampling full circuit shots at a rate of 1 kHz. Stim uses a stabilizer tableau representation, similar to Aaronson and Gottesman’s CHP simulator, but with three main improvements. First, Stim improves the asymptotic complexity of deterministic measurement from quadratic to linear by tracking the inverse of the circuit’s stabilizer tableau. Second, Stim improves the constant factors of the algorithm by using a cache-friendly data layout and 256 bit wide SIMD instructions. Third, Stim only uses expensive stabilizer tableau simulation to create an initial reference sample. Further samples are collected in bulk by using that sample as a reference for batches of Pauli frames propagating through the circuit.

------------------------------------------------------------------------

2\. · 100% match · 2004 · 69 cit/yr\
**Improved Simulation of Stabilizer Circuits** ([link](https://doi.org/10.1103/PhysRevA.70.052328))\
S. Aaronson and D. Gottesman\
*ArXiv* · Jun 25, 2004 · 1510 citations

> The Gottesman-Knill theorem says that a stabilizer circuit\char22{}that is, a quantum circuit consisting solely of controlled-NOT (CNOT), Hadamard, and phase gates\char22{}can be simulated efficiently on a classical computer. This paper improves that theorem in several directions. First, by removing the need for Gaussian elimination, we make the simulation algorithm much faster at the cost of a factor of 2 increase in the number of bits needed to represent a state. We have implemented the improved algorithm in a freely available program called CHP (CNOT-Hadamard-phase), which can handle thousands of qubits easily. Second, we show that the problem of simulating stabilizer circuits is complete for the classical complexity class $`\ensuremath{\bigoplus}\mathsf{L}`$, which means that stabilizer circuits are probably not even universal for classical computation. Third, we give efficient algorithms for computing the inner product between two stabilizer states, putting any $`n`$-qubit stabilizer circuit into a \`\`canonical form’’ that requires at most $`O({n}^{2}∕\mathrm{log}\phantom{\rule{0.2em}{0ex}}n)`$ gates, and other useful tasks. Fourth, we extend our simulation algorithm to circuits acting on mixed states, circuits containing a limited number of nonstabilizer gates, and circuits acting on general tensor-product initial states but containing only a limited number of measurements.

------------------------------------------------------------------------

3\. · 100% match · 1998 · 51 cit/yr\
**The Heisenberg Representation of Quantum Computers** ([link](https://www.semanticscholar.org/paper/4c7e6a179e0f096553c9388b372f021a5f457966))\
D. Gottesman\
Jun 24, 1998 · 1433 citations

> Since Shor\`s discovery of an algorithm to factor numbers on a quantum computer in polynomial time, quantum computation has become a subject of immense interest. Unfortunately, one of the key features of quantum computers–the difficulty of describing them on classical computers–also makes it difficult to describe and understand precisely what can be done with them. A formalism describing the evolution of operators rather than states has proven extremely fruitful in understanding an important class of quantum operations. States used in error correction and certain communication protocols can be described by their stabilizer, a group of tensor products of Pauli matrices. Even this simple group structure is sufficient to allow a rich range of quantum effects, although it falls short of the full power of quantum computation.

------------------------------------------------------------------------

4\. · 100% match · 2003 · 8.5 cit/yr\
**Clifford group, stabilizer states, and linear and quadratic operations over GF(2)** ([link](https://doi.org/10.1103/PhysRevA.68.042318))\
J. Dehaene and B. Moor\
*Physical Review A* · Apr 18, 2003 · 196 citations

> We describe stabilizer states and Clifford group operations using linear operations and quadratic forms over binary vector spaces. We show how the n-qubit Clifford group is isomorphic to a group with an operation that is defined in terms of a (2n+1)x(2n+1) binary matrix product and binary quadratic forms. As an application we give two schemes to efficiently decompose Clifford group operations into one- and two-qubit operations. We also show how the coefficients of stabilizer states and Clifford group operations in a standard basis expansion can be described by binary quadratic forms. Our results are useful for quantum error correction, entanglement distillation, and possibly quantum computing.

------------------------------------------------------------------------

5\. · 100% match · 2005 · 7.9 cit/yr\
**Fast simulation of stabilizer circuits using a graph-state representation** ([link](https://doi.org/10.1103/PhysRevA.73.022334))\
S. Anders and H. Briegel\
*Physical Review A* · Apr 15, 2005 · 167 citations

> According to the Gottesman-Knill theorem, a class of quantum circuits\char22{}namely, the so-called stabilizer circuits\char22{}can be simulated efficiently on a classical computer. We introduce an algorithm for this task, which is based on the graph-state formalism. It shows significant improvement in comparison to an existing algorithm, given by Gottesman and Aaronson, in terms of speed and of the number of qubits the simulator can handle. We also present an implementation.

------------------------------------------------------------------------

6\. · 100% match · 2008 · 0.8 cit/yr\
**Graphical description of Pauli measurements on stabilizer states** ([link](https://doi.org/10.1088/1751-8113/43/2/025301))\
M. Elliott, Bryan Eastin, and C. Caves\
*Journal of Physics A: Mathematical and Theoretical* · Jun 16, 2008 · 15 citations

> We use a graphical representation of stabilizer states to describe, simply and efficiently, the effect of measurements of Pauli products on stabilizer states. This work complements our earlier work (2008 Phys. Rev. A 77 042307), which described in graphical terms the action of Clifford operations on stabilizer states.

------------------------------------------------------------------------

7\. · 100% match · 2017 · 1.5 cit/yr\
**Discrete Wigner Function Derivation of the Aaronson-Gottesman Tableau Algorithm** ([link](https://doi.org/10.3390/e19070353))\
L. Kocia, Yifei Huang, and P. Love\
*Entropy* · Mar 14, 2017 · 14 citations

> The Gottesman–Knill theorem established that stabilizer states and Clifford operations can be efficiently simulated classically. For qudits with odd dimension three and greater, stabilizer states and Clifford operations have been found to correspond to positive discrete Wigner functions and dynamics. We present a discrete Wigner function-based simulation algorithm for odd-d qudits that has the same time and space complexity as the Aaronson–Gottesman algorithm for qubits. We show that the efficiency of both algorithms is due to harmonic evolution in the symplectic structure of discrete phase space. The differences between the Wigner function algorithm for odd-d and the Aaronson–Gottesman algorithm for qubits are likely due only to the fact that the Weyl–Heisenberg group is not in S U ( d ) for d = 2 and that qubits exhibit state-independent contextuality. This may provide a guide for extending the discrete Wigner function approach to qubits.

------------------------------------------------------------------------

8\. · 100% match · 2024\
**A Simple Method for Compiling Quantum Stabilizer Circuits** ([link](https://doi.org/10.1109/QCE60285.2024.00135))\
B. Reid\
*2024 IEEE International Conference on Quantum Computing and Engineering (QCE)* · Apr 30, 2024 · 0 citations

> Stabilizer circuits play an important role in quantum error correction protocols, and will be vital for ensuring fault tolerance in future quantum hardware. While stabilizer circuits are defined on the Clifford generating set, $`\{H, S, CX\}`$, not all of these gates are native to quantum hardware. As such they must be compiled into the native gateset, with the key difference across hardware archetypes being the native two-qubit gate. Here we introduce an intuitive and accessible method for Clifford gate compilation. While multiple open source solutions exist for quantum circuit compilation, these operate on arbitrary quantum gates. By restricting ourselves to Clifford gates, the compilation process becomes almost trivial and even large circuits can be compiled manually. The core idea is well known: if two Clifford circuits conjugate Paulis identically, they are equivalent. Compilation is then reduced to ensuring that the instantaneous Pauli conjugation is correct for each qubit at every timestep. This is Tableaux Manipulation, so called as we directly interrogate stabilizer tableaux to ensure correct Pauli conjugation. We provide a brief explanation of the process along with a worked example to build intuition; we finally show some comparisons for compiling large circuits to open source software, and highlight that this method ensures a minimal number of quantum gates are employed.

------------------------------------------------------------------------

9\. · 98% match · 2023 · 2.0 cit/yr\
**Fast algorithms for classical specifications of stabiliser states and Clifford gates** ([link](https://doi.org/10.22331/q-2025-01-08-1586))\
Nadish de Silva, Wilfred Salmon, and Ming Yin\
*ArXiv* · Nov 17, 2023 · 5 citations

> The stabiliser formalism plays a central role in quantum computing, error correction, and fault tolerance. Conversions between and verifications of different specifications of stabiliser states and Clifford gates are important components of many classical algorithms in quantum information, e.g. for gate synthesis, circuit optimisation, and simulating quantum circuits. These core functions are also used in the numerical experiments critical to formulating and testing mathematical conjectures on the stabiliser formalism.We develop novel mathematical insights concerning stabiliser states and Clifford gates that significantly clarify their descriptions. We then utilise these to provide ten new fast algorithms which offer asymptotic advantages over any existing implementations. We show how to rapidly verify that a vector is a stabiliser state, and interconvert between its specification as amplitudes, a quadratic form, and a check matrix. These methods are leveraged to rapidly check if a given unitary matrix is a Clifford gate and to interconvert between the matrix of a Clifford gate and its compact specification as a stabiliser tableau.For example, we extract the stabiliser tableau of a 2n×2n matrix, promised to be a Clifford gate, in O(n2n) time. Remarkably, it is not necessary to read all the elements of a Clifford gate matrix to extract its stabiliser tableau. This is an asymptotic speedup over the best-known method that is exponential in the number of qubits.We provide implementations of our algorithms in Python and C++ that exhibit vastly improved practical performance over existing algorithms in the cases where they exist.

------------------------------------------------------------------------

10\. · 92% match · 2025\
**Quantum Circuit Optimization and MBQC Scheduling With a Pauli Tracking Library** ([link](https://doi.org/10.1109/TQE.2025.3610112))\
Jannis Ruh and Simon Devitt\
*IEEE Transactions on Quantum Engineering* · 0 citations

> In this article, we present a generic framework for the scheduling of qubits in measurement-based quantum computing and error-corrected circuits that are implemented through Clifford circuits. The framework is based on the commutation of Pauli operators through quantum circuits, which is called Pauli tracking. We complement the framework by providing an independent software library to perform the Pauli tracking, specifically in the context of measurement-based quantum computing. Tracking Pauli operators allows one to reduce the number of Pauli gates that must be executed on quantum hardware and to capture the constraints on the order of measurements. The scheduling problem is numerically investigated on small scale.

------------------------------------------------------------------------

11\. · 88% match · 2014 · 1.9 cit/yr\
**Software-based Pauli tracking in fault-tolerant quantum circuits** ([link](https://doi.org/10.7873/DATE.2014.137))\
A. Paler, S. Devitt, K. Nemoto, and I. Polian\
*2014 Design, Automation & Test in Europe Conference & Exhibition (DATE)* · Mar 24, 2014 · 23 citations

------------------------------------------------------------------------

12\. · 84% match · 2017 · 5.2 cit/yr\
**Pauli frames for quantum computer architectures** ([link](https://doi.org/10.1145/3061639.3062300))\
L. Riesebos, Xiang Fu, Savvas Varsamopoulos, C. G. Almudéver, and K. Bertels\
*2017 54th ACM/EDAC/IEEE Design Automation Conference (DAC)* · Jun 18, 2017 · 46 citations

------------------------------------------------------------------------

13\. · 83% match · 2026\
**A Graphical Rule Book for Clifford Manipulations of Stabilizer States** ([link](https://doi.org/10.1109/TQE.2026.3653200))\
A. Patil and Saikat Guha\
*IEEE Transactions on Quantum Engineering* · 0 citations

> Stabilizer states, along with Clifford manipulations (unitary transformations and measurements) thereof—despite being efficiently simulable on a classical computer—are an important tool in quantum information processing, with applications to quantum computing, error correction, and networking. Graph states, defined on a graph, are a special class of stabilizer states that are central to measurement-based quantum computing, all-photonic quantum repeaters, distributed quantum computing, and entanglement distribution in a network. All stabilizer states are local-Clifford equivalent to graph states. In this article, we review the stabilizer framework and extend it by incorporating general stabilizer measurements such as multiqubit joint projections. We provide an explicit procedure—using Karnaugh maps from Boolean algebra—for converting arbitrary stabilizer gates into tableau operations of the cnot–Hadamard–Phase formalism for efficient stabilizer manipulations. We derive graphical rules for arbitrary stabilizer manipulations of graph states, including multiqubit stabilizer projections and unitaries. We implement the graphical rulebook resulting from above into a MATLAB simulator with a graphical user interface. A user of this tool, e.g., for research in quantum networks, will not require any background in quantum information or the stabilizer framework.

------------------------------------------------------------------------

14\. · 78% match · 2023 · 3.3 cit/yr\
**Clifford Manipulations of Stabilizer States: A graphical rule book for Clifford unitaries and measurements on cluster states, and application to photonic quantum computing** ([link](https://www.semanticscholar.org/paper/7d5acb74121c2484692160e3dbaff93ee1b23231))\
A. Patil and Saikat Guha\
Dec 4, 2023 · 8 citations

> Stabilizer states along with Clifford manipulations (unitary transformations and measurements) thereof – despite being efficiently simulable on a classical computer – are an important tool in quantum information processing, with applications to quantum computing, error correction and networking. Cluster states, defined on a graph, are a special class of stabilizer states that are central to measurement based quantum computing, all-photonic quantum repeaters, distributed quantum computing, and entanglement distribution in a network. All cluster states are local-Clifford equivalent to a stabilizer state. In this paper, we review the stabilizer framework, and extend it, by: incorporating general stabilizer measurements such as multi-qubit fusions, and providing an explicit procedure – using Karnaugh maps from Boolean algebra – for converting arbitrary stabilizer gates into tableau operations of the CHP formalism for efficient stabilizer manipulations. Using these tools, we develop a graphical rule-book and a MATLAB simulator with a graphical user interface for arbitrary stabilizer manipulations of cluster states, a user of which, e.g., for research in quantum networks, will not require any background in quantum information or the stabilizer framework. We extend our graphical rule-book to include dual-rail photonic-qubit cluster state manipulations with probabilistically-heralded linear-optical circuits for various rotated Bell measurements, i.e., fusions (including new \`Type-I’ fusions we propose, where only one of the two fused qubits is destructively measured), by incorporating graphical rules for their success and failure modes. Finally, we show how stabilizer descriptions of multi-qubit fusions can be mapped to linear optical circuits.

------------------------------------------------------------------------

15\. · 75% match · 2025 · 1.0 cit/yr\
**A Streamlined Demonstration That Stabilizer Circuits Simulation Reduces to Boolean Linear Algebra** ([link](https://doi.org/10.1134/S1995080225607842))\
V. Yashin\
*Lobachevskii Journal of Mathematics* · Apr 18, 2025 · 1 citations

> Gottesman–Knill theorem states that computations on stabilizer circuits can be simulated on a classical computer, conventional simulation algorithms extensively use linear algebra over bit strings. For instance, given a non-adaptive stabilizer circuit, the problem of computing the probability of a given outcome (strong simulation) is known to be log-space reducible to solving the system of linear equations over Boolean variables, which is commonly done by Gaussian elimination. This note aims to make the connection between stabilizer circuits and Boolean linear algebra even more explicit. To do this, we extend the stabilizer tableau formalism to include stabilizer tableau descriptions of arbitrary stabilizer operations (Clifford channels). Finding the tableau corresponding to the composition of two channels becomes a linear algebra problem. Any stabilizer circuit rewrites to a diagram with stabilizer tableaux on vertices, contracting an edge means to take the composition of channels, to compute the result of the circuit means to fully contract the diagram. Thus, simulating stabilizer circuits reduces to a sequence of Gaussian eliminations. This approach gives a new perspective on explaining the work of stabilizer tableau methods (reproducing the asymptotics) and creates opportunity for exploring various tensor-contraction techniques in stabilizer simulation.

------------------------------------------------------------------------

16\. · 72% match · 2023 · 2.3 cit/yr\
**Architecture-Aware Synthesis of Stabilizer Circuits from Clifford Tableaus** ([link](https://www.semanticscholar.org/paper/5858a45a548ec4f5db68ea489062665104ae2bf7))\
David Winderl, Qunsheng Huang, A. M. D. Griend, and Richie Yeung\
Sep 16, 2023 · 6 citations

> Since quantum computing is currently in the NISQ-Era, compilation strategies to reduce the number of gates executed on specific hardware are required. In this work, we utilize the concept of synthesis of a data structure called Clifford tableaus, focusing on applying CNOTs within the respective connectivity graph of the quantum device. We hence contribute to the field of compilation or, more precisely, synthesis by reducing the number of CNOTs in the synthesized quantum circuit. Upon convergence, our method shows to outperform other state-of-the-art synthesis techniques, when executed with respect to a specific hardware. Upon executing the resulting circuits on real hardware, our synthesized circuits tend to increase the final fidelity and reduce the overall execution times.

------------------------------------------------------------------------

17\. · 66% match · 2014 · 8.2 cit/yr\
**How to efficiently select an arbitrary Clifford group element** ([link](https://doi.org/10.1063/1.4903507))\
Robert Koenig and J. Smolin\
*Journal of Mathematical Physics* · Jun 9, 2014 · 97 citations

> We give an algorithm which produces a unique element of the Clifford group on n qubits ( Cn ) from an integer 0≤i\<Cn (the number of elements in the group). The algorithm involves O(n 3) operations and provides, in addition to a canonical mapping from the integers to group elements g, a factorization of g into a sequence of at most 4n symplectic transvections. The algorithm can be used to efficiently select random elements of Cn which are often useful in quantum information theory and quantum computation. We also give an algorithm for the inverse map, indexing a group element in time O(n 3).

------------------------------------------------------------------------

18\. · 63% match · 2004 · 6.4 cit/yr\
**Stabilizer states and Clifford operations for systems of arbitrary dimensions and modular arithmetic** ([link](https://doi.org/10.1103/PhysRevA.71.042315))\
E. Hostens, J. Dehaene, and B. Moor\
*Physical Review A* · Aug 31, 2004 · 139 citations

> We describe generalizations of the Pauli group, the Clifford group, and stabilizer states for qudits in a Hilbert space of arbitrary dimension d. We examine a link with modular arithmetic, which yields an efficient way of representing the Pauli group and the Clifford group with matrices over Z{sub d}. We further show how a Clifford operation can be efficiently decomposed into one and two-qudit operations. We also focus in detail on standard basis expansions of stabilizer states.

------------------------------------------------------------------------

19\. · 60% match · 2011 · 2.9 cit/yr\
**A linearized stabilizer formalism for systems of finite dimension** ([link](https://doi.org/10.26421/QIC13.1-2-6))\
N. D. Beaudrap\
*Quantum Inf. Comput.* · Feb 16, 2011 · 44 citations

> The stabilizer formalism is a scheme, generalizing well-known techniques developed by Gottesman \[1\] in the case of qubits, to efficiently simulate a class of transformations (stabilizer circuits, which include the quantum Fourier transform and highly entangling operations) on standard basis states of d-dimensional qudits. To determine the state of a simulated system, existing treatments involve the computation of cumulative phase factors which involve quadratic dependencies. We present a simple formalism in which Pauli operators are represented using displacement operators in discrete phase space, expressing the evolution of the state via linear transformations modulo D ≤ 2d. We thus obtain a simple proof that simulating stabilizer circuits on n qudits, involving any constant number of measurement rounds, is complete for the complexity class coModdL and may be simulated by O(log(n)2)-depth circuits for any constant d ≥ 2.

------------------------------------------------------------------------

20\. · 58% match · 1997 · 101 cit/yr\
**Stabilizer Codes and Quantum Error Correction** ([link](https://doi.org/10.7907/RZR7-DT72.))\
D. Gottesman\
*arXiv: Quantum Physics* · May 28, 1997 · 2918 citations

> Controlling operational errors and decoherence is one of the major challenges facing the field of quantum computation and other attempts to create specified many-particle entangled states. The field of quantum error correction has developed to meet this challenge. A group-theoretical structure and associated subclass of quantum codes, the stabilizer codes, has proved particularly fruitful in producing codes and in understanding the structure of both specific codes and classes of codes. I will give an overview of the field of quantum error correction and the formalism of stabilizer codes. In the context of stabilizer codes, I will discuss a number of known codes, the capacity of a quantum channel, bounds on quantum codes, and fault-tolerant quantum computation.

------------------------------------------------------------------------

21\. · 56% match · 2000 · 5.5 cit/yr\
**An Introduction to Quantum Error Correction** ([link](https://doi.org/10.1090/psapm/058/1922900))\
Daniel Gottesman\
*arXiv: Quantum Physics* · Apr 18, 2000 · 142 citations

> Quantum states are very delicate, so it is likely some sort of quantum error correction will be necessary to build reliable quantum computers. The theory of quantum error-correcting codes has some close ties to and some striking differences from the theory of classical error-correcting codes. Many quantum codes can be described in terms of the stabilizer of the codewords. The stabilizer is a finite Abelian group, and allows a straightforward characterization of the error-correcting properties of the code. The stabilizer formalism for quantum codes also illustrates the relationships to classical coding theory, particularly classical codes over GF(4), the finite field with four elements.

------------------------------------------------------------------------

22\. · 52% match · 2009 · 39 cit/yr\
**An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation** ([link](https://doi.org/10.1090/psapm/068/2762145))\
D. Gottesman\
*arXiv: Quantum Physics* · Apr 16, 2009 · 667 citations

> Quantum states are very delicate, so it is likely some sort of quantum error correction will be necessary to build reliable quantum computers. The theory of quantum error-correcting codes has some close ties to and some striking differences from the theory of classical error-correcting codes. Many quantum codes can be described in terms of the stabilizer of the codewords. The stabilizer is a finite Abelian group, and allows a straightforward characterization of the error-correcting properties of the code. The stabilizer formalism for quantum codes also illustrates the relationships to classical coding theory, particularly classical codes over GF(4), the finite field with four elements. To build a quantum computer which behaves correctly in the presence of errors, we also need a theory of fault-tolerant quantum computation, instructing us how to perform quantum gates on qubits which are encoded in a quantum error-correcting code. The threshold theorem states that it is possible to create a quantum computer to perform an arbitrary quantum computation provided the error rate per physical gate or time step is below some constant threshold value.

------------------------------------------------------------------------

23\. · 50% match · 2023 · 3.5 cit/yr\
**Quantum Circuits for Stabilizer Error Correcting Codes: A Tutorial** ([link](https://doi.org/10.1109/MCAS.2024.3349668))\
Arijit Mondal and Keshab K. Parhi\
*IEEE Circuits and Systems Magazine* · Sep 21, 2023 · 9 citations

> Quantum computers have the potential to provide exponential speedups over their classical counterparts. Quantum principles are being applied to fields such as communications, information processing, and artificial intelligence to achieve quantum advantage. However, quantum bits are extremely noisy and prone to decoherence. Thus, keeping the qubits error free is extremely important toward reliable quantum computing. Quantum error correcting codes have been studied for several decades and methods have been proposed to import classical error correcting codes to the quantum domain. Along with the exploration into novel and more efficient quantum error correction codes, it is also essential to design circuits for practical realization of these codes. This article serves as a tutorial on designing and simulating quantum encoder and decoder circuits for stabilizer codes. We first describe Shor’s 9-qubit code which was the first quantum error correcting code. We discuss the stabilizer formalism along with the design of encoding and decoding circuits for stabilizer codes such as the five-qubit code and Steane code. We also design nearest neighbor compliant circuits for the above codes. The circuits were simulated and verified using IBM Qiskit.

------------------------------------------------------------------------

24\. · 47% match · 2024 · 0.5 cit/yr\
**On Classical Simulation of Quantum Circuits Composed of Clifford Gates** ([link](https://doi.org/10.12743/quanta.v13i1.265))\
G. Biswas\
*Quanta* · May 22, 2024 · 1 citations

> The Gottesman–Knill theorem asserts that quantum circuits composed solely of Clifford gates can be efficiently simulated classically. This theorem hinges on the fact that Clifford gates map Pauli strings to other Pauli strings, thereby allowing for a structured simulation process using classical computations. In this work, we break down the step-by-step procedure of the Gottesman–Knill theorem in a beginner-friendly manner, leveraging concepts such as matrix products, tensor products, commutation, anti-commutation, eigenvalues, and eigenvectors of quantum mechanical operators. Through detailed examples illustrating superposition and entanglement phenomena, we aim to provide a clear understanding of the classical simulation of Clifford gate-based quantum circuits. While we do not provide a formal proof of the theorem, we offer intuitive physical insights at each stage where necessary, empowering readers to grasp the fundamental principles underpinning this intriguing aspect of quantum computation.Quanta 2024; 13: 20–27.

------------------------------------------------------------------------

25\. · 42% match · 2015 · 0.2 cit/yr\
**Stabilizer Formalism and Its Applications** ([link](https://doi.org/10.1007/978-981-287-996-7_2))\
K. Fujii\
2 citations

------------------------------------------------------------------------

26\. · 40% match · 2014\
**1 Quantum Error Correction and The Stabilizer Formalism ( Revision )** ([link](https://www.semanticscholar.org/paper/9ca582bd85bf395daadeebf658449d9e458a9c59))\
D. Browne\
0 citations

------------------------------------------------------------------------

27\. · 38% match · 2024 · 2.5 cit/yr\
**Lecture notes on quantum entanglement: From stabilizer states to stabilizer channels** ([link](https://doi.org/10.1007/s11467-024-1397-4))\
Amir R. Arab\
*Frontiers of Physics* · Apr 16, 2024 · 5 citations

------------------------------------------------------------------------

28\. · 35% match · 2024 · 0.6 cit/yr\
**Qudit Quantum Programming with Projective Cliffords** ([link](https://doi.org/10.1145/3776646))\
Jennifer Paykin and Sam Winnick\
*Proceedings of the ACM on Programming Languages* · Jul 23, 2024 · 1 citations

> This paper introduces a novel abstraction for programming quantum operations, specifically projective Cliffords, as functions over the qudit Pauli group. Generalizing the idea behind Pauli tableaux, we introduce a type system and lambda calculus for projective Cliffords called LambdaPC that captures well-formed Clifford operations via a Curry-Howard correspondence with a particular encoding of the Clifford and Pauli groups. In LambdaPC, users write functions that encode projective Cliffords P ↦ U P U†, and such functions are compiled to circuits executable on modern quantum computers that transform quantum states \|ϕ⟩ into U \|ϕ⟩, up to a global phase. Importantly, the language captures not just qubit operations, but qudit operations for any dimension d. Throughout the paper we explore what it means to program with projective Cliffords through a number of examples and a case study focusing on stabilizer error correcting codes.

------------------------------------------------------------------------

29\. · 32% match · 2024 · 0.6 cit/yr\
**Condensed encodings of projective Clifford operations in arbitrary dimension** ([link](https://doi.org/10.1063/5.0239974))\
Sam Winnick and Jennifer Paykin\
*Journal of Mathematical Physics* · Jul 23, 2024 · 1 citations

> We provide a careful analysis of the structure theorem for the $`n`$-qudit projective Clifford group and various encoding schemes for its elements. In particular, we derive formulas for evaluation, composition, and inversion. Our results apply to all integers $`d\geq2`$, most notably the even case.

------------------------------------------------------------------------

30\. · 31% match · 2026\
**PauliEngine: High-Performant Symbolic Arithmetic for Quantum Operations** ([link](https://doi.org/10.48550/arXiv.2601.02233))\
Leon Müller, Adelina Bärligea, Alexander Knapp, and Jakob S. Kottmann\
*ArXiv* · Jan 5, 2026 · 0 citations

> Quantum computation is inherently hybrid, and fast classical manipulation of qubit operators is necessary to ensure scalability in quantum software. We introduce PauliEngine, a high-performance C++ framework that provides efficient primitives for Pauli string multiplication, commutators, symbolic phase tracking, and structural transformations. Built on a binary symplectic representation and optimized bit-wise operations, PauliEngine supports both numerical and symbolic coefficients and is accessible through a Python interface. Runtime benchmarks demonstrate substantial speedups over state-of-the-art implementations. PauliEngine provides a scalable backend for operator-based quantum software tools and simulations.

------------------------------------------------------------------------

31\. · 30% match · 2023 · 2.0 cit/yr\
**SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits** ([link](https://doi.org/10.1145/3649329.3655902))\
Wang Fang and Mingsheng Ying\
*2024 61st ACM/IEEE Design Automation Conference (DAC)* · Nov 7, 2023 · 5 citations

> This paper proposes an efficient stabilizer circuit simulation algorithm that only traverses the circuit forward once. We introduce phase symbolization into stabilizer generators, which allows possible Pauli faults in the circuit to be accumulated explicitly as symbolic expressions in the phases of stabilizer generators. This way, the measurement outcomes are also symbolic expressions, and we can sample them by substituting the symbolic variables with concrete values, without traversing the circuit repeatedly. We show how to integrate symbolic phases into the stabilizer tableau and maintain them efficiently using bit-vector encoding. A new data layout of the stabilizer tableau in memory is proposed, which improves the performance of our algorithm (and other stabilizer simulation algorithms based on the stabilizer tableau). We implement our algorithm and data layout in a Julia package named SymPhase.jl, and compare it with Stim, the state-of-the-art simulator, on several benchmarks. We show that SymPhase.jl has superior performance in terms of sampling time, which is crucial for generating a large number of samples for further analysis.

------------------------------------------------------------------------

32\. · 29% match · 2015 · 2.0 cit/yr\
**Simulation of Quantum Circuits via Stabilizer Frames** ([link](https://doi.org/10.1109/TC.2014.2360532))\
Héctor J. García and I. Markov\
*IEEE Transactions on Computers* · Aug 1, 2015 · 22 citations

> Generic quantum-circuit simulation appears intractable for conventional computers and may be unnecessary because useful quantum circuits exhibit significant structure that can be exploited during simulation. For example, Gottesman and Knill identified an important subclass, called stabilizer circuits, which can be simulated efficiently using group-theory techniques and insights from quantum physics. Realistic circuits enriched with quantum error-correcting codes and fault-tolerant procedures are dominated by stabilizer subcircuits and contain a relatively small number of non-Clifford components. Therefore, we develop new data structures and algorithms that facilitate parallel simulation of such circuits. Stabilizer frames offer more compact storage than previous approaches but require more sophisticated bookkeeping. Our implementation, called Quipu, simulates certain quantum arithmetic circuits (e.g., reversible ripple-carry adders) in polynomial time and space for equal superpositions of n-qubits. On such instances, known linear-algebraic simulation techniques, such as the (state-of-the-art) BDD-based simulator QuIDDPro, take exponential time. We simulate quantum Fourier transform and quantum fault-tolerant circuits using Quipu, and the results demonstrate that our stabilizer-based technique empirically outperforms QuIDDPro in all cases. While previous high-performance, structure-aware simulations of quantum circuits were difficult to parallellize, we demonstrate that Quipu can be parallelized with a nontrivial computational speedup.

------------------------------------------------------------------------

33\. · 28% match · 2012 · 2.1 cit/yr\
**Efficient Inner-product Algorithm for Stabilizer States** ([link](https://www.semanticscholar.org/paper/8bb6b59914487f634d3f840c66ab4dfa1882b657))\
Héctor J. García, I. Markov, and Andrew W. Cross\
*ArXiv* · Oct 24, 2012 · 28 citations

> Large-scale quantum computation is likely to require massive quantum error correction (QEC). QEC codes and circuits are described via the stabilizer formalism, which represents stabilizer states by keeping track of the operators that preserve them. Such states are obtained by stabilizer circuits (consisting of CNOT, Hadamard and Phase only) and can be represented compactly on conventional computers using Omega(n^2) bits, where n is the number of qubits. Although techniques for the efficient simulation of stabilizer circuits have been studied extensively, techniques for efficient manipulation of stabilizer states are not currently available. To this end, we design new algorithms for: (i) obtaining canonical generators for stabilizer states, (ii) obtaining canonical stabilizer circuits, and (iii) computing the inner product between stabilizer states. Our inner-product algorithm takes O(n^3) time in general, but observes quadratic behavior for many practical instances relevant to QECC (e.g., GHZ states). We prove that each n-qubit stabilizer state has exactly 4(2^n - 1) nearest-neighbor stabilizer states, and verify this claim experimentally using our algorithms. We design techniques for representing arbitrary quantum states using stabilizer frames and generalize our algorithms to compute the inner product between two such frames.
