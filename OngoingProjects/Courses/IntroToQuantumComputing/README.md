Cloud Object TOC: [link](https://www.wolframcloud.com/obj/mohammadb/Published/Introduction-to-quantum-computing.nb)

# Introduction to Quantum Computing


![image](https://www.wolframcloud.com/obj/dd249149-18cb-4eb2-9688-afed1808ade7)


Follow the links for lesson notebooks

## Table of Contents

1. [What Is Quantum Computation?](#1-what-is-quantum-computation)
   - [1.1 Classical Computation](#11-classical-computation)
   - [1.2 Quantum Circuits and Qubits](#12-quantum-circuits-and-qubits)

2. [Quantum States: Single Particle Case](#2-quantum-states-single-particle-case)
   - [2.1 Bloch Sphere](#21-bloch-sphere)
   - [2.2 Bra-Ket Notation](#22-bra-ket-notation)
   - [2.3 State Vectors and Amplitudes](#23-state-vectors-and-amplitudes)
   - [2.4 Pauli Operators and Bloch Sphere](#24-pauli-operators-and-bloch-sphere)
   - [2.5 Higher Dimensions](#25-higher-dimensions)

3. [Probability and Measurement](#3-probability-and-measurement)
   - [3.1 Frequency vs Probability](#31-frequency-vs-probability)
   - [3.2 Joint vs Marginal Probabilities](#32-joint-vs-marginal-probabilities)

4. [Quantum States: Multi-Particle Cases](#4-quantum-states-multi-particle-cases)
   - [4.1 Tensor Product Structure](#41-tensor-product-structure)
   - [4.2 Separable vs Nonseparable Objects](#42-separable-vs-nonseparable-objects)

5. [Superposition and Entanglement](#5-superposition-and-entanglement)
   - [5.1 Superposition](#51-superposition)
   - [5.2 Entanglement](#52-entanglement)

6. [A Quick Tour of Quantum Gates](#6-a-quick-tour-of-quantum-gates)
   - [6.1 Single Qubit Gates](#61-single-qubit-gates)
   - [6.2 Multi-Qubit Gates](#62-multi-qubit-gates)

7. [Deutsch-Jozsa Algorithm](#7-deutsch-jozsa-algorithm)
   - [7.1 Statement of the Deutsch-Jozsa Problem](#71-statement-of-the-deutsch-jozsa-problem)
   - [7.2 Deutsch-Jozsa Oracles](#72-deutsch-jozsa-oracles)
   - [7.3 Deutsch-Jozsa Algorithm](#73-deutsch-jozsa-algorithm)

8. [Bernstein-Vazirani Algorithm](#8-bernstein-vazirani-algorithm)
   - [8.1 Statement of the Bernstein-Vazirani Problem](#81-statement-of-the-bernstein-vazirani-problem)
   - [8.2 Quantum Bernstein-Vazirani Algorithm](#82-quantum-bernstein-vazirani-algorithm)
   - [8.3 Effects of Noise](#83-effects-of-noise)

9. [Grover Search Algorithm](#9-grover-search-algorithm)
   - [9.1 Statement of the Grover Search Problem](#91-statement-of-the-grover-search-problem)
   - [9.2 Grover Oracle](#92-grover-oracle)
   - [9.3 Grover's Algorithm](#93-grovers-algorithm)

10. [Quantum Arithmetic](#10-quantum-arithmetic)
    - [10.1 Discrete Fourier Transform and QFT](#101-discrete-fourier-transform-and-qft)
    - [10.2 Phase Basis Encoding](#102-phase-basis-encoding)
    - [10.3 Quantum Phase Adders](#103-quantum-phase-adders)
    - [10.4 Generic Addition](#104-generic-addition)
    - [10.5 Multiplication](#105-multiplication)

11. [Quantum Phase Estimation Algorithm](#11-quantum-phase-estimation-algorithm)
    - [11.1 Eigenvectors and Eigenvalues](#111-eigenvectors-and-eigenvalues)
    - [11.2 Quantum Phase Estimation: Simple Example](#112-quantum-phase-estimation-simple-example)
    - [11.3 Quantum Phase Estimation: Complicated Example](#113-quantum-phase-estimation-complicated-example)

12. [Shor's Algorithm for Factorization](#12-shors-algorithm-for-factorization)
    - [12.1 Factorization Through Order Finding](#121-factorization-through-order-finding)
    - [12.2 Order Finding with a Quantum Circuit](#122-order-finding-with-a-quantum-circuit)
