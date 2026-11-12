```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Five-qubit code: fifteen single-qubit errors, fifteen distinct syndromes

Five qubits store one logical qubit in the joint $+1$ eigenspace of four commuting checks. One qubit suffers an $X$, $Y$ or $Z$ error, and the four checks are measured. Do the four outcomes tell all fifteen errors apart, for either logical basis state?

The checks are the cyclic shifts of $XZZXI$, and the logical operators are

$$
g_1=XZZXI,\quad g_2=IXZZX,\quad g_3=XIXZZ,\quad g_4=ZXIXZ,\qquad \bar Z=ZZZZZ,\quad \bar X=XXXXX .
$$

After an error $E$ on a codeword $|\psi\rangle$, check $g_k$ reads the bit $s_k=\left(1-\langle\psi|E^\dagger g_kE|\psi\rangle\right)/2$, which is 1 exactly when $E$ anticommutes with $g_k$.

`PauliStabilizer` stores $|0_L\rangle$ as the tableau of $g_1,\dots,g_4,\bar Z$; QF applies each error to $|0_L\rangle$ and to $|1_L\rangle=\bar X|0_L\rangle$ and reads the four bits (black is 1):

```wl
fiveCode = PauliStabilizer["5QubitCode"]; fiveChecks = Most[fiveCode["Stabilizers"]];
fiveErrors = Flatten[Outer[Rule, {"X", "Y", "Z"}, Range[5]]];
fivePanel[label_, ps_] := With[{s = Table[(1 - ps[{e}]["Expectation", #] & /@ fiveChecks)/2, {e, fiveErrors}]},
   Labeled[ArrayPlot[Transpose[s], Mesh -> All, Frame -> True, ImageSize -> 520,
     FrameTicks -> {{Transpose[{Range[4], fiveChecks}], None}, {Transpose[{Range[15], First[#] <> ToString[Last[#]] & /@ fiveErrors}], None}}],
    Row[{label, ":  ", Length[DeleteCases[Union[s], {0, 0, 0, 0}]], " distinct nonzero syndromes"}], Top]];
Column[{fivePanel[Row[{"|", Subscript[0, "L"], "\[RightAngleBracket]"}], fiveCode], fivePanel[Row[{"|", Subscript[1, "L"], "\[RightAngleBracket]"}], fiveCode[Thread["X" -> Range[5]]]]}, Spacings -> 2]
```

The two tables are identical, as they must be: an error flips exactly the checks it anticommutes with, independently of the codeword. Every column is a different nonzero four-bit word, and four bits have only fifteen nonzero values, so the checks have no room to spare; this is why the five-qubit code is called perfect, and why a decoder can invert the table.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### A thousand-qubit GHZ state: entanglement entropy as a rank over F₂

Take $n=1000$ qubits in the GHZ state $(|0\cdots0\rangle+|1\cdots1\rangle)/\sqrt2$, and for comparison in the one-dimensional cluster state, made by putting every qubit in $|+\rangle$ and applying $CZ$ to each neighboring pair. How much entanglement separates a set $A$ of $k$ qubits from the rest, when $A$ is a block at the end of the chain and when $A$ is scattered at random?

A stabilizer state is fixed by $n$ commuting Pauli generators: $X^{\otimes n}$ and $Z_iZ_{i+1}$ for GHZ, $Z_{i-1}X_iZ_{i+1}$ for the cluster. Write each generator as a binary vector $(x|z)\in\mathbb F_2^{2n}$. Then the entanglement entropy, in bits, is (Fattal et al. 2004)

$$
S(A)=\operatorname{rank}_{\mathbb F_2}G_A-|A|,
$$

Stack the generators as the rows of an $n\times2n$ matrix $G$; $G_A$ keeps only the $x$ and $z$ columns of the qubits in $A$.

`PauliStabilizer` builds both states as tableaux and evaluates this rank for $k=50,100,\dots,950$, with one random set $A$ per $k$; no $2^{1000}$-entry vector is ever formed:

```wl
ghzChain = PauliStabilizer[1000][Prepend[Table["CNOT" -> {q, q + 1}, {q, 999}], "H"]];
clusterChain = PauliStabilizer[1000][Join[Table["H" -> q, {q, 1000}], Table["CZ" -> {q, q + 1}, {q, 999}]]];
entropyPanel[label_, ps_] := ListLogPlot[{Table[{k, ps["Entropy", Range[k]]}, {k, 50, 950, 50}],
    Table[{k, ps["Entropy", RandomSample[Range[1000], k]]}, {k, 50, 950, 50}]}, Joined -> True, PlotMarkers -> "OpenMarkers", PlotRange -> {0.7, 500},
   PlotLegends -> {"block {1, \[Ellipsis], k}", "k qubits at random"}, Frame -> True, FrameLabel -> {"k", "S(A)  (bits)"}, PlotLabel -> label, ImageSize -> 340];
SeedRandom[1]; entropyPanels = {entropyPanel["GHZ", ghzChain], entropyPanel["cluster", clusterChain]};
Legended[Row[First /@ entropyPanels, "  "], Last[First[entropyPanels]]]
```

The GHZ state carries exactly one bit across every cut, a block or a scattered set of any size: its entanglement is a single shared parity. The cluster state also carries one bit across a block boundary, but far more when the same number of qubits is scattered, so its entanglement depends on the shape of the cut, not only on its size.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Stabilizer measurement kept symbolic: random outcomes as formal bits

Two qubits start in $|+\rangle$, and two CNOTs copy their parity onto a third qubit. Then an $X$ on qubit 1, an $X$ on qubit 3, or a $Z$ on qubit 1 may strike, and all three qubits are measured in the $Z$ basis. Which outcomes are random, and how is the third fixed by the first two?

Before the error the state and its stabilizer generators are

$$
\mathrm{CNOT}_{23}\,\mathrm{CNOT}_{13}\,H_1H_2|000\rangle=\tfrac12\sum_{m_1,m_2}|m_1,m_2,m_1\veebar m_2\rangle,\qquad \langle X_1X_3,\;X_2X_3,\;Z_1Z_2Z_3\rangle .
$$

A measurement whose outcome is random does not branch: it allocates a formal bit $s_k$, which enters its generator as the sign $(-1)^{s_k}=1-2s_k$. A measurement whose outcome is already determined is stored as a Boolean function of the earlier bits.

QF measures each state once with formal bits; the title shows how many formal bits were allocated and the stored third outcome $m_3$, and the bars count 1000 shots drawn by substituting random values of $s_1,s_2$:

```wl
parityState = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
parityPanel[label_, ps_] := With[{ms = ps["SymbolicMeasure", {1, 2, 3}]},
   With[{s = Union[Cases[ms["Phase"], _\[FormalS], Infinity]], m3 = First[Values[ms["Outcomes"]]]},
    BarChart[Lookup[Counts[First[Keys[#["M", {1, 2, 3}]]] & /@ ms["SampleOutcomes", 1000]], Tuples[{0, 1}, 3], 0],
     ChartLabels -> (StringJoin[ToString /@ #] & /@ Tuples[{0, 1}, 3]), Frame -> True, FrameLabel -> {Row[Subscript["m", #] & /@ Range[3]], "shots"},
     PlotLabel -> Row[{label, ":  ", Length[s], " formal bits,  ", Subscript["m", 3], " = ", TraditionalForm[BooleanFunction[
          Thread[Tuples[{0, 1}, 2] -> (m3 /. Thread[s -> #] & /@ Tuples[{0, 1}, 2])], {Subscript["s", 1], Subscript["s", 2]}, "ANF"]]}], ImageSize -> 260]]];
SeedRandom[1];
Grid[{{parityPanel["no error", parityState], parityPanel["X on qubit 1", parityState[{"X" -> 1}]]},
  {parityPanel["X on qubit 3", parityState[{"X" -> 3}]], parityPanel["Z on qubit 1", parityState[{"Z" -> 1}]]}}]
```

Each state allocates two formal bits, and the shots split evenly over four strings: two fair coins and one dependent bit. Without an error $m_3=s_1\veebar s_2$ and every shot has even parity. An $X$ on qubit 1 or on qubit 3 negates the stored function and moves every shot to odd parity, which is how a parity check detects a bit flip; a $Z$ error changes nothing that a $Z$ measurement sees.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### T gates as a sum of stabilizer states: how many terms survive compression

The gate $T=\mathrm{diag}(1,e^{i\pi/4})$ is not a Clifford gate, so a state made with $t$ T gates is stored as a sum of stabilizer states. On $n$ qubits, two circuits use $t$ T gates: in the first, qubits $1,\dots,t$ each get $H$ then $T$, and a CNOT ladder then entangles the chain; in the second, the ladder first makes a GHZ state and the T gates follow. How many terms does the frame keep after compression as $n$ grows?

As an operator, $T$ on qubit $q$ is a sum of two Pauli operators,

$$
T_q=\frac{1+e^{i\pi/4}}{2}\,\mathbb 1+\frac{1-e^{i\pi/4}}{2}\,Z_q,\qquad |\psi\rangle=\sum_{i=1}^{2^t}c_i\,|s_i\rangle ,
$$

so each T gate doubles the list of stabilizer states $|s_i\rangle$ in the frame. The overlaps $\langle s_i|s_j\rangle$ follow from the tableaux, and from them `"Compress"` keeps a largest linearly independent subset of the $|s_i\rangle$ and recomputes the $c_i$ exactly, so $|\psi\rangle$ is unchanged.

QF builds each state as a `StabilizerFrame` and counts its terms after `"Compress"`, against the $2^n$ amplitudes of the state vector (dashed):

```wl
magicLadder[n_] := Table["CNOT" -> {i, i + 1}, {i, n - 1}];
magicTerms[n_, t_, tFirst_] := PauliStabilizer[n][If[tFirst, Join[Thread["H" -> Range[t]], Thread["T" -> Range[t]], magicLadder[n]],
      Join[{"H" -> 1}, magicLadder[n], Thread["T" -> Range[t]]]]]["Compress"]["Length"];
magicPanel[label_, tFirst_] := ListLogPlot[Append[Table[{n, magicTerms[n, t, tFirst]}, {t, 5}, {n, 5, 30, 5}], Table[{n, 2^n}, {n, 5, 30, 5}]],
   Joined -> True, PlotMarkers -> "OpenMarkers", PlotStyle -> {Automatic, Automatic, Automatic, Automatic, Automatic, Dashed}, AspectRatio -> 1, PlotRange -> {{4, 31}, All},
   PlotLegends -> Append[Row[{"t = ", #}] & /@ Range[5], Row[{Superscript[2, "n"], " amplitudes"}]],
   Frame -> True, FrameLabel -> {"n", "stabilizer terms"}, PlotLabel -> label, ImageSize -> 320];
magicPanels = {magicPanel["H, T on qubits 1..t, then ladder", True], magicPanel["ladder, then T on qubits 1..t", False]};
Legended[Row[First /@ magicPanels, "  "], Last[First[magicPanels]]]
```

When each T gate acts on its own $|+\rangle$, no term is redundant: the frame keeps $2^t$ terms at every $n$, while the state vector grows as $2^n$. When the same T gates act on the GHZ state, they only change the relative phase of $|0\cdots0\rangle$ and $|1\cdots1\rangle$, and two terms suffice for every $t$. The count is the dimension spanned by the frame, an upper bound on the stabilizer rank, not the rank itself.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Repetition code: the logical error rate is a polynomial in the flip probability

One bit is stored in an odd number $n$ of qubits as $|0\cdots0\rangle$. Each qubit is flipped independently with probability $p$, and the bit is decoded by majority vote. What is the probability $P_L$ that the decoded bit is wrong, and for which $p$ does a longer code help?

Each qubit passes through the bit-flip channel

$$
\mathcal E_p(\rho)=(1-p)\,\rho+p\,X\rho X,
$$

and majority vote fails when more than $n/2$ qubits flipped.

QF applies `QuantumChannel["BitFlip"[p]]` to each qubit of the stabilizer state $|0\cdots0\rangle$; every application splits each branch into an unflipped and a flipped tableau with weights $1-p$ and $p$, and qubit $k$ reads 1 when its generator $Z_k$ has sign $-1$. Summing the weights of the failing branches gives $P_L(p)$ exactly:

```wl
repFlip[branches_, q_] := Join @@ (Function[{w, ps}, {w #1, #2} & @@@ QuantumChannel["BitFlip"[p], {q}][ps]] @@@ branches);
repRate[n_] := Expand[Total[Cases[Fold[repFlip, {{1, PauliStabilizer[n]}}, Range[n]], {w_, ps_} /; Count[ps["StabilizerSigns"], -1] > n/2 :> w]]];
repRates = repRate /@ {1, 3, 5, 7};
repPanels = {Plot[Evaluate[repRates], {p, 0, 1}, PlotLegends -> (Row[{"n = ", #}] & /@ {1, 3, 5, 7}), Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]},
    PlotLabel -> Row[{Subscript["P", "L"], "(5) < ", Subscript["P", "L"], "(3)  \[DoubleLongLeftRightArrow]  ", Reduce[repRates[[3]] < repRates[[2]] && 0 < p < 1, p]}], ImageSize -> 320],
   LogLogPlot[Evaluate[repRates], {p, 10^-3, 1/2}, Frame -> True, FrameLabel -> {"p", Subscript["P", "L"]},
    PlotLabel -> Row[{Subscript["P", "L"], "(3) = ", repRates[[2]]}], ImageSize -> 320]};
Legended[Row[{First[First[repPanels]], Last[repPanels]}, "  "], Last[First[repPanels]]]
```

For $n=3$ the rate has no linear term, so at small $p$ the code fails far less often than a bare qubit; on the log-log axes each code becomes a straight line at small $p$, with slope equal to the smallest number of flips that defeats the majority. All curves meet at $p=1/2$, and `Reduce` makes the threshold exact: below one half the longer code is better, above it worse.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Variational ground state: plain against natural gradient descent

A variational quantum eigensolver looks for the ground state of $H=X$ by adjusting the two angles of a one-qubit circuit to lower $\langle H\rangle$. Plain gradient descent steps against the gradient; quantum natural gradient descent first rescales it by the geometry of the state space. Which reaches the ground energy, and how fast, from different starting angles and learning rates $\eta$?

The circuit, its energy and the update rules are

$$
|\psi(\theta,\phi)\rangle=P(2\phi)\,R_Y(2\theta)|0\rangle=\cos\theta\,|0\rangle+e^{2i\phi}\sin\theta\,|1\rangle,\qquad E(\theta,\phi)=\langle\psi|X|\psi\rangle=\sin2\theta\cos2\phi,
$$

$$
\boldsymbol\theta_{k+1}=\boldsymbol\theta_k-\eta\,\nabla E\quad\text{(plain)},\qquad \boldsymbol\theta_{k+1}=\boldsymbol\theta_k-\eta\,g^{-1}\nabla E\quad\text{(natural)},
$$

with the Fubini-Study metric $g_{ij}=\mathrm{Re}\left(\langle\partial_i\psi|\partial_j\psi\rangle-\langle\partial_i\psi|\psi\rangle\langle\psi|\partial_j\psi\rangle\right)=\mathrm{diag}(1,\sin^22\theta)$. It vanishes in the $\phi$ direction at the poles $\theta=0,\pi/2$, where $\phi$ no longer changes the state.

QF builds the state from the circuit and `FubiniStudyMetricTensor` gives $g$; both optimizers run 40 steps:

```wl
Needs["Wolfram`QuantumFramework`QuantumOptimization`"];
vqeState = QuantumCircuitOperator[{"RY"[2 \[FormalTheta]], "P"[2 \[FormalPhi]]}][];
vqeEnergy[th_, ph_] = ComplexExpand[vqeState["Dagger"][QuantumOperator["X"][vqeState]]["Scalar"]] /. {\[FormalTheta] -> th, \[FormalPhi] -> ph};
vqeMetric = FubiniStudyMetricTensor[vqeState, "Parameters" -> {\[FormalTheta], \[FormalPhi]}];
vqePanel[\[Eta]_, start_] := Show[ListLinePlot[{vqeEnergy @@@ Re[QuantumNaturalGradientDescent[vqeEnergy, vqeMetric, "InitialPoint" -> start, "LearningRate" -> \[Eta], "MaxIterations" -> 40]],
      vqeEnergy @@@ GradientDescent[vqeEnergy, start, "LearningRate" -> \[Eta], "MaxIterations" -> 40]}, DataRange -> {0, 40}, PlotRange -> All,
     PlotLegends -> {"natural gradient", "plain gradient"}],
    PlotRange -> {-1.1, 1.1}, Axes -> False, Frame -> True, FrameLabel -> {"step", "E"}, PlotLabel -> Row[{"\[Eta] = ", \[Eta], ",   (\[Theta], \[Phi]) = (", Row[start, ", "], ")"}], ImageSize -> 280];
vqePanels = Outer[vqePanel, {0.05, 0.2}, {{\[Pi]/12, \[Pi]/12}, {\[Pi]/3, \[Pi]/12}, {0.01, 1}}, 1];
Legended[Grid[Map[First, vqePanels, {2}]], Last[vqePanels[[1, 1]]]]
```

At the smaller rate the natural gradient reaches the ground energy in fewer steps from every start. At the larger rate its rescaled steps first push the energy up from the other two starts, and plain descent gets there sooner. From $(\pi/3,\pi/12)$ plain descent, at either rate, stalls at $E=0$ at a pole, where the state no longer depends on $\phi$; the natural gradient divides the vanishing $\phi$ gradient by the vanishing metric and escapes.

#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Quantum linear solver: HHL on four qubits

Solve $A\boldsymbol x=\boldsymbol b$ for $A=\begin{pmatrix}1&-1/3\\-1/3&1\end{pmatrix}$ and the unit vector $\boldsymbol b=(\cos\beta,\sin\beta)$. The Harrow-Hassidim-Lloyd (HHL) circuit loads $\boldsymbol b$ into one qubit, writes the eigenvalues of $A$ into a two-qubit clock by phase estimation, rotates an ancilla by an angle set by each eigenvalue, clears the clock, and keeps only the runs where the ancilla reads 1. Does the kept qubit give the probabilities $|x_j|^2/|\boldsymbol x|^2$ of the solution, and how often is a run kept?

$A$ has eigenvalues $\lambda_j=2/3,\,4/3$ with eigenvectors $|u_j\rangle$, and $|b\rangle=\sum_jb_j|u_j\rangle$. Phase estimation of $U=e^{iAt}$ with $t=3\pi/4$ writes $\tilde\lambda_j=2^2\lambda_jt/2\pi=1,2$ into the clock, exactly, because the eigenvalue ratio is 2. The controlled rotation $R_Y(\theta_j)$ with $\sin(\theta_j/2)=1/\tilde\lambda_j$, followed by clearing the clock, gives

$$
\sum_jb_j|u_j\rangle|0\rangle_a\;\longrightarrow\;\sum_jb_j|u_j\rangle\left(\sqrt{1-\tilde\lambda_j^{-2}}\,|0\rangle_a+\tilde\lambda_j^{-1}|1\rangle_a\right),
$$

and the branch with the ancilla in $|1\rangle$ is proportional to $A^{-1}|b\rangle$.

Each clock value sets one bit, so each rotation needs a single control. QF assembles the circuit from the named `"PhaseEstimation"` circuit (without its initial $X$ and final measurements), measures ancilla and solution qubit, and compares with `LinearSolve` (line):

```wl
hhlA = {{1, -1/3}, {-1/3, 1}};
hhlQPE = QuantumCircuitOperator["PhaseEstimation"[QuantumOperator[MatrixExp[I 3 \[Pi]/4 hhlA]], 2]][[2 ;; -3]];
hhlProbs[\[Beta]_] := QuantumCircuitOperator[{QuantumState[{Cos[\[Beta]], Sin[\[Beta]]}] -> {4}, hhlQPE -> {2, 3, 4},
     Splice[MapThread["C"["RY"[#1]] -> {#2, 1} &, {2 ArcSin[1/{1, 2}], {3, 2}}]], hhlQPE["Dagger"] -> {2, 3, 4}, {1}, {4}}][]["ProbabilitiesList"];
hhlPoints = Table[With[{pr = hhlProbs[\[Beta]]}, {{\[Beta], pr[[4]]/(pr[[3]] + pr[[4]])}, {\[Beta], pr[[3]] + pr[[4]]}}], {\[Beta], 0, \[Pi], \[Pi]/8}];
ListPlot[{Table[{\[Beta], Last[Normalize[LinearSolve[hhlA, {Cos[\[Beta]], Sin[\[Beta]]}]]]^2}, {\[Beta], 0, \[Pi], \[Pi]/100}], hhlPoints[[All, 1]], hhlPoints[[All, 2]]},
 Joined -> {True, False, False}, PlotMarkers -> {None, Automatic, Automatic}, Frame -> True, FrameLabel -> {"\[Beta]", "probability"},
 FrameTicks -> {{Automatic, None}, {Range[0, \[Pi], \[Pi]/4], None}}, PlotRange -> {0, 1.05},
 PlotLegends -> {Row[{"LinearSolve:  ", Subsuperscript["x", 2, 2], "/|x|\.b2"}], "HHL:  P(solution qubit 1 | ancilla 1)", "HHL:  P(ancilla 1)"}, ImageSize -> 480]
```

At every $\beta$ the HHL points sit on the `LinearSolve` curve: the kept qubit carries $|x_2|^2/|\boldsymbol x|^2$. The price is the post-selection. A run is always kept when $\boldsymbol b$ lies along the eigenvector of the smaller eigenvalue, whose rotation turns the ancilla fully to $|1\rangle$, and only a quarter of the time along the other.
