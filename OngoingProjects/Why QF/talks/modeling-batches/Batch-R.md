```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
```
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Six photons that the suppression law does not darken

One photon enters each port of a six-mode lossless network whose transfer matrix is the discrete Fourier matrix. The suppression law of Tichy and coworkers says that an output pattern $\{t_1,\dots,t_6\}$ is dark whenever $\sum_j (j-1)\,t_j$ is not a multiple of 6. Is every pattern the law allows actually bright?

The network maps $a_k^\dagger\to\sum_j W_{jk}\,a_j^\dagger$ with $W_{jk}=\omega^{(j-1)(k-1)}/\sqrt6$ and $\omega=e^{2\pi i/6}$. With one photon per input port, a pattern $T$ has probability

$$
P(T)=\frac{\big|\operatorname{Per}W_{T}\big|^2}{\prod_j t_j!},
$$

where $W_T$ uses every column of $W$ once and repeats row $j$ $t_j$ times.

QF computes the probability of all $\binom{11}{6}$ output patterns from $|1,1,1,1,1,1\rangle$ as exact rational numbers, so a zero is exact, and sorts the patterns by $\sum_j (j-1)\,t_j \bmod 6$:

```wl
sixPhotonProbs = KeyMap[#["Name"] &, QuantumInterferometer[FourierMatrix[6]]["Probabilities", ConstantArray[1, 6]]];
sixPhotonResidue = Mod[Range[0, 5] . #, 6] &;
sixPhotonCounts = GroupBy[Keys[sixPhotonProbs], {sixPhotonResidue, sixPhotonProbs[#] === 0 &}, Length];
sixPhotonAllowedDark = Select[Keys[sixPhotonProbs], sixPhotonResidue[#] == 0 && sixPhotonProbs[#] === 0 &];
Labeled[BarChart[Table[Lookup[sixPhotonCounts[r], {False, True}, 0], {r, 0, 5}], ChartLayout -> "Stacked",
    ChartLabels -> {Range[0, 5], None}, ChartLegends -> {"bright", "exactly dark"}, Frame -> True,
    FrameLabel -> {Row[{"\[Sum] (j \[Minus] 1) ", Subscript["t", "j"], "  mod 6"}], "output patterns"}, ImageSize -> 480],
  Style[Column[{Row[{"allowed by the law, yet exactly dark (", Subscript["t", "1"], "\[Ellipsis]", Subscript["t", "6"], "):"}], Grid[Partition[Row /@ sixPhotonAllowedDark, 6], Spacings -> {2, 0.5}]}, Alignment -> Center], 13, FontFamily -> "Helvetica"], Top]
```

Every pattern whose sum is not a multiple of 6 is exactly dark, as the law requires. Most patterns with sum 0 are bright, but some are exactly dark as well (top label), each with two ports holding two photons, two holding one and two empty. The law is sufficient for darkness, not necessary.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Fully bunched photons are n! times likelier

$n$ photons enter the first $n$ ports of a six-mode lossless network with transfer matrix $W$, one per port, and all of them leave through the same port $j$. How much likelier is this for identical photons than for distinguishable ones, whatever $W$ is?

With $W_T$ the $n\times n$ matrix that takes columns $1,\dots,n$ of $W$ and repeats row $j$ of $W$ $t_j$ times for an output pattern $T$, identical (B) and distinguishable (D) photons give

$$
P_{\rm B}(T)=\frac{\big|\operatorname{Per}W_T\big|^2}{\prod_j t_j!},\qquad P_{\rm D}(T)=\frac{\operatorname{Per}\big(|(W_T)_{ab}|^2\big)}{\prod_j t_j!},
$$

and for the fully bunched pattern $W_T$ is $n$ copies of one row, so $\operatorname{Per}W_T=n!\prod_{k\le n}W_{jk}$.

The cell tests $P_{\rm B}=n!\,P_{\rm D}$ for every fully bunched pattern with all 36 entries of $W$ symbolic, $n=2,\dots,6$ (top label), and plots $P_{\rm B}/P_{\rm D}$ for every output pattern of one Haar-random network, with $P_{\rm B}$ from QF and each strip sorted by the ratio:

```wl
nfactRows[t_] := Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, t];
nfactInput[n_] := PadRight[ConstantArray[1, n], 6];
nfactBoson[w_, n_, t_] := Abs[Permanent[w[[nfactRows[t], Range[n]]]]]^2/Times @@ (t!);
nfactDist[w_, n_, t_] := Permanent[Abs[w[[nfactRows[t], Range[n]]]]^2]/Times @@ (t!);
nfactTheorem = With[{w = Array[\[FormalW], {6, 6}]},
   AllTrue[Range[2, 6], n |-> AllTrue[n IdentityMatrix[6], FullSimplify[nfactBoson[w, n, #] == n! nfactDist[w, n, #]] &]]];
nfactNetwork = BlockRandom[SeedRandom[61]; QuantumInterferometer["Random"[6]]];
nfactStrip[n_] := With[{r = SortBy[KeyValueMap[{#1["Name"], #2/nfactDist[nfactNetwork["Unitary"], n, #1["Name"]]} &,
       nfactNetwork["Probabilities", nfactInput[n]]], Last]},
   GroupBy[MapIndexed[{Max[First[#1]] == n, {n - 0.3 + 0.6 (First[#2] - 1)/(Length[r] - 1), Last[#1]}} &, r], First -> Last]];
nfactStrips = nfactStrip /@ Range[2, 6];
Labeled[ListLogPlot[{Catenate[Lookup[nfactStrips, False]], Catenate[Lookup[nfactStrips, True]],
     Catenate[Table[{{n - 0.4, n!}, {n + 0.4, n!}, {n + 0.5, None}}, {n, 2, 6}]]},
    Joined -> {False, False, True}, PlotStyle -> {Automatic, Automatic, Dashed}, PlotMarkers -> {{Automatic, 5}, {Automatic, 11}, None},
    PlotLegends -> {"other patterns", "fully bunched patterns", "n!"}, PlotRange -> All, Frame -> True,
    FrameTicks -> {{Automatic, None}, {Range[2, 6], None}}, FrameLabel -> {"photons n", Row[{Subscript["P", "B"], " / ", Subscript["P", "D"]}]}, ImageSize -> 480],
  Style[Row[{Subscript["P", "B"], " = n! ", Subscript["P", "D"], " for every fully bunched pattern, W fully symbolic, n = 2 to 6:  ", nfactTheorem}], 13, FontFamily -> "Helvetica"], Top]
```

In the sampled network every fully bunched pattern sits on $n!$, and every other pattern lies below it, spread over orders of magnitude. The symbolic test makes the factor $n!$ hold for every network, not only the sampled one.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### How often identical photons share an output port

$n$ identical photons enter $n$ of the $m$ ports of a lossless network drawn from the Haar measure on $U(m)$, one per port. Averaged over networks, how likely is it that some port receives two or more photons, as $m$ grows with $n$?

The Haar average of $U|S\rangle\langle S|U^\dagger$ commutes with every network, and the $n$-photon states, $\binom{n+m-1}{n}$ of them, carry an irreducible representation of $U(m)$, so by Schur's lemma the average is the identity over that dimension: every output pattern is equally likely. Only $\binom mn$ patterns have at most one photon per port, so

$$
\overline{P}_{\rm bunch}(n,m)=1-\binom mn\Big/\binom{n+m-1}{n}.
$$

The cell evaluates this law for $m=n$, $m=n^2$ and $m=4n^2$ (lines), takes its limit $n\to\infty$ in closed form for $m=n$ and for $m=c\,n^2$ (dashed, top label), and checks it against the average of QF's exact output probabilities over 400 Haar-random networks (markers):

```wl
haarBunchLaw[n_, m_] := 1 - Binomial[m, n]/Binomial[n + m - 1, n];
haarBunchLimits = {Limit[haarBunchLaw[n, n], n -> Infinity], Limit[haarBunchLaw[n, c n^2], n -> Infinity, Assumptions -> c > 0]};
haarBunchModes = {# &, #^2 &, 4 #^2 &};
haarBunchSample[n_, m_] := BlockRandom[SeedRandom[100 n + m];
   Mean @ Table[Total @ KeySelect[QuantumInterferometer["Random"[m]]["Probabilities", PadRight[ConstantArray[1, n], m]], Max[#["Name"]] >= 2 &], 400]];
Labeled[Show[
   ListLinePlot[Table[{n, haarBunchLaw[n, f[n]]}, {f, haarBunchModes}, {n, 2, 12}], PlotRange -> All,
    PlotLegends -> {"m = n", Row[{"m = ", Superscript["n", 2]}], Row[{"m = 4", Superscript["n", 2]}]}],
   ListLinePlot[Table[{{2, l}, {12, l}}, {l, {First[haarBunchLimits], Last[haarBunchLimits] /. c -> 1, Last[haarBunchLimits] /. c -> 4}}], PlotStyle -> Dashed, PlotRange -> All],
   ListPlot[MapThread[Table[{n, haarBunchSample[n, #1[n]]}, {n, 2, #2}] &, {haarBunchModes, {4, 3, 2}}], PlotMarkers -> {Automatic, 11}, PlotRange -> All],
   PlotRange -> {{1.6, 12.4}, {0, 1.02}}, Frame -> True, FrameLabel -> {"photons n", "Haar-averaged probability of a shared port"}, ImageSize -> 480],
  Style[Row[{"n \[Rule] \[Infinity]:   m = n gives ", First[haarBunchLimits], ",   m = c ", Superscript["n", 2], " gives ", TraditionalForm[Last[haarBunchLimits]]}], 13, FontFamily -> "Helvetica"], Top]
```

With as many ports as photons a shared port becomes certain as $n$ grows. With $m=c\,n^2$ the probability settles at $1-e^{-1/c}$, which is small for large $c$, so with $m$ a large multiple of $n^2$ most detection events stay free of collisions at every $n$. The Haar averages of QF's probabilities fall on the law.
#| style: SlideShowNavigationBar
#| tags: SlideShowHeader

&nbsp;

### Twin Fock states on a balanced beam splitter

$n$ photons enter each port of a balanced beam splitter. What is the probability of finding $t_1$ photons in the first output, for every $n$?

The splitter maps $a_1^\dagger\to(a_1^\dagger+a_2^\dagger)/\sqrt2$ and $a_2^\dagger\to(a_2^\dagger-a_1^\dagger)/\sqrt2$, so each input pair becomes $a_1^\dagger a_2^\dagger\to(a_2^{\dagger2}-a_1^{\dagger2})/2$ and

$$
|n,n\rangle=\frac{(a_1^\dagger a_2^\dagger)^n}{n!}|0,0\rangle\;\to\;\frac{(a_2^{\dagger2}-a_1^{\dagger2})^n}{2^n\,n!}|0,0\rangle .
$$

Only even $t_1=2k$ occur, and the binomial expansion gives

$$
P(2k)=\binom{2k}{k}\binom{2n-2k}{n-k}\Big/4^n,\qquad P(2k+1)=0 .
$$

The cell sums the law over $k$ for symbolic $n$ (top label) and compares it with QF's exact probabilities of every output $t_1=0,\dots,2n$:

```wl
twinFockLaw[n_, k_] := Binomial[2 k, k] Binomial[2 (n - k), n - k]/4^n;
twinFockNorm = Sum[twinFockLaw[n, k], {k, 0, n}, Assumptions -> Element[n, Integers] && n >= 0];
twinFockSplitter = QuantumInterferometer[{"BS"[\[Pi]/4]}];
twinFockPanel[n_] := ListPlot[{KeyValueMap[{First[#1["Name"]], #2} &, twinFockSplitter["Probabilities", {n, n}]], Table[{2 k, twinFockLaw[n, k]}, {k, 0, n}]},
   Joined -> {False, True}, Filling -> {1 -> Axis}, PlotStyle -> {Automatic, Dashed}, PlotMarkers -> {{Automatic, 9}, None}, Axes -> False,
   PlotLegends -> {"QF, every output", Row[{"law at ", Subscript["t", "1"], " = 2k"}]}, PlotRange -> {{-0.5, 2 n + 0.5}, {-0.012, All}}, Frame -> True,
   FrameLabel -> {Subscript["t", "1"], Row[{"P(", Subscript["t", "1"], ")"}]}, PlotLabel -> Row[{"n = ", n}], ImageSize -> 300];
twinFockPanels = twinFockPanel /@ {2, 5, 10};
Labeled[Legended[Row[First /@ twinFockPanels], Last[First[twinFockPanels]]],
  Style[Row[{"\[Sum]", Subscript["", "k"], " P(2k) for symbolic n:  ", twinFockNorm}], 13, FontFamily -> "Helvetica"], Top]
```

Every odd output is exactly zero, and QF's even outputs lie on the law for each $n$. The distribution is largest at the edges and smallest near the middle, so the likeliest outcome sends all $2n$ photons out of one port, and the closed form is normalized for every $n$.
