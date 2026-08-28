(* backbone-entanglement-sudden-death-checks.wl *)
(* The exact-decision battery behind the tails of backbone-entanglement-sudden-death.md. *)
(* Run: wolframscript -file backbone-entanglement-sudden-death-checks.wl                 *)
(* Every line prints a label and a verdict; the page's tails compress these calls.       *)
(* Conventions: |0> = {1,0} is the ground level; basis order |00>,|01>,|10>,|11>;        *)
(* sigma_minus = |0><1|; X-state populations a,b,c,d, outer coherence w, inner z,        *)
(* both real and nonnegative by local phase rotations; damped fraction p = 1 - E^(-g t). *)

ws[label_, val_] := WriteString[$Output, label, ": ", ToString[val, InputForm], "\n"];

id2 = IdentityMatrix[2];
sm = {{0, 1}, {0, 0}};
sz = {{1, 0}, {0, -1}};
sx = {{0, 1}, {1, 0}};
zero4 = ConstantArray[0, {4, 4}];
dis[op_][rh_] := op.rh.ConjugateTranspose[op] - (1/2) (ConjugateTranspose[op].op.rh + rh.ConjugateTranspose[op].op);
smA = KroneckerProduct[sm, id2]; smB = KroneckerProduct[id2, sm];
szA = KroneckerProduct[sz, id2]; szB = KroneckerProduct[id2, sz];
xstate = {{a, 0, 0, w}, {0, b, z, 0}, {0, z, c, 0}, {w, 0, 0, d}};
statecons = a >= 0 && b >= 0 && c >= 0 && d >= 0 && w >= 0 && z >= 0 && w^2 <= a*d && z^2 <= b*c;

(* ---- A. the local amplitude-damping channel (S1) ---- *)
k0 = {{1, 0}, {0, Sqrt[1 - p]}}; k1 = {{0, Sqrt[p]}, {0, 0}};
kraus = Flatten[Table[KroneckerProduct[ki, kj], {ki, {k0, k1}}, {kj, {k0, k1}}], 1];
chan[rh_] := Total[(#.rh.Transpose[#]) & /@ kraus];
mapClaim = {{a + p*b + p*c + p^2*d, 0, 0, (1 - p)*w}, {0, (1 - p)*(b + p*d), (1 - p)*z, 0}, {0, (1 - p)*z, (1 - p)*(c + p*d), 0}, {(1 - p)*w, 0, 0, (1 - p)^2*d}};
ws["A1 six-parameter map on the X class", Simplify[chan[xstate] - mapClaim] === zero4];
rhot = mapClaim /. p -> 1 - Exp[-g*t];
ws["A2 the map solves the local Lindblad equation", FullSimplify[D[rhot, t] - g*(dis[smA][rhot] + dis[smB][rhot]), Assumptions -> g > 0 && t >= 0] === zero4];
chan1[rh_] := k0.rh.Transpose[k0] + k1.rh.Transpose[k1];
oneq = chan1[{{al^2, al*be}, {al*be, be^2}}];
ws["A3 lone-qubit coherence carries Sqrt[1-p]", Simplify[oneq[[1, 2]] - al*be*Sqrt[1 - p], Assumptions -> 0 <= p < 1] === 0];
rphi = mapClaim /. {a -> al^2, b -> 0, c -> 0, d -> be^2, w -> al*be, z -> 0};
rpsi = mapClaim /. {a -> 0, b -> al^2, c -> be^2, d -> 0, w -> 0, z -> al*be};
redA = {{rphi[[1, 1]] + rphi[[2, 2]], rphi[[1, 3]] + rphi[[2, 4]]}, {rphi[[3, 1]] + rphi[[4, 2]], rphi[[3, 3]] + rphi[[4, 4]]}};
ws["A4 reduced state of either qubit on evolved Phi", Simplify[redA - {{al^2 + be^2*p, 0}, {0, be^2*(1 - p)}}] === {{0, 0}, {0, 0}}];
ws["A5 excited population of qubit A is (1-p)(c+d) on the whole X class", Simplify[(mapClaim[[3, 3]] + mapClaim[[4, 4]]) - (1 - p)*(c + d)] === 0];
ws["A6 excited population of qubit B is (1-p)(b+d) on the whole X class", Simplify[(mapClaim[[2, 2]] + mapClaim[[4, 4]]) - (1 - p)*(b + d)] === 0];

(* ---- B. Wootters concurrence on X states (S1) ---- *)
sy = {{0, -I}, {I, 0}};
yy = KroneckerProduct[sy, sy];
rtil = yy.xstate.yy;
evClaim = {(Sqrt[a*d] + w)^2, (Sqrt[a*d] - w)^2, (Sqrt[b*c] + z)^2, (Sqrt[b*c] - z)^2};
ws["B1 spectrum of rho.rhotilde on X states", Simplify[Expand[CharacteristicPolynomial[xstate.rtil, x] - Times @@ ((x - #) & /@ evClaim)]] === 0];
ws["B2 max form of the concurrence", Resolve[ForAll[{sa, sb, wv, zv}, sa >= wv >= 0 && sb >= zv >= 0, Max[sa + wv, sa - wv, sb + zv, sb - zv] - sa - sb == Max[wv - sb, zv - sa]], Reals]];
wootters[rh_] := Block[{rt, lam}, rt = yy.Conjugate[rh].yy; lam = Sort[Sqrt[Eigenvalues[rh.rt]], Greater]; Max[0, lam[[1]] - lam[[2]] - lam[[3]] - lam[[4]]]];
(* wootters sorts with Greater: safe on exact NUMERIC instances only, never on symbolic entries *)
xnum[av_, bv_, cv_, dv_, wv_, zv_] := {{av, 0, 0, wv}, {0, bv, zv, 0}, {0, zv, cv, 0}, {wv, 0, 0, dv}};
spot[av_, bv_, cv_, dv_, wv_, zv_] := RootReduce[wootters[xnum[av, bv, cv, dv, wv, zv]] - 2*Max[0, wv - Sqrt[bv*cv], zv - Sqrt[av*dv]]] === 0;
ws["B3 exact instance 1 (outer-entangled)", spot[1/3, 1/6, 1/6, 1/3, 1/4, 1/8]];
ws["B4 exact instance 2 (outer-entangled)", spot[1/2, 1/8, 1/8, 1/4, 1/3, 1/9]];
ws["B5 exact instance 3 (inner-entangled)", spot[1/8, 3/8, 3/8, 1/8, 1/10, 1/3]];
ws["B6 at most one live branch", Resolve[ForAll[{a, b, c, d, w, z}, statecons, ! (w^2 > b*c && z^2 > a*d)], Reals]];

(* ---- C. death and survival under loss (S2) ---- *)
cons = al^2 + be^2 == 1 && 0 < al < 1 && 0 < be < 1;
lamPhi = (1 - p)*be*(al - be*p);
ws["C1 branch identity on evolved Phi", Simplify[rphi[[1, 4]] - Sqrt[rphi[[2, 2]]*rphi[[3, 3]]] - lamPhi, Assumptions -> 0 <= p < 1 && 0 < al < 1 && 0 < be < 1] === 0];
ws["C2 death region is exactly p >= al/be", Resolve[ForAll[{al, be, p}, cons && 0 <= p < 1, Equivalent[lamPhi <= 0, p >= al/be]], Reals]];
ws["C3 death at finite time iff be > al", Resolve[ForAll[{al, be}, cons, Equivalent[Exists[p, 0 <= p < 1 && lamPhi <= 0], be > al]], Reals]];
ws["C4 z branch never fires on Phi", Resolve[ForAll[{al, be, p}, cons && 0 <= p < 1, 0 - Sqrt[rphi[[1, 1]]*rphi[[4, 4]]] <= 0], Reals]];
ws["C5 Psi family never dies", Resolve[ForAll[{al, be, p}, cons && 0 <= p < 1, rpsi[[2, 3]] - Sqrt[rpsi[[1, 1]]*rpsi[[4, 4]]] > 0], Reals]];
ws["C6 no doubly excited population is ever created", Simplify[rpsi[[4, 4]]] === 0];
ws["C7 transversal crossing at the death point", Simplify[(D[lamPhi, p] /. p -> al/be) + be*(be - al), Assumptions -> cons] === 0];
ws["C8 evolved Phi stays a state", Resolve[ForAll[{al, be, p}, cons && 0 <= p <= 1, rphi[[1, 1]] >= 0 && rphi[[2, 2]] >= 0 && rphi[[4, 4]] >= 0 && rphi[[1, 1]]*rphi[[4, 4]] - rphi[[1, 4]]^2 >= 0], Reals]];
ws["C9 death time formula", Simplify[1 - Exp[-g*(Log[be/(be - al)]/g)] - al/be, Assumptions -> 0 < al < be < 1 && g > 0] === 0];
genconsOuter = a >= 0 && b >= 0 && c >= 0 && d >= 0 && w > Sqrt[b*c] && w^2 <= a*d;
ws["C10 outer death criterion at the common fraction", Resolve[ForAll[{a, b, c, d, w}, genconsOuter, Equivalent[Exists[p, 0 <= p < 1 && w - Sqrt[(b + p*d)*(c + p*d)] <= 0], w^2 < (b + d)*(c + d)]], Reals]];
ws["C11 the outer crossing is permanent", Resolve[ForAll[{b, c, d, w, p1, p2}, b >= 0 && c >= 0 && d >= 0 && w >= 0 && 0 <= p1 < p2 < 1, (w - Sqrt[(b + p2*d)*(c + p2*d)]) <= (w - Sqrt[(b + p1*d)*(c + p1*d)])], Reals]];
ws["C12 balanced Bell state has C = (1-p)^2", Simplify[2*(lamPhi /. {al -> 1/Sqrt[2], be -> 1/Sqrt[2]}) - (1 - p)^2] === 0];
ws["C13 evolved Psi inner coherence is al be (1-p)", Simplify[rpsi[[2, 3]] - al*be*(1 - p)] === 0];
deriv = D[w - Sqrt[(b + pp*d)*(c + pp*d)], pp] /. pp -> p;
ws["C14 general transversality at any crossing with d > 0", Resolve[ForAll[{b, c, d, w, p}, b >= 0 && c >= 0 && d > 0 && w > 0 && 0 <= p < 1 && w == Sqrt[(b + p*d)*(c + p*d)], deriv < 0], Reals]];
innercons = a >= 0 && b >= 0 && c >= 0 && d >= 0 && a + b + c + d == 1 && z >= 0 && z^2 > a*d && z^2 <= b*c;
ws["C15 inner death criterion z^2 < d (trace one)", Resolve[ForAll[{a, b, c, d, z}, innercons, Equivalent[Exists[p, 0 <= p < 1 && z^2 <= d*(a + p*(b + c) + p^2*d)], z^2 < d]], Reals]];
ws["C16 inner radicand nondecreasing", Resolve[ForAll[{a, b, c, d, p1, p2}, a >= 0 && b >= 0 && c >= 0 && d >= 0 && 0 <= p1 < p2 < 1, d*(a + p1*(b + c) + p1^2*d) <= d*(a + p2*(b + c) + p2^2*d)], Reals]];
wit = xnum[0, 9/20, 9/20, 1/10, 0, 9/20];
ws["C17 survivor witness is a valid inner-entangled state", And[Tr[wit] === 1, Min[Eigenvalues[wit]] >= 0, (9/20)^2 <= (9/20)*(9/20), (9/20)^2 > 0*(1/10)]];
ws["C18 witness survives forever despite d > 0", Resolve[ForAll[p, 0 <= p < 1, (9/20)^2 > (1/10)*(0 + p*(9/10) + p^2/10)], Reals]];
ws["C19 sigma_x x sigma_x exchanges the superposition weights", KroneckerProduct[sx, sx].{al, 0, 0, be} === {be, 0, 0, al}];

(* ---- C'. unequal damping fractions (S2's rate-free clause) ---- *)
k0p[q_] := {{1, 0}, {0, Sqrt[1 - q]}}; k1p[q_] := {{0, Sqrt[q]}, {0, 0}};
kraus12 = Flatten[Table[KroneckerProduct[ki, kj], {ki, {k0p[p1], k1p[p1]}}, {kj, {k0p[p2], k1p[p2]}}], 1];
chan12[rh_] := Total[(#.rh.Transpose[#]) & /@ kraus12];
x12 = Simplify[chan12[xstate]];
outer12 = x12[[1, 4]] - Sqrt[x12[[2, 2]]*x12[[3, 3]]];
claimOuter12 = Sqrt[(1 - p1)*(1 - p2)]*(w - Sqrt[(b + p1*d)*(c + p2*d)]);
ws["Cu1 outer branch with unequal fractions", Simplify[outer12 - claimOuter12, Assumptions -> 0 <= p1 < 1 && 0 <= p2 < 1 && a >= 0 && b >= 0 && c >= 0 && d >= 0 && w >= 0 && z >= 0] === 0];
inner12 = x12[[2, 3]] - Sqrt[x12[[1, 1]]*x12[[4, 4]]];
claimInner12 = Sqrt[(1 - p1)*(1 - p2)]*(z - Sqrt[d*(a + p2*b + p1*c + p1*p2*d)]);
ws["Cu2 inner branch with unequal fractions", Simplify[inner12 - claimInner12, Assumptions -> 0 <= p1 < 1 && 0 <= p2 < 1 && a >= 0 && b >= 0 && c >= 0 && d >= 0 && w >= 0 && z >= 0] === 0];
ws["Cu3 rate-free outer death criterion", Resolve[ForAll[{b, c, d, w}, b >= 0 && c >= 0 && d >= 0 && w >= 0 && w^2 > b*c, Equivalent[Exists[{p1, p2}, 0 <= p1 < 1 && 0 <= p2 < 1 && w^2 <= (b + p1*d)*(c + p2*d)], w^2 < (b + d)*(c + d)]], Reals]];
ws["Cu4 outer radicand monotone in each fraction", Resolve[ForAll[{b, c, d, p1, p2, q1, q2}, b >= 0 && c >= 0 && d >= 0 && 0 <= p1 <= q1 < 1 && 0 <= p2 <= q2 < 1, (b + p1*d)*(c + p2*d) <= (b + q1*d)*(c + q2*d)], Reals]];
ws["Cu5 inner radicand dominated by its diagonal value", Resolve[ForAll[{a, b, c, d, p1, p2, q}, a >= 0 && b >= 0 && c >= 0 && d >= 0 && 0 <= p1 <= q && 0 <= p2 <= q && q < 1, d*(a + p2*b + p1*c + p1*p2*d) <= d*(a + q*b + q*c + q^2*d)], Reals]];
(* Cu5 + C15: an unequal-fraction witness reduces to a diagonal one at q = Max[p1, p2], so the inner criterion is rate-free *)

(* ---- D. the PPT route (S1, C1) ---- *)
ptB[rh_] := ArrayFlatten[Map[Transpose, Partition[rh, {2, 2}], {2}]];
ws["D1 partial transpose swaps the two coherences", ptB[xstate] === {{a, 0, 0, z}, {0, b, w, 0}, {0, w, c, 0}, {z, 0, 0, d}}];
evPT = Join[Eigenvalues[{{a, z}, {z, d}}], Eigenvalues[{{b, w}, {w, c}}]];
ws["D2 concurrence and PPT vanish together on X states", Resolve[ForAll[{a, b, c, d, w, z}, statecons, Equivalent[w^2 > b*c || z^2 > a*d, Or @@ ((# < 0) & /@ evPT)]], Reals]];
evPTphi = Eigenvalues[ptB[rphi]];
ws["D3 PPT crossing lands on the death point", Resolve[ForAll[{al, be, p}, cons && 0 <= p < 1, Equivalent[Or @@ ((# < 0) & /@ evPTphi), lamPhi > 0]], Reals]];

(* ---- E. the dephasing direction (S2') ---- *)
q0 = Sqrt[(1 + s)/2]*id2; q1 = Sqrt[(1 - s)/2]*sz;
dkraus = Flatten[Table[KroneckerProduct[qi, qj], {qi, {q0, q1}}, {qj, {q0, q1}}], 1];
dchan[rh_] := Total[(#.rh.Transpose[#]) & /@ dkraus];
dmapClaim = {{a, 0, 0, s^2*w}, {0, b, s^2*z, 0}, {0, s^2*z, c, 0}, {s^2*w, 0, 0, d}};
ws["E1 dephasing map on the X class", Simplify[dchan[xstate] - dmapClaim] === zero4];
dxt = dmapClaim /. s -> Exp[-h*t];
ws["E2 the dephasing map solves its Lindblad equation", FullSimplify[D[dxt, t] - (h/2)*(dis[szA][dxt] + dis[szB][dxt]), Assumptions -> h > 0 && t >= 0] === zero4];
ws["E3 no pure superposition dies under dephasing", Resolve[ForAll[{al, be, sv}, cons && 0 < sv <= 1, sv^2*al*be > 0], Reals]];
rF = {{f/2, 0, 0, (f/2)*sv^2}, {0, (1 - f)/2, 0, 0}, {0, 0, (1 - f)/2, 0}, {(f/2)*sv^2, 0, 0, f/2}};
ws["E4 the mixed witness is a state", Resolve[ForAll[f, 0 <= f <= 1, And @@ ((# >= 0) & /@ Eigenvalues[rF /. sv -> 1])], Reals]];
ws["E5 witness entangled at t = 0 iff F > 1/2", Resolve[ForAll[f, 0 < f < 1, Equivalent[Max[f/2 - (1 - f)/2, 0 - f/2] > 0, f > 1/2]], Reals]];
ws["E6 witness death point (F strictly below one)", Resolve[ForAll[{f, sv}, 1/2 < f < 1 && 0 < sv <= 1, Equivalent[(f/2)*sv^2 - (1 - f)/2 <= 0, sv^2 <= (1 - f)/f]], Reals]];
ws["E7 that death point is at finite time", Resolve[ForAll[f, 1/2 < f < 1, 0 < (1 - f)/f < 1], Reals]];
ws["E8 populations are frozen by dephasing", Simplify[{dmapClaim[[1, 1]] - a, dmapClaim[[2, 2]] - b, dmapClaim[[3, 3]] - c, dmapClaim[[4, 4]] - d}] === {0, 0, 0, 0}];
ws["E9 inner sector with ad = 0 never dies under dephasing", Resolve[ForAll[{zv, sv}, zv > 0 && 0 < sv <= 1, sv^2*zv - 0 > 0], Reals]];
ws["E10 sector death criterion (z = 0, bc > 0)", Resolve[ForAll[{b, c, w, sv}, b > 0 && c > 0 && w > Sqrt[b*c] && 0 < sv <= 1, Equivalent[sv^2*w - Sqrt[b*c] <= 0, sv^2 <= Sqrt[b*c]/w]], Reals]];
ws["E11 that death point is at finite time", Resolve[ForAll[{b, c, w}, b > 0 && c > 0 && w > Sqrt[b*c], 0 < Sqrt[b*c]/w < 1], Reals]];
ws["E12 bc = 0 outer sector never dies under dephasing", Resolve[ForAll[{w, sv}, w > 0 && 0 < sv <= 1, sv^2*w > 0], Reals]];

(* ---- F. the collective corner (S0') ---- *)
jm = smA + smB;
psim = {0, 1, -1, 0}/Sqrt[2];
rm = Outer[Times, psim, psim];
ws["F1 the singlet is annihilated by the collective lowering operator", jm.psim === {0, 0, 0, 0}];
ws["F2 the singlet is stationary under collective loss", Simplify[dis[jm][rm]] === zero4];
hsF = -(om1/2)*szA - (om2/2)*szB;
lfullF = -I*(hsF.rm - rm.hsF) + gg*dis[jm][rm];
ws["F3 singlet stationary with splittings iff om1 == om2", Reduce[And @@ Thread[Flatten[Simplify[lfullF]] == 0], {om1, om2}]];

(* ---- G. local energy splittings (N) ---- *)
rr4 = Array[rv, {4, 4}];
ham = -(om1/2)*szA - (om2/2)*szB;
lh[rh_] := -I*(ham.rh - rh.ham);
ld[rh_] := dis[smA][rh] + dis[smB][rh];
ws["G1 the splitting commutes with the loss generator", Simplify[lh[ld[rr4]] - ld[lh[rr4]]] === zero4];
diagvals = Diagonal[ham];
uu = DiagonalMatrix[Exp[-I*t*diagvals]]; uui = DiagonalMatrix[Exp[I*t*diagvals]];
rHt = uu.rhot.uui;
ws["G2 rotated damped solution solves the full equation", FullSimplify[D[rHt, t] - (lh[rHt] + g*(dis[smA][rHt] + dis[smB][rHt])), Assumptions -> g > 0 && t >= 0] === zero4];
ws["G3 coherence moduli unchanged by the splitting", Simplify[{rHt[[1, 4]]*rHt[[4, 1]] - rhot[[1, 4]]*rhot[[4, 1]], rHt[[2, 3]]*rHt[[3, 2]] - rhot[[2, 3]]*rhot[[3, 2]]}] === {0, 0}];
ws["G4 excited level sits above the ground level by om1", Simplify[({0, 1}.(-(om1/2)*sz).{0, 1}) - ({1, 0}.(-(om1/2)*sz).{1, 0})] === om1];

(* ---- X. review-round boundaries and limits (S2, H) ---- *)
ws["X1 boundary state is valid and entangled", And[1/2 + 1/20 + 1/10 + 7/20 === 1, (1/Sqrt[10])^2 <= (1/2)*(7/20), (1/Sqrt[10])^2 > (1/20)*(1/10)]];
ws["X2 boundary state satisfies the death criterion", (1/10 < (1/20 + 7/20)*(1/10 + 7/20)) === True];
ws["X3 boundary state dies when both fractions run to one", Resolve[Exists[{p1, p2}, 0 <= p1 < 1 && 0 <= p2 < 1 && 1/10 <= (1/20 + p1*(7/20))*(1/10 + p2*(7/20))], Reals]];
ws["X4 boundary state survives forever with the second bath off", Resolve[ForAll[p1, 0 <= p1 < 1, 1/10 > (1/20 + p1*(7/20))*(1/10)], Reals]];
ws["X5 outer radicand limit is (b+d)(c+d)", Simplify[((b + p1*d)*(c + p2*d) /. {p1 -> 1, p2 -> 1}) - (b + d)*(c + d)] === 0];
ws["X6 inner radicand limit is d at trace one", Simplify[(d*(a + p2*b + p1*c + p1*p2*d) /. {p1 -> 1, p2 -> 1}) - d*(a + b + c + d)] === 0];
ws["X7 evolved X marginals are diagonal at all times", {mapClaim[[1, 3]] + mapClaim[[2, 4]], mapClaim[[1, 2]] + mapClaim[[3, 4]]} === {0, 0}];
derivInner = D[z - Sqrt[d*(a + pp*(b + c) + pp^2*d)], pp] /. pp -> p;
ws["X8 inner transversality at any crossing with d > 0 (positivity kept)", Resolve[ForAll[{a, b, c, d, z, p}, a >= 0 && b >= 0 && c >= 0 && d > 0 && z > 0 && z^2 <= b*c && 0 <= p < 1 && z == Sqrt[d*(a + p*(b + c) + p^2*d)], derivInner < 0], Reals]];
ws["X9 the more entangled Phi state can be the dying one", Resolve[ForAll[{b1, b2}, 0 < b1 < 1/Sqrt[2] && 1/Sqrt[2] < b2 < Sqrt[1 - b1^2], 2*b2*Sqrt[1 - b2^2] > 2*b1*Sqrt[1 - b1^2]], Reals]];
uphase = KroneckerProduct[DiagonalMatrix[{1, Exp[I*ph1]}], DiagonalMatrix[{1, Exp[I*ph2]}]];
xgen = {{a, 0, 0, w*Exp[I*thw]}, {0, b, z*Exp[I*thz], 0}, {0, z*Exp[-I*thz], c, 0}, {w*Exp[-I*thw], 0, 0, d}};
rot = uphase.xgen.ConjugateTranspose[uphase];
ws["X10 local phases set both coherences real nonnegative", Simplify[{rot[[1, 4]] /. {ph1 -> (thw + thz)/2, ph2 -> (thw - thz)/2}, rot[[2, 3]] /. {ph1 -> (thw + thz)/2, ph2 -> (thw - thz)/2}} - {w, z}, Assumptions -> Element[{thw, thz}, Reals]] === {0, 0}];
ws["X11 XX conjugation permutes the X parameters a<->d, b<->c", Simplify[KroneckerProduct[sx, sx].xstate.KroneckerProduct[sx, sx] - {{d, 0, 0, w}, {0, c, z, 0}, {0, z, b, 0}, {w, 0, 0, a}}] === zero4];
spA = KroneckerProduct[Transpose[sm], id2]; spB = KroneckerProduct[id2, Transpose[sm]];
ws["X12 XX conjugates the raising generator into the lowering one", Simplify[KroneckerProduct[sx, sx].(dis[spA][KroneckerProduct[sx, sx].xstate.KroneckerProduct[sx, sx]] + dis[spB][KroneckerProduct[sx, sx].xstate.KroneckerProduct[sx, sx]]).KroneckerProduct[sx, sx] - (dis[smA][xstate] + dis[smB][xstate])] === zero4];
rhoColl = {{1 - be^2*Exp[-2*g*t] - 2*be^2*g*t*Exp[-2*g*t], 0, 0, al*be*Exp[-g*t]}, {0, be^2*g*t*Exp[-2*g*t], be^2*g*t*Exp[-2*g*t], 0}, {0, be^2*g*t*Exp[-2*g*t], be^2*g*t*Exp[-2*g*t], 0}, {al*be*Exp[-g*t], 0, 0, be^2*Exp[-2*g*t]}};
ws["X13 shared-bath closed form solves the collective equation on Phi", FullSimplify[D[rhoColl, t] - g*dis[jm][rhoColl], Assumptions -> g > 0 && t >= 0] === zero4];
ws["X14 shared-bath trace stays one", FullSimplify[Tr[rhoColl]] === 1];
ws["X15 shared-bath outer branch factor", FullSimplify[(rhoColl[[1, 4]] - Sqrt[rhoColl[[2, 2]]*rhoColl[[3, 3]]]) - be*(al*Exp[g*t] - be*g*t)*Exp[-2*g*t], Assumptions -> g > 0 && t >= 0 && 0 < al < 1 && 0 < be < 1] === 0];
ws["X16 the maximum of u Exp[-u] is 1/E at u = 1 (derivative sign pieces)", {Simplify[D[u*Exp[-u], u] - (1 - u)*Exp[-u]] === 0, Simplify[(1 - u)*Exp[-u] > 0, Assumptions -> u < 1] === True, Simplify[(1 - u)*Exp[-u] < 0, Assumptions -> u > 1] === True, 1*Exp[-1] === 1/E}];
collConc[tv_] := 2*Max[0, (rhoColl[[1, 4]] - Sqrt[rhoColl[[2, 2]]*rhoColl[[3, 3]]]) /. {al -> 1/10, be -> Sqrt[99]/10, g -> 1, t -> tv}, (rhoColl[[2, 3]] - Sqrt[rhoColl[[1, 1]]*rhoColl[[4, 4]]]) /. {al -> 1/10, be -> Sqrt[99]/10, g -> 1, t -> tv}];
ws["X17 dark interval then revival at one exact member", {collConc[3/10] === 0, collConc[1] === 0, collConc[8] > 0}];
