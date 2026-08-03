#!/usr/bin/env wolframscript
(* Verification probe (NOT the deliverable): confirm 4 load-bearing claims about
   the notebook's 𝒯 before rewriting it. Pure symbolic/structural; no QF. *)

out[s_] := WriteString[$Output, s, "\n"];
out["=== (a) DIMENSION: IdentityMatrix[2] vs d=4 Hamiltonian ==="];
d4H = RandomReal[1, {4, 4}];
badHeff = IdentityMatrix[2] - 0.1 d4H;              (* current code: hardcoded 2 *)
fixHeff = IdentityMatrix[Length[d4H]] - 0.1 d4H;    (* proposed fix: Length[Hmat] *)
out["  hardcoded IdentityMatrix[2] gives a valid 4x4 matrix? " <> ToString[MatrixQ[badHeff] && Dimensions[badHeff] === {4, 4}]];
out["  IdentityMatrix[Length[Hmat]] gives a valid 4x4 matrix? " <> ToString[MatrixQ[fixHeff] && Dimensions[fixHeff] === {4, 4}]];

out[""];
out["=== (b) MIDDLE-OPTIONAL ARGUMENT BINDING of T[ρ,H,L,η:None,V:None,δt,tf] ==="];
ClearAll[g];
g[r_, h_, l_, eta_ : None, v_ : None, dt_, tf_] := {r, h, l, eta, v, dt, tf};
b5 = g[r, h, l, dt, tf];          (* 5 args: the notebook's call shape *)
b6 = g[r, h, l, eta, dt, tf];     (* 6 args: + efficiency *)
b7 = g[r, h, l, eta, v, dt, tf];  (* 7 args: + unmonitored V *)
out["  5-arg  T[ρ,H,L,δt,tf]      -> " <> ToString[b5, InputForm]];
out["    expected {r,h,l,None,None,dt,tf}: " <> ToString[b5 === {r, h, l, None, None, dt, tf}]];
out["  6-arg  T[ρ,H,L,η,δt,tf]    -> " <> ToString[b6, InputForm]];
out["    expected {r,h,l,eta,None,dt,tf}:  " <> ToString[b6 === {r, h, l, eta, None, dt, tf}]];
out["  7-arg  T[ρ,H,L,η,V,δt,tf]  -> " <> ToString[b7, InputForm]];
out["    expected {r,h,l,eta,v,dt,tf}:     " <> ToString[b7 === {r, h, l, eta, v, dt, tf}]];

out[""];
out["=== (c) M_code vs M_paper (single monitored L, general η): expect diff = -1/2 (1-η) δt L.L ==="];
Lm = {{l11, l12}, {l21, l22}};
Hm = {{h11, h12}, {Conjugate[h12], h22}};
Ld = ConjugateTranspose[Lm];
(* Code's construction (single operator, no V): Heff + Υ + 1/2 Υ.Υ, Υ = √η δy L *)
HeffCode = IdentityMatrix[2] - dt (I Hm + 1/2 (Lm.Lm + Ld.Lm));
Ups = Sqrt[eta] dy Lm;
Mcode = HeffCode + Ups + 1/2 Ups.Ups;
(* Paper Eq (Mn2), single L: I - δt(iH + ½L†L) + √η L δy + ½ η L² (δy² - δt) *)
Mpaper = IdentityMatrix[2] - dt (I Hm + 1/2 Ld.Lm) + Sqrt[eta] dy Lm + 1/2 eta (dy^2 - dt) Lm.Lm;
diffC = FullSimplify[Mcode - Mpaper];
predC = FullSimplify[-(1/2) (1 - eta) dt Lm.Lm];
out["  M_code - M_paper simplifies to -1/2 (1-η) δt L.L ? " <> ToString[FullSimplify[diffC - predC] === {{0, 0}, {0, 0}}]];
out["  (so at η=1 they are identical: diff at η->1 = " <> ToString[FullSimplify[diffC /. eta -> 1], InputForm] <> ")"];

out[""];
out["=== (d) Spurious V^2 term in Heff (single unmonitored V) ==="];
Vm = {{v11, v12}, {v21, v22}};
Vd = ConjugateTranspose[Vm];
(* Code adds to Heff: -δt (1/2)(V.V + V†.V). Paper's M has only -δt (1/2) V†V *)
codeVpart = -dt (1/2) (Vm.Vm + Vd.Vm);
paperVpart = -dt (1/2) Vd.Vm;
diffV = FullSimplify[codeVpart - paperVpart];
out["  code V-part - paper V-part = -1/2 δt V.V ? " <> ToString[FullSimplify[diffV - (-(1/2) dt Vm.Vm)] === {{0, 0}, {0, 0}}]];

out[""];
out["=== (e) MULTI-OPERATOR double-sum collapse, TWO NON-COMMUTING monitored ops, general eta ==="];
L1 = {{a11, a12}, {a21, a22}};
L2 = {{b11, b12}, {b21, b22}};
noncomm = FullSimplify[L1 . L2 - L2 . L1];
S2 = Sqrt[et1] dy1 L1 + Sqrt[et2] dy2 L2;                       (* code's S *)
Heff2 = IdentityMatrix[2] - dt (I Hm + 1/2 (ConjugateTranspose[L1] . L1 + ConjugateTranspose[L2] . L2));
Lcorr2 = (dt/2) (et1 L1 . L1 + et2 L2 . L2);
Mcode2 = Heff2 + S2 + 1/2 S2 . S2 - Lcorr2;                     (* what the code builds *)
(* paper Eq (Mn2): the raw double sum 1/2 Sum_{r,s} Sqrt[et_r et_s] L_r L_s (dy_r dy_s - delta_rs dt) *)
dsum2 = 1/2 (et1 (dy1^2 - dt) L1 . L1 + et2 (dy2^2 - dt) L2 . L2
   + Sqrt[et1] Sqrt[et2] dy1 dy2 (L1 . L2 + L2 . L1));   (* Sqrt[et_r et_s] = Sqrt[et_r] Sqrt[et_s], et>0 *)
Mpaper2 = (IdentityMatrix[2] - dt (I Hm + 1/2 (ConjugateTranspose[L1] . L1 + ConjugateTranspose[L2] . L2))
   + Sqrt[et1] dy1 L1 + Sqrt[et2] dy2 L2 + dsum2);
diff2 = FullSimplify[Mcode2 - Mpaper2];
out["  L1, L2 chosen non-commuting (L1.L2 != L2.L1)? " <> ToString[noncomm =!= {{0, 0}, {0, 0}}]];
out["  code's (Heff+S+1/2 S.S-Lcorr) === paper's raw double sum, keeping L1 L2 order? " <> ToString[diff2 === {{0, 0}, {0, 0}}]];
out["  numeric residual at a random point: " <> ToString[Max[Abs[Flatten[(Mcode2 - Mpaper2) /. {a11 -> 1.2, a12 -> 0.3 + 0.1 I, a21 -> -0.4 + 0.2 I, a22 -> 0.7, b11 -> 0.2, b12 -> -0.5 I, b21 -> 0.6, b22 -> 0.1, et1 -> 0.7, et2 -> 0.4, dy1 -> 0.13, dy2 -> -0.09, dt -> 0.01, h11 -> 1., h12 -> 0.2, h22 -> -0.5}]]]]];

out[""];
out["=== SUMMARY ==="];
out["  (a) dimension fix needed: " <> ToString[!(MatrixQ[badHeff] && Dimensions[badHeff] === {4,4}) && (MatrixQ[fixHeff] && Dimensions[fixHeff]==={4,4})]];
out["  (b) arg-binding as claimed: " <> ToString[b5 === {r,h,l,None,None,dt,tf} && b6 === {r,h,l,eta,None,dt,tf} && b7 === {r,h,l,eta,v,dt,tf}]];
out["  (c) η discrepancy confirmed: " <> ToString[FullSimplify[diffC - predC] === {{0,0},{0,0}} && FullSimplify[diffC /. eta->1] === {{0,0},{0,0}}]];
out["  (d) spurious V^2 confirmed: " <> ToString[FullSimplify[diffV - (-(1/2) dt Vm.Vm)] === {{0,0},{0,0}}]];
out["  (e) two-op non-commuting collapse holds: " <> ToString[noncomm =!= {{0,0},{0,0}} && diff2 === {{0,0},{0,0}}]];
