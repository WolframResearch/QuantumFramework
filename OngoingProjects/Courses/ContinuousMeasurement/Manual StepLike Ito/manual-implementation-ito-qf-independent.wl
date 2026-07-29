(* ::Package:: *)

(* Positivity-preserving integrator for the diffusive stochastic master
   equation, independent of QuantumFramework.

   Scheme: Rouchon and Ralph, "Efficient Quantum Filtering for Quantum Feedback
   Control", Phys. Rev. A 91, 012118 (2015); arXiv:1410.5345, Eqs. (3)-(4).

     rho_c(n+1) = N[ M rho M^dag  +  dt Sum_j V_j rho V_j^dag
                     +  dt Sum_r (1-eta_r) L_r rho L_r^dag ] ,
     M = I - (i H + 1/2 Sum_j V_j^dag V_j + 1/2 Sum_r L_r^dag L_r) dt
           + Sum_r Sqrt[eta_r] L_r dy_r
           + 1/2 Sum_{r,s} Sqrt[eta_r eta_s] L_r L_s (dy_r dy_s - delta_rs dt) ,
     dy_r = Sqrt[eta_r] Tr[(L_r + L_r^dag) rho] dt + dW_r ,

   where N[.] is trace normalization. Because the numerator is a sum of terms of
   the form A rho A^dag, the conditioned state stays positive semidefinite for
   any step size dt.

   Public symbols (script letters, so this drops into the source notebook once
   operators are supplied as plain matrices):
     \[ScriptCapitalR]  per-step readout increments dy   (list, one entry per monitored operator)
     \[ScriptCapitalT]  conditional trajectory            (list of density matrices, seed included)

   Nothing here touches QuantumFramework. States and operators are ordinary
   dense complex matrices. Extract them once from whatever source you like
   (for instance qo["Computational"]["Matrix"]) and pass them in.

   \[ScriptCapitalT][rho0, H, Ls, eta, Vs, dt, tf]
     rho0  d x d density matrix (initial conditioned state)
     H     d x d Hamiltonian
     Ls    list of d x d monitored jump operators (already scaled, e.g. Sqrt[gamma] L)
     eta   list of detector efficiencies in (0,1], or None -> all 1 (optional)
     Vs    list of d x d unmonitored Lindblad operators, or None (optional)
     dt    time increment
     tf    final time; the trajectory has Floor[tf/dt]+1 states

   The per-step readout vector dy is Sow-ed, so Reap[\[ScriptCapitalT][...]] returns
   {trajectory, {records}} with records[[n]] the dy of step n. *)

ClearAll[\[ScriptCapitalR], \[ScriptCapitalT], rouchonReadout, rouchonTrajectory];

(* Readout increments over one step, matrices in, list out. Kept as a standalone
   for API parity with the source notebook; \[ScriptCapitalT] inlines the same expression. *)
\[ScriptCapitalR][rho_, Ls_List, dt_, dw_List, eta_List] :=
  MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {Ls, dw, eta}];
rouchonReadout = \[ScriptCapitalR];

\[ScriptCapitalT][rho0_?MatrixQ, ham_?MatrixQ, Ls_List, eta_ : None, Vs_ : None, dt_, tf_] :=
 Module[
  {d = Length[ham], nL = Length[Ls], etav, sqeta, etaL, lsum, ldl, vdv,
   dtN, heff, lcorr, ineff, vchan, combined, extra, nsteps, dw, step},

  etav = If[eta === None, ConstantArray[1, nL], eta];
  dtN  = N[dt];

  (* --- step-invariant pieces, computed once --- *)
  sqeta = Sqrt[N[etav]];                                       (* Sqrt[eta_r] *)
  etaL  = N @ MapThread[#1 #2 &, {sqeta, Ls}];                 (* Sqrt[eta_r] L_r *)
  lsum  = N @ Map[# + ConjugateTranspose[#] &, Ls];            (* L_r + L_r^dag *)
  ldl   = Total[Map[ConjugateTranspose[#] . # &, Ls]];         (* Sum_r L_r^dag L_r *)
  vdv   = If[Vs === None, 0, Total[Map[ConjugateTranspose[#] . # &, Vs]]];
  heff  = N[IdentityMatrix[d] - dtN (I ham + (1/2) ldl + (1/2) vdv)];
  (* Ito correction 1/2 Sum_r eta_r L_r^2 dt (the -delta_rs dt piece of M).
     Weighted by eta_r, matching the paper for inefficient detection. *)
  lcorr = N[(dtN/2) Total[MapThread[#1 (#2 . #2) &, {etav, Ls}]]];

  (* dissipation restored for unmonitored (V) and unmonitored fraction (1-eta) *)
  ineff = MapThread[
     If[#2 < 1, {#1, ConjugateTranspose[#1], dtN (1 - #2)}, Nothing] &, {Ls, etav}];
  vchan = If[Vs === None, {}, Map[{#, ConjugateTranspose[#], dtN} &, Vs]];
  combined = N @ Join[ineff, vchan];
  extra = If[combined === {}, (0 &),
     With[{cc = combined},
      Function[r, Total[Map[(#[[3]] (#[[1]] . r . #[[2]])) &, cc]]]]];

  (* --- Wiener increments for the whole run --- *)
  nsteps = Floor[tf/dt];
  dw = RandomVariate[NormalDistribution[0, Sqrt[dtN]], {nL, nsteps}];

  (* --- one update, all constants captured (no pattern matching in the loop) --- *)
  step = With[
    {heff = heff, etaL = etaL, lsum = lsum, sqeta = sqeta,
     lcorr = lcorr, extra = extra, dtN = dtN},
    Function[{rho, dwv},
     Module[{dy, s, m, num},
      dy  = Sow @ MapThread[#2 Tr[#1 . rho] dtN + #3 &, {lsum, sqeta, dwv}];
      s   = dy . etaL;                              (* Sum_r dy_r Sqrt[eta_r] L_r *)
      m   = heff + s + (1/2) s . s - lcorr;
      num = m . rho . ConjugateTranspose[m] + extra[rho];
      num / Tr[num]
      ]]];

  FoldList[step, N[rho0], Transpose[dw]]
  ];
rouchonTrajectory = \[ScriptCapitalT];
