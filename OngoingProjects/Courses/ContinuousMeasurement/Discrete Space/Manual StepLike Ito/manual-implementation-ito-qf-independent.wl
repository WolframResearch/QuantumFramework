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

   Writing S = Sum_r dy_r Sqrt[eta_r] L_r, the double sum collapses to
   1/2 S^2 - 1/2 dt Sum_r eta_r L_r^2, so M = Heff + S + 1/2 S^2 - Lcorr with
   Heff = I - (iH + 1/2 Sum V^dag V + 1/2 Sum L^dag L) dt and the Ito correction
   Lcorr = 1/2 dt Sum_r eta_r L_r^2. That collapse (including the non-commuting
   L_r L_s ordering) is what the code below builds; it is checked symbolically
   for one and two non-commuting operators in the companion verify script.

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
   {trajectory, {records}} with records[[n]] the dy of step n. \[ScriptCapitalR]
   computes that same dy standalone (matrices in, list out). *)

ClearAll[\[ScriptCapitalR], \[ScriptCapitalT], rouchonReadout, rouchonTrajectory];

(* Readout increments over one step. dy_r = Sqrt[eta_r] Tr[(L_r+L_r^dag) rho] dt + dw_r.
   \[ScriptCapitalT] inlines this same expression (with L_r+L_r^dag precomputed);
   a harness test pins the two forms together so they cannot drift. *)
\[ScriptCapitalR][rho_?MatrixQ, Ls_List, dt_?NumericQ, dw_List, eta_List] :=
  MapThread[Sqrt[#3] Tr[(#1 + ConjugateTranspose[#1]) . rho] dt + #2 &, {Ls, dw, eta}];
rouchonReadout = \[ScriptCapitalR];

\[ScriptCapitalT][rho0_?MatrixQ, ham_?MatrixQ, Ls_List, eta : (None | _List) : None, Vs : (None | _List) : None, dt_?NumericQ, tf_?NumericQ] :=
 Block[
  {d = Length[ham], nL = Length[Ls], rho0n, hamn, Lsn, Vsn, etav, sqeta, etaL,
   lsum, ldl, vdv, dtN, heff, lcorr, ineff, vchan, combined, extra, nsteps, dw, step},

  (* numericalize the inputs once, up front *)
  {rho0n, hamn, Lsn} = N@{rho0, ham, Ls};
  Vsn = If[Vs === None, None, N@Vs];
  dtN = N[dt];
  etav = N@If[eta === None, ConstantArray[1, nL], eta];

  (* step-invariant pieces (Times is Listable: sqeta Lsn is {Sqrt[eta_r] L_r}) *)
  sqeta = Sqrt[etav];
  etaL  = sqeta Lsn;                                       (* Sqrt[eta_r] L_r *)
  lsum  = # + ConjugateTranspose[#] & /@ Lsn;              (* L_r + L_r^dag *)
  ldl   = Total[ConjugateTranspose[#] . # & /@ Lsn];       (* Sum_r L_r^dag L_r *)
  vdv   = If[Vsn === None, 0, Total[ConjugateTranspose[#] . # & /@ Vsn]];
  heff  = IdentityMatrix[d] - dtN (I hamn + (1/2) ldl + (1/2) vdv);
  lcorr = (dtN/2) Total[etav (# . # & /@ Lsn)];            (* 1/2 dt Sum_r eta_r L_r^2 *)

  (* dissipation restored for unmonitored (V) and the unmonitored fraction (1-eta) *)
  ineff = MapThread[
     If[#2 < 1, {#1, ConjugateTranspose[#1], dtN (1 - #2)}, Nothing] &, {Lsn, etav}];
  vchan = If[Vsn === None, {}, {#, ConjugateTranspose[#], dtN} & /@ Vsn];
  combined = Join[ineff, vchan];
  extra = If[combined === {}, (0 &),
     With[{cc = combined}, Function[r, Total[#[[3]] (#[[1]] . r . #[[2]]) & /@ cc]]]];

  nsteps = Floor[tf/dt];
  (* Draw {nL, nsteps} then Transpose (not {nsteps, nL} directly): this fixes the
     random-draw ORDER to the notebook's, so a given SeedRandom reproduces the
     notebook's trajectory bit for bit. One O(nL nsteps) transpose per call. *)
  dw = RandomVariate[NormalDistribution[0, Sqrt[dtN]], {nL, nsteps}];

  (* one update, all constants captured, no pattern matching or Module in the loop *)
  step = With[
    {heff = heff, etaL = etaL, lsum = lsum, sqeta = sqeta,
     lcorr = lcorr, extra = extra, dtN = dtN},
    Function[{rho, dwv},
     With[{dy = Sow @ MapThread[#2 Tr[#1 . rho] dtN + #3 &, {lsum, sqeta, dwv}]},
      With[{s = dy . etaL},                                (* Sum_r dy_r Sqrt[eta_r] L_r *)
       With[{m = heff + s + (1/2) s . s - lcorr},
        With[{num = m . rho . ConjugateTranspose[m] + extra[rho]},
         num/Tr[num]]]]]]];

  FoldList[step, rho0n, Transpose[dw]]
  ];
rouchonTrajectory = \[ScriptCapitalT];
