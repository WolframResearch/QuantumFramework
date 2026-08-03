(* ::Package:: *)

(* Prototypes accompanying SDE-schemes-literature-audit.md. QuantumFramework-
   independent throughout: states and operators are plain dense complex vectors
   and matrices, so the file loads and runs from a bare wolframscript.

   Read against the density-matrix Rouchon-Ralph baseline in
   manual-implementation-ito-qf-independent.wl (public symbol \[ScriptCapitalT]).

   Three integrators for the diffusive stochastic master equation:

     \[ScriptCapitalP]  pure-state (state-vector) Rouchon map. Valid when detection is
            perfect (eta = 1), every dissipation channel is monitored, and the
            initial state is pure. It produces the SAME conditioned trajectory as
            \[ScriptCapitalT] (a rank-one density matrix at every step), because
              M rho M^dag / Tr[M rho M^dag]  with rho = |psi><psi|
            equals |phi><phi| for phi = M psi / ||M psi||. The step never forms the
            d x d operator M: it applies the pieces to the ket,
              M psi = Heff psi + (S psi) + 1/2 S (S psi) - Lcorr psi ,
            a fixed number of matrix-vector products (O(d^2) each), instead of the
            matrix-matrix products M rho M^dag (O(d^3)) of the density-matrix step;
            and it stores d numbers instead of d^2. Positivity is automatic: a
            normalized ket is a valid pure state at any dt. This is the
            asymptotically-optimal choice for the efficient-filtering regime the
            Rouchon-Ralph title names.

     \[ScriptCapitalF]  the same pure-state map, form-M variant: it forms M once (one
            d x d matmul, S.S) and applies it to the ket, rather than the pure O(d^2)
            matrix-vector chain of \[ScriptCapitalP]. Same conditioned trajectory. One
            O(d^3) product per step, but far fewer Wolfram Language evaluator
            dispatches, so it is faster than \[ScriptCapitalP] at qubit / qudit scale
            (where dispatch, not arithmetic, dominates); \[ScriptCapitalP] overtakes
            it only at large d. This is the recommended production map for
            finite-dimensional monitored systems.

     \[ScriptCapitalE]  Euler-Maruyama step of the nonlinear (trace-normalized) Ito SME,
            the first-order foil Rouchon-Ralph is measured against. First order like
            \[ScriptCapitalT], but NOT positivity preserving: the update is
            rho + drift dt + diffusion dW with no Kraus structure, so the state can
            acquire a negative eigenvalue at moderate dt. Included to exhibit the
            accuracy-vs-structure trade, never as a recommendation.

   Conventions match the baseline: monitored operators are already scaled
   (e.g. Sqrt[gamma] L); the readout dy_r = Tr[(L_r + L_r^dag) rho] dt + dW_r is
   Sow-ed each step; dW is drawn as RandomVariate[..., {nL, nsteps}] then
   Transposed, so a shared SeedRandom reproduces \[ScriptCapitalT]'s noise
   realization exactly and the maps can be compared trajectory by
   trajectory. *)

ClearAll[\[ScriptCapitalP], \[ScriptCapitalF], \[ScriptCapitalE],
   pureStateTrajectory, pureStateFormMTrajectory, eulerMaruyamaTrajectory, rouchonHeffLcorr];

(* step-invariant pieces {Heff, Lcorr} of the Rouchon map at eta=1, shared by the two
   pure-state maps so the Ito correction is written once. *)
rouchonHeffLcorr[ham_, Ls_, dt_] := {
   IdentityMatrix[Length[ham]] - dt (I ham + (1/2) Total[ConjugateTranspose[#] . # & /@ Ls]),
   (dt/2) Total[# . # & /@ Ls]};                            (* {Heff, Lcorr}, eta=1 *)

(* ---- pure-state Rouchon map, O(d^2) matrix-vector form (eta = 1, pure seed) ----

   \[ScriptCapitalP][psi0, H, Ls, dt, tf]
     psi0  length-d initial ket (need not be normalized; it is normalized once)
     H     d x d Hamiltonian
     Ls    list of d x d monitored jump operators (already scaled)
     dt    time increment
     tf    final time; the trajectory has Floor[tf/dt]+1 kets

   M is never formed: every per-step operation is a matrix-vector product, O(d^2). This
   is the asymptotically-optimal large-d choice; at small d WL evaluator dispatch, not
   flops, dominates, and \[ScriptCapitalF] below (one matmul, fewer dispatches) is faster.  *)

\[ScriptCapitalP][psi0_?VectorQ, ham_?MatrixQ, Ls_List, dt_?NumericQ, tf_?NumericQ] :=
 Block[
  {nL = Length[Ls], psi0n, hamn, Lsn, dtN, heff, lcorr, nsteps, dw, step},
  {psi0n, hamn, Lsn} = N@{psi0, ham, Ls};
  dtN = N[dt];
  {heff, lcorr} = rouchonHeffLcorr[hamn, Lsn, dtN];
  nsteps = Floor[tf/dt];
  dw = RandomVariate[NormalDistribution[0, Sqrt[dtN]], {nL, nsteps}];
  (* dy_r = Tr[(L_r+L_r^dag) rho] dt + dW_r = 2 Re<psi|L_r|psi> dt + dW_r reuses L_r psi. *)
  step = With[
    {heff = heff, Lsn = Lsn, lcorr = lcorr, dtN = dtN},
    Function[{psi, dwv},
     With[{lpsi = # . psi & /@ Lsn},                       (* {L_r psi} : nL matrix-vector products *)
      With[{dy = Sow @ MapThread[2 Re[Conjugate[psi] . #1] dtN + #2 &, {lpsi, dwv}]},
       With[{sp = dy . lpsi},                              (* S psi = Sum_r dy_r L_r psi *)
        Normalize[heff . psi + sp + (1/2) (dy . (# . sp & /@ Lsn)) - lcorr . psi]]]]]];
  FoldList[step, Normalize[psi0n], Transpose[dw]]
  ];
pureStateTrajectory = \[ScriptCapitalP];

(* ---- pure-state Rouchon map, form-M variant (eta = 1, pure seed) ----

   \[ScriptCapitalF][psi0, H, Ls, dt, tf], same signature and same conditioned trajectory
   as \[ScriptCapitalP]. It forms M = Heff + S + 1/2 S.S - Lcorr once (one d x d matmul,
   S.S) then applies it to the ket. That is one O(d^3) product per step instead of
   \[ScriptCapitalP]'s pure O(d^2), but far fewer evaluator dispatches, so in Wolfram
   Language it is the faster pure-state map at qubit / qudit scale (small-to-moderate d),
   where dispatch dominates arithmetic; \[ScriptCapitalP] overtakes it only once d is
   large enough that flops dominate. This is the recommended production map for
   finite-dimensional monitored systems.                                             *)

\[ScriptCapitalF][psi0_?VectorQ, ham_?MatrixQ, Ls_List, dt_?NumericQ, tf_?NumericQ] :=
 Block[
  {nL = Length[Ls], psi0n, hamn, Lsn, dtN, heff, lcorr, lsum, nsteps, dw, step},
  {psi0n, hamn, Lsn} = N@{psi0, ham, Ls};
  dtN = N[dt];
  {heff, lcorr} = rouchonHeffLcorr[hamn, Lsn, dtN];
  lsum = # + ConjugateTranspose[#] & /@ Lsn;               (* L_r + L_r^dag *)
  nsteps = Floor[tf/dt];
  dw = RandomVariate[NormalDistribution[0, Sqrt[dtN]], {nL, nsteps}];
  step = With[
    {heff = heff, Lsn = Lsn, lsum = lsum, lcorr = lcorr, dtN = dtN},
    Function[{psi, dwv},
     With[{dy = Sow @ MapThread[Re[Conjugate[psi] . (#1 . psi)] dtN + #2 &, {lsum, dwv}]},
      With[{s = dy . Lsn},                                 (* S = Sum_r dy_r L_r *)
       Normalize[(heff + s + (1/2) s . s - lcorr) . psi]]]]];
  FoldList[step, Normalize[psi0n], Transpose[dw]]
  ];
pureStateFormMTrajectory = \[ScriptCapitalF];

(* ---- Euler-Maruyama, nonlinear Ito SME (eta = 1): first order, NOT positive ----

   drho = (-i[H,rho] + Sum_r D[L_r](rho)) dt + Sum_r H[L_r](rho) dW_r ,
     D[L](rho)  = L rho L^dag - 1/2 (L^dag L rho + rho L^dag L)     (Lindblad)
     H[L](rho)  = L rho + rho L^dag - Tr[(L+L^dag) rho] rho          (measurement)
   trace-renormalized each step so the demonstrated failure is negativity alone,
   not trace drift.                                                              *)

\[ScriptCapitalE][rho0_?MatrixQ, ham_?MatrixQ, Ls_List, dt_?NumericQ, tf_?NumericQ] :=
 Block[
  {nL = Length[Ls], rho0n, hamn, Lsn, dtN, lsum, ldl, nsteps, dw, step},
  {rho0n, hamn, Lsn} = N@{rho0, ham, Ls};
  dtN = N[dt];
  lsum = # + ConjugateTranspose[#] & /@ Lsn;
  ldl  = ConjugateTranspose[#] . # & /@ Lsn;               (* per-operator L^dag L *)
  nsteps = Floor[tf/dt];
  dw = RandomVariate[NormalDistribution[0, Sqrt[dtN]], {nL, nsteps}];
  step = With[
    {hamn = hamn, Lsn = Lsn, lsum = lsum, ldl = ldl, dtN = dtN},
    Function[{rho, dwv},
     With[{dy = Sow @ MapThread[Re[Tr[#1 . rho]] dtN + #2 &, {lsum, dwv}]},
      With[{
        drift = -I (hamn . rho - rho . hamn)
          + Total[MapThread[
             #1 . rho . ConjugateTranspose[#1] - (1/2) (#2 . rho + rho . #2) &,
             {Lsn, ldl}]],
        diff = Total[MapThread[
           (#1 . rho + rho . ConjugateTranspose[#1] - Re[Tr[#2 . rho]] rho) #3 &,
           {Lsn, lsum, dwv}]]},
       With[{st = rho + drift dtN + diff},
        st/Re[Tr[st]]]]]]];
  FoldList[step, rho0n, Transpose[dw]]
  ];
eulerMaruyamaTrajectory = \[ScriptCapitalE];
