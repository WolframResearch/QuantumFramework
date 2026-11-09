(* The exponential rule on K+, K0, K-: exact check with K- = d/dz, K0 = z d/dz + k, K+ = z^2 d/dz + 2 k z (symbolic k, coefficients, f), then the branch of g22^(-2k) against the exact flow up to t = 10. *)
out[s_, e_] := WriteString["stdout", s, " -> ", ToString[e, InputForm], "\n"];
(* holomorphic realization, symbolic k: Km = d/dz, K0 = z d/dz + k, Kp = z^2 d/dz + 2 k z *)
Kp[h_] := z^2 D[h, z] + 2 k z h; K0[h_] := z D[h, z] + k h; Km[h_] := D[h, z];
comm[A_, B_][h_] := A[B[h]] - B[A[h]];
out["[Km,Kp]-2K0", Simplify[comm[Km, Kp][f[z]] - 2 K0[f[z]]]];
out["[K0,Kp]-Kp", Simplify[comm[K0, Kp][f[z]] - Kp[f[z]]]];
out["[K0,Km]+Km", Simplify[comm[K0, Km][f[z]] + Km[f[z]]]];
out["Casimir - k(k-1)", Simplify[K0[K0[f[z]]] - K0[f[z]] - Kp[Km[f[z]]] - k (k - 1) f[z]]];
(* 2x2 matrix, g' = M g, det g = 1 *)
Mp = {{0, 1}, {0, 0}}; Mm = {{0, 0}, {-1, 0}}; M0 = {{1/2, 0}, {0, -1/2}};
out["2x2 [Km,Kp]=2K0", Mm.Mp - Mp.Mm == 2 M0];
M = la Mp + mu M0 + nu Mm;
g = {{g11[t], g12[t]}, {g21[t], g22[t]}};
p = g12[t]/g22[t]; r = -g21[t]/g22[t];
(* product action on f: e^{pKp} e^{qK0} e^{rKm} f with e^{qk} = g22^(-2k), e^q = g22^(-2) *)
u = (1 - p z)^(-2 k) g22[t]^(-2 k) f[g22[t]^(-2) z/(1 - p z) + r];
dg = Thread[Flatten[D[g, t]] -> Flatten[M.g]];
G[h_] := la Kp[h] + mu K0[h] + nu Km[h];
res = (D[u, t] - G[u]) /. dg /. g11[t] -> (1 + g12[t] g21[t])/g22[t];
out["d/dt u - G u (symbolic k, la, mu, nu, generic f)", Simplify[res/u /. f -> Function[x, Exp[c x]]]];
out["same, f = generic, ratio-free", Simplify[res]];
out["u at t=0 (g = 1)", u /. {g12[t] -> 0, g21[t] -> 0, g22[t] -> 1}];
Mp = {{0, 1}, {0, 0}}; Mm = {{0, 0}, {-1, 0}}; M0 = {{1/2, 0}, {0, -1/2}};
{la, mu, nu} = {-0.6 I, -2 I, -0.6 I}; (* -i t H, H = 2 K0 + 0.6 (K+ + K-): stable *)
M = la Mp + mu M0 + nu Mm;
f[x_] := 1 + x^2; z0 = 0.3 + 0.1 I;
prod[t_, k_] := With[{g = MatrixExp[t M]}, With[{p = g[[1, 2]]/g[[2, 2]], r = -g[[2, 1]]/g[[2, 2]]},
   (1 - p z0)^(-2 k) g[[2, 2]]^(-2 k) f[g[[2, 2]]^(-2) z0/(1 - p z0) + r]]];
exact[tmax_, k_] := First @ NDSolve[{Z'[s] == la Z[s]^2 + mu Z[s] + nu, Z[0] == z0, L'[s] == 2 k la Z[s] + mu k, L[0] == 0}, {Z, L}, {s, 0, tmax}, WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 20];
Do[With[{sol = exact[10, k]},
  out["k=" <> ToString[k], Table[{t, Round[Chop[prod[t, k]/(Exp[L[t]] f[Z[t]] /. sol), 10^-8], 0.001]}, {t, 0, 10, 1}]]], {k, {1/4, 3/4, 7/10}}];
out["g22 along t", Table[{t, Round[MatrixExp[t M][[2, 2]], 0.01]}, {t, 0, 10, 1}]];
