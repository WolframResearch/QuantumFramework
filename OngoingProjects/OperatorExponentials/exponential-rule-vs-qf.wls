PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName, 3], "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
(* QF ordered e^{x a^2 + y ad^2 + w ad a} against the rule on K+, K0, K- at k = 1/4: Fock elements m, n <= 3 at a 40-digit point. *)
out[s_, e_] := WriteString["stdout", s, " -> ", ToString[e, InputForm], "\n"];
v = {{ad}, {a}};
Mp = {{0, 1}, {0, 0}}; Mm = {{0, 0}, {-1, 0}}; M0 = {{1/2, 0}, {0, -1/2}};
(* G = x a^2 + y ad^2 + w ad a = 2x Km + 2y Kp + 2w K0 - w/2, Kp = ad^2/2, Km = a^2/2, K0 = (N + 1/2)/2 *)
abstract[x_, y_, w_] := With[{g = MatrixExp[2 y Mp + 2 w M0 + 2 x Mm]},
  With[{p = g[[1, 2]]/g[[2, 2]], r = -g[[2, 1]]/g[[2, 2]], g22 = g[[2, 2]]},
   Exp[-w/2] g22^(-1/2) Exp[p ad ** ad/2] ** Exp[Log[g22^(-2)]/2 ad ** a] ** Exp[r a ** a/2]]];
qf[x_, y_, w_] := BosonicExpOrder[Exp[x a ** a + y ad ** ad + w ad ** a], v];

pts = {N[{3/10 + I/5, -1/7 + I/3, 2/5 - I/4}, 40]};
Do[With[{x = pt[[1]], y = pt[[2]], w = pt[[3]]},
  out["max |QF - abstract| Fock m,n<=3 at " <> ToString[pt, InputForm],
   Max[Abs[N[Flatten[Table[BosonicMatrixElement[{m, n}, qf[x, y, w], v] - BosonicMatrixElement[{m, n}, abstract[x, y, w], v], {m, 0, 3}, {n, 0, 3}]], 30]]]]], {pt, pts}];
