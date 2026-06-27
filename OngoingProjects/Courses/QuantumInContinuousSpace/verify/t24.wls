\[Phi][x_] := Sech[x]/Sqrt[2];
\[Chi][x_] := Sech[x - 1] Exp[I x]/Sqrt[2];
(* Approach A: matrix-element difference as a single total-derivative integral
   <phi|p chi> - <p phi|chi> = -i Integrate[ d/dx (phi^* chi) ] = -i [phi^* chi]_{-inf}^{inf} *)
A = TimeConstrained[Integrate[-I D[Conjugate[\[Phi][x]] \[Chi][x], x], {x, -Infinity, Infinity}], 60, "TIMEOUT"];
WriteString[$Output, "A total-derivative integral => ", ToString[A, InputForm], "\n"];
(* Approach B: numeric confirmation of the two-integral difference *)
B = Chop[NIntegrate[Conjugate[\[Phi][x]] (-I D[\[Chi][x], x]), {x, -Infinity, Infinity}, WorkingPrecision -> 30] -
        NIntegrate[Conjugate[-I D[\[Phi][x], x]] \[Chi][x], {x, -Infinity, Infinity}, WorkingPrecision -> 30], 10^-15];
WriteString[$Output, "B numeric difference (Chop) => ", ToString[B, InputForm], "\n"];
(* Approach C: boundary term via Limits *)
C0 = Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> Infinity] - Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> -Infinity];
WriteString[$Output, "C boundary term => ", ToString[C0, InputForm], "\n"];
WriteString[$Output, "DONE\n"];
