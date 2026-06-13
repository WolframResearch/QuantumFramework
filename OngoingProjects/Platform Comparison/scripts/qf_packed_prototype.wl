PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

(* ---- packed helpers (prototype, standalone) ---- *)
packRows[m_] := FromDigits[Reverse[#], 2] & /@ Transpose[m];   (* {nq,2n}->len-2n ints, bit q-1 = qubit q *)
unpackRows[ints_, n_] := Transpose[Reverse /@ IntegerDigits[ints, 2, n]];  (* len-2n ints -> {nq,2n} *)
getPacked[ps_] := With[{t = ps["Tableau"]}, {Dimensions[t][[2]], packRows[t[[1]]], packRows[t[[2]]], ps["Signs"]}];
tableauFromPacked[n_, px_, pz_] := {unpackRows[px, n], unpackRows[pz, n]};
fromPacked[{n_, px_, pz_, s_}] := PauliStabilizer[<|"Signs" -> s, "Tableau" -> tableauFromPacked[n, px, pz]|>];

shiftBit[v_, from_, to_] := With[{d = to - from}, If[d >= 0, BitShiftLeft[v, d], BitShiftRight[v, -d]]];

pH[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]},
   With[{bx = BitAnd[px, m], bz = BitAnd[pz, m]},
     {n, px - bx + bz, pz - bz + bx, s (1 - 2 Sign[BitAnd[bx, bz]])}]];
pS[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]},
   With[{bx = BitAnd[px, m], bz = BitAnd[pz, m]},
     {n, px, BitXor[pz, bx], s (1 - 2 Sign[BitAnd[bx, bz]])}]];
pSdg[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]},
   With[{bx = BitAnd[px, m], bz = BitAnd[pz, m]},
     {n, px, BitXor[pz, bx], s (1 - 2 Sign[BitAnd[bx, BitXor[bz, m]]])}]];
pCNOT[{n_, px_, pz_, s_}, j_, k_] := With[{mj = BitShiftLeft[1, j - 1], mk = BitShiftLeft[1, k - 1]},
   With[{bxj = BitAnd[px, mj], bzk = BitAnd[pz, mk],
         a = Sign[BitAnd[px, mj]], b = Sign[BitAnd[pz, mk]], c = Sign[BitAnd[px, mk]], d = Sign[BitAnd[pz, mj]]},
     {n, BitXor[px, shiftBit[bxj, j, k]], BitXor[pz, shiftBit[bzk, k, j]],
      s (1 - 2 (a b (1 - BitXor[c, d])))}]];
pSWAP[{n_, px_, pz_, s_}, j_, k_] := With[{mj = BitShiftLeft[1, j - 1], mk = BitShiftLeft[1, k - 1]},
   {n, BitXor[px, (mj + mk) BitXor[Sign[BitAnd[px, mj]], Sign[BitAnd[px, mk]]]],
       BitXor[pz, (mj + mk) BitXor[Sign[BitAnd[pz, mj]], Sign[BitAnd[pz, mk]]]], s}];
pX[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]}, {n, px, pz, s (1 - 2 Sign[BitAnd[pz, m]])}];
pZ[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]}, {n, px, pz, s (1 - 2 Sign[BitAnd[px, m]])}];
pY[{n_, px_, pz_, s_}, j_] := With[{m = BitShiftLeft[1, j - 1]}, {n, px, pz, s (1 - 2 BitXor[Sign[BitAnd[px, m]], Sign[BitAnd[pz, m]]])}];

(* ---- equivalence check: random Clifford streams, packed vs canonical ---- *)
SeedRandom[1234];
canon[ps_, op_] := Switch[op[[1]],
   "H", ps["H", op[[2]]], "S", ps["S", op[[2]]], "Sdg", ps[SuperDagger["S"], op[[2]]],
   "X", ps["X", op[[2]]], "Y", ps["Y", op[[2]]], "Z", ps["Z", op[[2]]],
   "SWAP", ps["SWAP", op[[2]], op[[3]]], "CNOT", ps["CNOT", op[[2]], op[[3]]]];
pk[st_, op_] := Switch[op[[1]],
   "H", pH[st, op[[2]]], "S", pS[st, op[[2]]], "Sdg", pSdg[st, op[[2]]],
   "X", pX[st, op[[2]]], "Y", pY[st, op[[2]]], "Z", pZ[st, op[[2]]],
   "SWAP", pSWAP[st, op[[2]], op[[3]]], "CNOT", pCNOT[st, op[[2]], op[[3]]]];
randOp[n_] := With[{g = RandomChoice[{"H", "S", "Sdg", "X", "Y", "Z", "SWAP", "CNOT"}]},
   If[MatchQ[g, "CNOT" | "SWAP"],
     With[{a = RandomInteger[{1, n}]}, {g, a, With[{b = RandomInteger[{1, n}]}, If[b == a, Mod[a, n] + 1, b]]}],
     {g, RandomInteger[{1, n}]}]];

results = Table[
   Module[{n = nq, ops, ps0, psCanon, stPacked, psPacked, ok},
     ops = Table[randOp[n], {300}];
     ps0 = PauliStabilizer[n];
     psCanon = Fold[canon, ps0, ops];
     stPacked = Fold[pk, getPacked[ps0], ops];
     psPacked = fromPacked[stPacked];
     ok = (psCanon["Tableau"] === psPacked["Tableau"]) && (psCanon["Signs"] === psPacked["Signs"]);
     {n, ok}],
   {nq, {3, 5, 8, 16, 40, 100}}];
pr["equivalence_per_n", results];
pr["all_equiv", AllTrue[results, Last]];

(* round-trip pack/unpack sanity *)
Module[{ps = Fold[canon, PauliStabilizer[20], Table[randOp[20], {100}]]},
  pr["roundtrip_pack", getPacked[ps] // (fromPacked[#]["Tableau"] === ps["Tableau"] &)]];
Print["DONE"];
