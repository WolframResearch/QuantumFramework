(* Helpers shared by MidCircuitMeasurement.wlt and MeasurementSampling.wlt:
   random Clifford gates, a GHZ preparation, and the density matrix of a
   tableau's state in the computational basis. *)

stabMeasGates = {"H", "S", SuperDagger["S"], "X", "Y", "Z", "CNOT", "CZ", "SWAP"};

stabMeasRandomGate[n_] := With[{g = RandomChoice[stabMeasGates]},
    If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], g -> RandomSample[Range[n], 2], g -> RandomInteger[{1, n}]]
]

stabMeasGHZ[n_] := Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, n - 1}]]

stabMeasTableauDensity[ps_PauliStabilizer] := Normal[ps["State"]["Computational"]["DensityMatrix"]]

(* Exact equality of two exact arrays, entry by entry. *)
stabMeasExactSameQ[a_, b_] := AllTrue[Flatten[{a - b}], PossibleZeroQ]

(* The n-qubit Pauli string with Z on the given wires and I elsewhere. *)
stabMeasZ[n_, wires_List] := StringJoin @ ReplacePart[ConstantArray["I", n], Thread[wires -> "Z"]]
