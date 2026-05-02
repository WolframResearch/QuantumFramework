Package["Wolfram`QuantumFramework`"]

PackageScope[$PauliStabilizerNames]



(* ============================================================================ *)
(* Catalog of named codes.                                                      *)
(* Phase 1 cleanup: dropped duplicate "SteaneCode" entry                        *)
(* (OngoingProjects/Stabilizer/paulistabilizer-source-audit.md \[Section]2.13).           *)
(* ============================================================================ *)

$PauliStabilizerNames = {
    "5QubitCode", "5QubitCode1",
    "SteaneCode", "7QubitCode", "7QubitCode1", "SteaneCode1",
    "9QubitCode", "9QubitCode1",
    "Random"
}


(* ============================================================================ *)
(* Named codes: 5-qubit, Steane (= 7-qubit), Shor (= 9-qubit) and their         *)
(* logical |1_L> variants.                                                      *)
(* References: Got97 \[Section]3.5 (5-qubit cyclic code), Got00 \[Section]4    *)
(* (Steane/Shor CSS construction).                                              *)
(* ============================================================================ *)

PauliStabilizer["5QubitCode"]  := PauliStabilizer[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"}]
PauliStabilizer["5QubitCode1"] := PauliStabilizer[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "-XXXXX"}]

PauliStabilizer["7QubitCode" | "SteaneCode"]   := PauliStabilizer[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ", "XXXXXXX"}]
PauliStabilizer["7QubitCode1" | "SteaneCode1"] := PauliStabilizer[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ", "-XXXXXXX"}]

PauliStabilizer["9QubitCode"]  := PauliStabilizer[{"ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX", "XXXXXXXXX"}]
PauliStabilizer["9QubitCode1"] := PauliStabilizer[{"ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX", "-XXXXXXXXX"}]
