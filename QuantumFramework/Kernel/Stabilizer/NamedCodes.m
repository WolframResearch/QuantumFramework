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

(* Patch P1 (2026-05-07, Formula_Test/findings-report.md F1): replace the n-th *)
(* generator from logical X (X-bar) to logical Z (Z-bar) so the named-code     *)
(* states match the api.md docstring "|0_L>" / "|1_L>" instead of returning    *)
(* "|+_L>" / "|-_L>". The 5Q / Steane / 9Q codespaces are stabilized by their *)
(* (n-k) real parity stabilizers; |0_L> is the +1 eigenstate of Z-bar within   *)
(* that codespace (Got97 §3.5, Got00 §4, Nielsen-Chuang §10.5).                *)

PauliStabilizer["5QubitCode"]  := PauliStabilizer[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ",  "ZZZZZ"}]
PauliStabilizer["5QubitCode1"] := PauliStabilizer[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "-ZZZZZ"}]

PauliStabilizer["7QubitCode" | "SteaneCode"]   := PauliStabilizer[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ",  "ZZZZZZZ"}]
PauliStabilizer["7QubitCode1" | "SteaneCode1"] := PauliStabilizer[{"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ", "-ZZZZZZZ"}]

PauliStabilizer["9QubitCode"]  := PauliStabilizer[{"ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX",  "ZZZZZZZZZ"}]
PauliStabilizer["9QubitCode1"] := PauliStabilizer[{"ZZIIIIIII", "IZZIIIIII", "IIIZZIIII", "IIIIZZIII", "IIIIIIZZI", "IIIIIIIZZ", "XXXXXXIII", "IIIXXXXXX", "-ZZZZZZZZZ"}]
