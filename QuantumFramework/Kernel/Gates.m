Package["Wolfram`QuantumFramework`Gates`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport[QuantumGateQ]

PackageExport[XGate]
PackageExport[YGate]
PackageExport[ZGate]
PackageExport[HGate]
PackageExport[SGate]
PackageExport[TGate]
PackageExport[CXGate]
PackageExport[CYGate]
PackageExport[CZGate]
PackageExport[SWAPGate]
PackageExport[CCXGate]
PackageExport[CSWAPGate]
PackageExport[PhaseGate]
PackageExport[SXGate]
PackageExport[SXdgGate]
PackageExport[IdentityGate]
PackageExport[RXGate]
PackageExport[RYGate]
PackageExport[RZGate]
PackageExport[RXXGate]
PackageExport[RYYGate]
PackageExport[RZZGate]
PackageExport[RZXGate]
PackageExport[U1Gate]
PackageExport[U2Gate]
PackageExport[U3Gate]
PackageExport[UGate]



XGate::usage = "XGate[] — Pauli-X (NOT) gate: flips |0⟩ and |1⟩.";
YGate::usage = "YGate[] — Pauli-Y gate: bit and phase flip.";
ZGate::usage = "ZGate[] — Pauli-Z gate: phase flip.";
HGate::usage = "HGate[] — Hadamard gate: creates superposition.";
SGate::usage = "SGate[] — S gate: phase gate (π/2 rotation around Z).";
TGate::usage = "TGate[] — T gate: π/4 phase gate.";
PhaseGate::usage = "PhaseGate[phi] — Phase shift by angle phi.";
SXGate::usage = "SXGate[] — Square root of X gate.";
SXdgGate::usage = "SXdgGate[] — Adjoint of sqrt(X) gate.";
IdentityGate::usage = "IdentityGate[] — identity gate.";

CXGate::usage = "CXGate[] — Controlled-X (CNOT) gate.";
CYGate::usage = "CYGate[] — Controlled-Y gate.";
CZGate::usage = "CZGate[] — Controlled-Z gate.";
SWAPGate::usage = "SWAPGate[] — Swaps two qubits.";
CCXGate::usage = "CCXGate[] — Toffoli (CCX) gate: controlled-controlled-X.";
CSWAPGate::usage = "CSWAPGate[] — Controlled-SWAP (Fredkin) gate.";

RXGate::usage = "RXGate[theta] — Rotation around X axis by angle theta.";
RYGate::usage = "RYGate[theta] — Rotation around Y axis by angle theta.";
RZGate::usage = "RZGate[theta] — Rotation around Z axis by angle theta.";
RXXGate::usage = "RXXGate[theta] — Two-qubit XX rotation by angle theta.";
RYYGate::usage = "RYYGate[theta] — Two-qubit YY rotation by angle theta.";
RZZGate::usage = "RZZGate[theta] — Two-qubit ZZ rotation by angle theta.";
RZXGate::usage = "RZXGate[theta] — Two-qubit ZX rotation by angle theta.";

U1Gate::usage = "U1Gate[lambda] — U1 single-qubit phase gate.";
U2Gate::usage = "U2Gate[phi, lambda] — U2 single-qubit gate.";
U3Gate::usage = "U3Gate[theta, phi, lambda] — U3 general single-qubit gate.";
UGate::usage = "UGate[theta, phi, lambda] — General single-qubit gate (U3).";


QuantumGateQ[_XGate | _YGate | _ZGate | _HGate | _SGate | _TGate | _CXGate | _CYGate | _CZGate |
  _SWAPGate | _CCXGate | _CSWAPGate | _PhaseGate | _SXGate | _SXdgGate | _IdentityGate |
  _RXGate | _RYGate | _RZGate | _RXXGate | _RYYGate | _RZZGate | _RZXGate |
  _U1Gate | __U2Gate | _U3Gate | _UGate
] := True
QuantumGateQ[_] := False


(* All gate QuantumOperator rules delegate to canonical symbolic names, supporting all usage patterns. Only internal gate parameters (like angles) are included in the head, never order or options. *)
QuantumOperator[XGate[d_ : 2], args___] := QuantumOperator["X"[d], args]
QuantumOperator[YGate[d_ : 2], args___] := QuantumOperator["Y"[d], args]
QuantumOperator[ZGate[d_ : 2], args___] := QuantumOperator["Z"[d], args]
QuantumOperator[HGate[d_ : 2], args___] := QuantumOperator["H"[d], args]
QuantumOperator[SGate[], args___] := QuantumOperator["S", args]
QuantumOperator[TGate[], args___] := QuantumOperator["T", args]
QuantumOperator[PhaseGate[phi_ : Pi], args___] := QuantumOperator[{"Phase", phi}, args]
QuantumOperator[SXGate[], args___] := QuantumOperator["SX", args]
QuantumOperator[SXdgGate[], args___] := QuantumOperator["SX", args]["Dagger"]
QuantumOperator[IdentityGate[d_ : 2], args___] := QuantumOperator["I"[d], args]

QuantumOperator[CXGate[], args___] := QuantumOperator["CX", args]
QuantumOperator[CYGate[], args___] := QuantumOperator["CY", args]
QuantumOperator[CZGate[], args___] := QuantumOperator["CZ", args]
QuantumOperator[SWAPGate[], args___] := QuantumOperator["SWAP", args]
QuantumOperator[CCXGate[], args___] := QuantumOperator["C"["CX"], args]
QuantumOperator[CSWAPGate[], args___] := QuantumOperator["CSWAP", args]

QuantumOperator[RXGate[theta_ : 0], args___] := QuantumOperator["RX"[theta], args]
QuantumOperator[RYGate[theta_ : 0], args___] := QuantumOperator["RY"[theta], args]
QuantumOperator[RZGate[theta_ : 0], args___] := QuantumOperator["RZ"[theta], args]
QuantumOperator[RXXGate[theta_ : 0], args___] := QuantumOperator["R"[theta, "XX"], args]
QuantumOperator[RYYGate[theta_ : 0], args___] := QuantumOperator["R"[theta, "YY"], args]
QuantumOperator[RZZGate[theta_ : 0], args___] := QuantumOperator["R"[theta, "ZZ"], args]
QuantumOperator[RZXGate[theta_ : 0], args___] := QuantumOperator["R"[theta, "YX"], args]

QuantumOperator[U1Gate[lambda_ : Pi], args___] := QuantumOperator["U1"[lambda], args]
QuantumOperator[U2Gate[phi_ : 0, lambda_ : Pi], args___] := QuantumOperator["U2"[phi, lambda], args]
QuantumOperator[(UGate | U3Gate)[theta_, phi_, lambda_], args___] := QuantumOperator["U"[theta, phi, lambda], args]

