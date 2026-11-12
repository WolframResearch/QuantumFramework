Package["Wolfram`QuantumFramework`"]

PackageExport["QiskitCircuit"]
PackageExport["ImportQASMCircuit"]
PackageExport["ImportQPY"]

PackageScope["QuantumCircuitOperatorToQiskit"]
PackageScope["qiskitQASM"]
PackageScope["qiskitPrimitiveSubmit"]
PackageScope["qiskitReorderByQubit"]
PackageScope["$qiskitVirtualCircuitPython"]


(* The qiskit name of a gate shortcut (QuantumShortcut's output for a unitary operator). The
   circuit reaches this only through circuitInstructions, which has already turned measurements,
   resets, injections and permutations into gates, resets and computational measurements. *)
shortcutToGate = Replace[
    {
        "R"[angle_, {name_ -> order_ ? orderQ}] :> shortcutToGate["R"[angle, name]] -> order,
        (name_ -> (order_ ? orderQ)) :> With[{shortcut = shortcutToGate[name]}, If[shortcut === Nothing, Nothing, shortcut -> order]],
        "I" -> Nothing,
        "X" | "NOT" -> "XGate",
        "Y" -> "YGate",
        "Z" | "1" -> "ZGate",
        "H" -> "HGate",
        "S" -> "SGate",
        "T" -> "TGate",
        "V" -> "SXGate",
        "SWAP" -> "SwapGate",
        "Diagonal"[x_] :> If[ListQ[x], {"Diagonal", NumericArray @ N[x]}, {"GlobalPhaseGate", Chop[N[Arg[x]]]}],
        "U2"[a_, b_] :> {"U2Gate", N[a], N[b]},
        "U"[a_, b_, c_] :> {"U3Gate", N[a], N[b], N[c]},
        "GlobalPhase"[phase_] :> {"GlobalPhaseGate", Chop[Re[N[phase]]]},
        "R"[angle_, "X"] :> {"RXGate", N[angle]},
        "R"[angle_, "Y"] :> {"RYGate", N[angle]},
        "R"[angle_, "Z"] :> {"RZGate", N[angle]},
        "P"[phase_] :> {"PhaseGate", N[phase]},
        "PhaseShift"[k_] :> {"PhaseGate", N[Sign[k] 2 Pi / 2 ^ Abs[k]]},
        "C"[name_, controls___] :> {"Control", shortcutToGate[name], controls},
        SuperDagger[name_] :> {"Dagger", shortcutToGate[name]},
        "Delay"[delay_] :> {"Delay", delay},
        (* a gate given by its matrix: the matrix and its wires both from the sorted operator *)
        Labeled[arr_ /; ArrayQ[arr] || NumericArrayQ[arr] -> order_, label_] :>
            {"Unitary", NumericArray[Normal @ N @ arr, "ComplexReal64"], ToString[label], order},
        shortcut_ :> With[{op = QuantumOperator[shortcut]["Sort"]},
            {"Unitary", NumericArray[Normal @ N @ op["Matrix"], "ComplexReal64"], ToString[shortcut], op["Order"]}
        ]
    }
]

qiskitGates[{"Gate", qo_}] := shortcutToGate /@ QuantumShortcut[qo]
qiskitGates[{"Reset", w_}] := {"Reset" -> {w}}
qiskitGates[{"Measure", w_, k_}] := {{"Measure", w, k}}
qiskitGates[{"Barrier", wires_}] := {"Barrier" -> wires}

QuantumCircuitOperatorToQiskit[qco_QuantumCircuitOperator] := Enclose @ Block[{
    instructions = Confirm @ circuitInstructions[qco], gates, arity
},
    gates = Catenate[qiskitGates /@ instructions];
    ConfirmAssert[FreeQ[gates, _Failure]];
    arity = {qco["Max"], Count[instructions, {"Measure", __}]};
    Confirm @ PythonEvaluate[Context[arity], "
from wolframclient.language import wl
from wolframclient.language.expression import WLFunction

from qiskit import QuantumCircuit
from qiskit.circuit.gate import Gate
from qiskit.circuit import Delay, Measure
from qiskit.circuit.library import *

import pickle

circuit = QuantumCircuit(
    * <* arity *>
)

def make_gate(gate_spec):
    order = []
    if type(gate_spec) == WLFunction and gate_spec.head == wl.Rule:
        gate, order = make_gate(gate_spec[0])
        order = list(gate_spec[1])
        return gate, order
    elif isinstance(gate_spec, tuple):
        name, *args = gate_spec
    else:
        name = gate_spec
        args = []
    if name == 'Control':
        base_gate, base_order = make_gate(args[0])
        size1 = len(args[1])
        size0 = len(args[2])
        gate = base_gate.control(size1 + size0, ctrl_state='0' * size0 + '1' * size1)
        order = list(args[1]) + list(args[2])
        if base_gate.name != 'global_phase':
            order = order + base_order
    elif name == 'Dagger':
        base_gate = make_gate(args[0])[0]
        gate = base_gate.inverse()
    elif name == 'Unitary':
        gate = UnitaryGate(args[0], label=args[1], check_input=False)
        assert(all(i == j for i, j in zip(args[2][0], args[2][1])))
        # QF writes the matrix with its first wire most significant, qiskit with its first qubit least
        order = list(args[2][0])[::-1]
    elif name == 'Barrier':
        # a bare barrier (no explicit qubits) spans the whole circuit
        if len(order) == 0:
            order = list(range(1, circuit.num_qubits + 1))
        gate = Barrier(len(order))
    elif name == 'Delay':
        gate = Delay(args[0])
    else:
        gate = getattr(qiskit.circuit.library, name)(*args)
    return gate, order

def add_gate(circuit, gate_spec):
    if isinstance(gate_spec, tuple) and gate_spec[0] == 'Measure':
        circuit.append(Measure(), [gate_spec[1] - 1], [gate_spec[2]])
        return
    gate, order = make_gate(gate_spec)
    order = [o - 1 for o in order]
    # an uncontrolled global phase belongs to the circuit, not to an instruction (OpenQASM 2 has
    # no zero-qubit gate, and qiskit cannot read back its own OpenQASM 3 for one)
    if gate.name == 'global_phase':
        circuit.global_phase += gate.params[0]
        return
    if gate.name == 'barrier':
        if len(order) == 0:
            order = list(range(circuit.num_qubits))
        gate = Barrier(len(order))
    if gate.name == 'reset':
        for o in order:
            circuit.append(gate, (o, ))
    else:
        circuit.append(gate, tuple(order))

for gate_spec in <* gates *>:
    add_gate(circuit, gate_spec)

wl.Wolfram.QuantumFramework.QiskitCircuit(pickle.dumps(circuit))
    "]
]

QiskitCircuit[qc_QuantumCircuitOperator] := qc["Qiskit"]

QiskitCircuit[bytes_ByteArray]["Bytes"] := bytes

QiskitCircuit[s_String /; qasmStringQ[s]] := QuantumCircuitOperator[s]["Qiskit"]

QiskitCircuit::badarg = "QiskitCircuit wraps a pickled ByteArray. Build one from a circuit with QuantumCircuitOperator[\[Ellipsis]][\"Qiskit\"]; a plain string is read only when it is OpenQASM source."

QiskitCircuit[_String] := (Message[QiskitCircuit::badarg]; $Failed)

qc_QiskitCircuit["Eval", attr_String, args___, kwargs : OptionsPattern[]] :=
    PythonEvalAttr[{qc["Bytes"], attr, args, kwargs}]

qc_QiskitCircuit["EvalBytes", attr_String, args___, kwargs : OptionsPattern[]] :=
    PythonEvalAttr[{qc["Bytes"], attr, args, kwargs}, "ReturnBytes" -> True]


qiskitGraphics[qc_QiskitCircuit, OptionsPattern[{"Scale" -> 5}]] := Enclose @ With[{
    latex = (
        PythonEvaluate["
import os
if 'texbin' not in os.environ['PATH']:
    os.environ['PATH'] += os.pathsep + os.pathsep.join(['/Library/TeX/texbin'])
"];
        qc["Eval", "draw", "output" -> "latex_source", "scale" -> OptionValue["Scale"]]
    )
},
    Confirm @ Check[Needs["MaTeX`"], ResourceFunction["MaTeXInstall"][]; Needs["MaTeX`"]];
    MaTeX`MaTeX[
        First @ StringCases[latex, "\\begin{document}" ~~ s___ ~~ "\\end{document}" :> s],
        "Preamble" -> Prepend[StringCases[latex, s : ("\\usepackage" ~~ Shortest[___] ~~ EndOfLine) :> s], "\\standaloneconfig{border=-16 -8 -16 -16}"]
    ]
]

Options[qiskitDiagram] = {"Scale" -> 5}

qiskitDiagram[qc_QiskitCircuit, opts : OptionsPattern[]] := Enclose @ Block[{
    $pythonBytes = qc["Bytes"],
    scale = OptionValue["Scale"],
    image
},
    image = Confirm @ PythonEvaluate[Context[scale], "
import pickle
import PIL
import matplotlib
matplotlib.use('Agg')
fig = pickle.loads(<* $pythonBytes *>).draw(output='mpl', style='iqp', scale=<* scale *>)
fig.canvas.draw()
image = PIL.Image.frombytes('RGBA', fig.canvas.get_width_height(), fig.canvas.buffer_rgba())
matplotlib.pyplot.close()
image
"];
    ImageCrop[ConfirmBy[image, ImageQ]]
]


qc_QiskitCircuit["Graphics", opts : OptionsPattern[]] := qiskitGraphics[qc, opts]

qc_QiskitCircuit["Diagram", opts : OptionsPattern[]] := qiskitDiagram[qc, opts]

qc_QiskitCircuit["Qubits"] := qc["Eval", "num_qubits"]

qc_QiskitCircuit["Depth"] := qc["Eval", "depth"]

qc_QiskitCircuit["Ops"] := qc["Eval", "count_ops"]

qc_QiskitCircuit["Clbits"] := Block[{
    $pythonBytes = qc["Bytes"]
},
    PythonEvaluate[Context[$pythonBytes], "
import pickle
qc = pickle.loads(<* $pythonBytes *>)
[qc.find_bit(clbit).index for clbit in qc.clbits]
"]
]


(* A transpiled circuit acts on physical qubits and records, in its layout, where each virtual
   qubit started (initial layout) and where it ended (final permutation). qf_virtual_circuit
   rewrites it on the virtual qubits: qubit v of the result is qubit v of the circuit before
   transpiling, and SWAP gates appended at the end bring each qubit back to its own wire. A
   physical qubit the layout adds (an ancilla) is kept only when a gate uses it. With basis_gates
   the SWAPs are written in those gates, without optimization (which would remove them again).
   Every reader of a transpiled circuit goes through this one function: the gate-set QASM export,
   the decoder to QF and the "Matrix" property. *)
$qiskitVirtualCircuitPython = "
def qf_virtual_circuit(circ, basis_gates=None):
    from qiskit import QuantumCircuit, QuantumRegister, transpile
    lay = getattr(circ, 'layout', None)
    if lay is None:
        return circ
    init = list(lay.initial_index_layout(filter_ancillas=False))
    n = len(list(lay.initial_index_layout(filter_ancillas=True)))
    used = set()
    for inst in circ.data:
        used.update(circ.find_bit(q).index for q in inst.qubits)
    phys_to_virt = {}
    size = n
    for v, p in enumerate(init):
        if v < n:
            phys_to_virt[p] = v
    for v, p in enumerate(init):
        if v >= n and p in used:
            phys_to_virt[p] = size
            size += 1
    out = QuantumCircuit(QuantumRegister(size, 'q'), *circ.cregs)
    out.global_phase = circ.global_phase
    for inst in circ.data:
        out.append(inst.operation, [phys_to_virt[circ.find_bit(q).index] for q in inst.qubits], [circ.find_bit(c).index for c in inst.clbits])
    final = list(lay.final_index_layout(filter_ancillas=True))
    # where[v]: the wire now holding virtual qubit v; held[w]: the virtual qubit on wire w
    where = [phys_to_virt[p] for p in final]
    if where != list(range(n)):
        held = [None] * size
        for v, w in enumerate(where):
            held[w] = v
        back = QuantumCircuit(size)
        for v in range(n):
            w = where[v]
            if w != v:
                back.swap(v, w)
                u = held[v]
                held[w], held[v] = u, v
                where[v] = v
                if u is not None:
                    where[u] = w
        if basis_gates is not None:
            back = transpile(back, basis_gates=basis_gates, optimization_level=0)
        out.compose(back, qubits=list(range(size)), inplace=True)
    return out
";


(* The decoder from qiskit to QF. Python only describes each instruction (its name, parameters,
   qubits and classical bits, and for a gate outside the gate table its controls, definition or
   matrix); the QF circuit is built here, from the gate table the OpenQASM importer uses. *)
$qiskitDescribePython = "
import numpy as np
from qiskit.circuit import ParameterExpression, ControlledGate

def qf_param(x):
    if isinstance(x, ParameterExpression):
        try:
            return float(x)
        except Exception:
            return wl.ToExpression(str(x))
    if isinstance(x, (int, float, complex)):
        return x
    return None

def qf_describe(op, qs, cs, known):
    d = {'name': op.name, 'qubits': qs, 'clbits': cs}
    if getattr(op, 'condition', None) is not None:
        d['conditioned'] = True
    params = [qf_param(x) for x in getattr(op, 'params', [])]
    if all(p is not None for p in params):
        d['params'] = params
    if op.name in known or op.name in ('measure', 'reset', 'barrier', 'delay', 'global_phase') or op.name.startswith('save_'):
        return d
    if isinstance(op, ControlledGate) and op.base_gate is not None:
        k = op.num_ctrl_qubits
        d['ctrl_bits'] = [(op.ctrl_state >> i) & 1 for i in range(k)]
        d['base'] = qf_describe(op.base_gate, qs[k:], [], known)
        return d
    if op.name != 'unitary' and getattr(op, 'definition', None) is not None:
        sub = op.definition
        d['sub'] = [qf_describe(i.operation, [qs[sub.find_bit(q).index] for q in i.qubits], [cs[sub.find_bit(c).index] for c in i.clbits], known) for i in sub.data]
        d['sub_phase'] = qf_param(sub.global_phase)
        return d
    try:
        m = np.asarray(op.to_matrix(), dtype=complex)
        d['re'] = m.real.tolist()
        d['im'] = m.imag.tolist()
    except Exception:
        pass
    return d

def qf_describe_circuit(circ, known):
    return {
        'num_qubits': circ.num_qubits,
        'global_phase': qf_param(circ.global_phase),
        'instructions': [qf_describe(i.operation, [circ.find_bit(q).index + 1 for q in i.qubits], [circ.find_bit(c).index for c in i.clbits], known) for i in circ.data]
    }
";

(* an angle qiskit holds as a float is read back as a rational multiple of Pi when it is one to
   machine precision (pi/2, 3 pi/4, ...), otherwise kept as the float *)
qiskitParameter[x_Real] := With[{r = Rationalize[x / Pi, 10 ^ -12]}, If[Denominator[r] <= 64, r Pi, x]]
qiskitParameter[x_] := x

qiskitElements[d_Association] := Enclose @ With[{
    name = d["name"], qs = d["qubits"], params = qiskitParameter /@ Lookup[d, "params", {}]
},
    Which[
        KeyExistsQ[d, "conditioned"],
            Confirm @ Failure["Qiskit", <|"MessageTemplate" -> "The classically controlled operation `1` has no QF circuit element.", "MessageParameters" -> {name}|>],
        name === "measure", {QuantumMeasurementOperator[qs]},
        name === "reset", {"Reset" -> qs},
        name === "barrier", {"Barrier"[qs]},
        name === "delay" || StringStartsQ[name, "save_"], {},
        name === "global_phase", {QuantumOperator["GlobalPhase"[First[params]]]},
        MemberQ[qasmGateNames[], name] && KeyExistsQ[d, "params"],
            ConfirmBy[qasmGateElements[name, params, qs], ListQ],
        KeyExistsQ[d, "base"],
            With[{
                base = QuantumOperator @ QuantumCircuitOperator[Confirm @ qiskitElements[d["base"]]],
                ctrls = qs[[;; Length[d["ctrl_bits"]]]]
            },
                {QuantumOperator["C"[base, Pick[ctrls, d["ctrl_bits"], 1], Pick[ctrls, d["ctrl_bits"], 0]]]}
            ],
        KeyExistsQ[d, "sub"],
            {QuantumCircuitOperator[
                Join[
                    Catenate[Confirm @ qiskitElements[#] & /@ d["sub"]],
                    If[TrueQ[Chop[N[Lookup[d, "sub_phase", 0]]] == 0], {}, {QuantumOperator["GlobalPhase"[qiskitParameter @ d["sub_phase"]]]}]
                ],
                name
            ]},
        KeyExistsQ[d, "re"],
            (* qiskit writes a gate's matrix with its first qubit least significant *)
            {QuantumOperator[d["re"] + I d["im"], Reverse[qs], "Label" -> name]},
        True,
            Confirm @ Failure["Qiskit", <|"MessageTemplate" -> "The qiskit operation `1` has no QF circuit element.", "MessageParameters" -> {name}|>]
    ]
]

qc_QiskitCircuit["QuantumCircuit" | "QuantumCircuitOperator"] := Enclose @ Block[{
    $pythonBytes = qc["Bytes"], $known = qasmGateNames[], desc
},
    desc = ConfirmBy[PythonEvaluate[Context[$pythonBytes], $qiskitVirtualCircuitPython <> $qiskitDescribePython <> "
import pickle
from wolframclient.language import wl
qf_describe_circuit(qf_virtual_circuit(pickle.loads(<* $pythonBytes *>)), set(<* $known *>))
"], AssociationQ];
    QuantumCircuitOperator[
        Join[
            Catenate[Confirm @ qiskitElements[#] & /@ desc["instructions"]],
            If[TrueQ[Chop[N[desc["global_phase"]]] == 0], {}, {QuantumOperator["GlobalPhase"[qiskitParameter @ desc["global_phase"]]]}]
        ]
    ]
]

QiskitCircuit /: QuantumCircuitOperator[qc_QiskitCircuit] := qc["QuantumCircuitOperator"]


(* the circuit's operator on its virtual qubits (final measurements removed) *)
qiskitMatrix[qc_QiskitCircuit] := Block[{$pythonBytes = qc["Bytes"]},
    PythonEvaluate[Context[$pythonBytes], $qiskitVirtualCircuitPython <> "
from qiskit.quantum_info import Operator

import pickle

qc = qf_virtual_circuit(pickle.loads(<* $pythonBytes *>))
qc.remove_final_measurements()
Operator(qc).data
"]
]

$IBMProvider = "IBM" | "IBMQ" | "IBMProvider" | "IBMRuntime"
$AWSProvider = "AWS" | "AWSBraket" | "Braket"

Options[qiskitInitBackend] = {"Provider" -> None, "Backend" -> Automatic, "FireOpal" -> False, "DensityMatrix" -> False}

qiskitInitBackend[qc_QiskitCircuit, OptionsPattern[]] := Enclose @ Block[{
    $pythonBytes = qc["Bytes"],
    provider, params,
    $backendName = Replace[OptionValue["Backend"], Automatic -> Null],
    $fireOpal = TrueQ[OptionValue["FireOpal"]],
    $token = Null,
    $matrixQ = TrueQ[OptionValue["DensityMatrix"]],
    env
},
    {provider, params} = Replace[OptionValue["Provider"], {
        {name_, params : OptionsPattern[]} | name_ :> {name, Flatten[{params}]}
    }];
    If[provider === None && $fireOpal, provider = "IBMRuntime"];

    Switch[
        provider
        ,
        $IBMProvider,
        $token = Lookup[
            params,
            "Token",
            Enclose[ConfirmBy[SystemCredential["IBMQuantumPlatform_APIKEY"], StringQ], Null &]
        ];
        env = "ibm"
        ,
        $AWSProvider,
        env = "braket",
        _,
        env = "default"
    ];

    If[$fireOpal, env = "qctrl"];

    Confirm @ PythonEvaluate[Context[$pythonBytes], Switch[provider,
    $IBMProvider,
"
from qiskit_ibm_runtime import QiskitRuntimeService
try:
    assert(isinstance(provider, QiskitRuntimeService))
except:
    token = <* $token *>
    if token is not None:
        QiskitRuntimeService.save_account(channel='ibm_quantum_platform', token=token, overwrite=True)
    provider = QiskitRuntimeService()
",
    $AWSProvider,
"
from qiskit_braket_provider import AWSBraketProvider
try:
    assert(isinstance(provider, AWSBraketProvider))
except:
    provider = AWSBraketProvider()
",
    _,
"
from qiskit.providers.basic_provider import BasicProvider
provider = BasicProvider()
"], env];
    Confirm @ PythonEvaluate[Context[$pythonBytes],
"
import pickle

qc = pickle.loads(<* $pythonBytes *>)
backend_name = <* $backendName *>
", env];
    Confirm @ PythonEvaluate[Context[$pythonBytes], Switch[provider,
    $IBMProvider,
"
if backend_name is None:
    backend = provider.backends()[0]
else:
    backend = provider.backend(backend_name)
if <* $fireOpal *>:
    import fireopal
    from fireopal.credentials import make_credentials_for_ibmq, make_credentials_for_braket
    if isinstance(provider, QiskitRuntimeService):
        fireopal_credentials = make_credentials_for_ibmq(provider._account.token, 'open', 'ibm-q', 'main')
    else:
        raise ValueError(f'Unsupported FireOpal provider: {provider}')
",
    $AWSProvider,
"
if backend_name is None:
    backend = provider.get_backend('SV1')
else:
    backend = provider.get_backend('basic_simulator')
",
    (* A representative fake device with a populated error-aware Target (per-instruction
       error/duration, readout error, coupling map). No credentials and no network: used to
       exercise the error-aware transpile / target-extraction paths offline. Backend is the
       qubit count. *)
    "GenericV2",
"
from qiskit.providers.fake_provider import GenericBackendV2
backend = GenericBackendV2(num_qubits=int(backend_name) if backend_name is not None else 5, seed=42)
",
    _,
"
from qiskit_aer import AerSimulator
# the simulation method only; qiskitApply adds the save instruction to its own copy of the circuit
if qc.num_clbits == 0:
    method = 'density_matrix' if <* $matrixQ *> else 'statevector'
else:
    method = None

backend = AerSimulator(method=method)
" 
], env];
    env
]

(* qiskit get_counts gives string keys "c_{n-1}...c_0"; reverse to the WL bit-list outcome *)
qiskitBitCounts[result_Association] := KeyMap[
    Reverse[Characters[StringDelete[#, Whitespace]]] /. {"0" -> 0, "1" -> 1} &
] @ result

(* IBM / qiskit return each shot as a classical-register bit list in classical-bit order,
   which follows the order the measurements were emitted to qiskit (the measurement's target
   order). QuantumMeasurement instead reports outcomes in ascending qubit order, so a circuit
   whose measurement targets are not already ascending (a permuted or interleaved measurement,
   or any backend whose transpiled measure -> clbit map is not the identity) needs the bits
   reordered. measuredQubits[[c+1]] is the original (pre-transpile) qubit index read into
   classical bit c, captured at transpile time; sorting the {clbit -> bit} pairs by that qubit
   index restores the ascending-qubit order. A missing (Automatic) or length-mismatched map is
   a no-op, leaving the classical-bit order untouched. *)
qiskitReorderByQubit[bits_List, measuredQubits_List] /; Length[measuredQubits] === Length[bits] :=
    Values @ KeySort @ AssociationThread[measuredQubits -> bits]
qiskitReorderByQubit[bits_List, _] := bits

qiskitMeasurementFromBitCounts[bitCounts_Association] := Block[{size = Length @ First @ Keys @ bitCounts},
    Enclose[
        ConfirmAssert[size <= 8];
        ConfirmBy[
            QuantumMeasurement[
                Join[Association[# -> 0 & /@ IntegerDigits[Range[2 ^ size] - 1, 2, size]], bitCounts]
            ],
            QuantumMeasurementQ
        ],
        bitCounts &
    ]
]

qiskitCountsToMeasurement[result_Association] := qiskitMeasurementFromBitCounts @ qiskitBitCounts @ result

Options[qiskitApply] = Join[{"Shots" -> 1024, "Validate" -> False}, Options[qiskitInitBackend]]

qiskitApply[qc_QiskitCircuit, qs: _ ? QuantumStateQ | Automatic, opts : OptionsPattern[]] := Enclose @ Block[{
    $state = Replace[qs, {Automatic -> Null, _ :> If[qs["Dimension"] == 1, Null, NumericArray[Normal @ N @ qs["Reverse"]["StateVector"], "ComplexReal64"]]}],
    $shots = OptionValue["Shots"],
    $fireOpal = TrueQ[OptionValue["FireOpal"]],
    $validate = TrueQ[OptionValue["Validate"]],
    $matrixQ = TrueQ[OptionValue["DensityMatrix"]],
    env, result
},
    If[ QuantumStateQ[qs],
        ConfirmAssert[qs["InputDimensions"] == {}];
        ConfirmAssert[AllTrue[qs["OutputDimensions"], EqualTo[2]]]
    ];

    env = Confirm @ qiskitInitBackend[qc, FilterRules[{opts}, Options[qiskitInitBackend]]];

    result = PythonEvaluate[Context[$state], "
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

if isinstance(backend, AerSimulator):
    clbits = qc.num_clbits
else:
    clbits = qc.num_clbits if qc.num_clbits > 0 else qc.num_qubits

circuit = QuantumCircuit(qc.num_qubits, clbits)

state = <* $state *>
if state is not None:
    circuit.initialize(state)
circuit = circuit.compose(qc)

if isinstance(backend, AerSimulator) and circuit.num_clbits == 0:
    if <* $matrixQ *>:
        circuit.save_density_matrix()
    else:
        circuit.save_statevector()

if not isinstance(backend, AerSimulator) and clbits == 0:
    circuit.measure_all()

try:
    from qiskit.transpiler import generate_preset_pass_manager
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    circuit = pm.run(circuit)
except:
    pass
if <* $fireOpal *>:
    from qiskit import qasm2
    qasm = qasm2.dumps(circuit)
    if <* $validate *>:
        validate_results = fireopal.validate(
            circuits=[qasm], credentials=fireopal_credentials
        )
        assert validate_results['results'] == [], validate_results['results'][0]['error_message']
    result = fireopal.execute(
        circuits=[qasm],
        shot_count=<* $shots *>,
        credentials=fireopal_credentials
    )['results'][0]
    result = {k: v for k, v in result.items()}
elif isinstance(backend, AerSimulator) and circuit.num_clbits == 0:
    import numpy as np
    result = backend.run(circuit).result()
    if <* $matrixQ *>:
        result = result.data()['density_matrix']
    else:
        result = result.get_statevector()
    result = np.array(result)
else:
    from qiskit_ibm_runtime import SamplerV2 as Sampler
    sampler = Sampler(mode=backend)
    _job = sampler.run([circuit], shots = <* $shots *>)
    try:
        _jid = _job.job_id()
    except Exception:
        _jid = None
    # classical bit -> original (pre-transpile) qubit, so the counts can be reordered into
    # ascending-qubit order regardless of the transpiled measure -> clbit layout
    _c2p = {}
    for _ins in circuit.data:
        if _ins.operation.name == 'measure':
            _c2p[circuit.find_bit(_ins.clbits[0]).index] = circuit.find_bit(_ins.qubits[0]).index
    _ncl = circuit.num_clbits
    _measured = None
    _layout = getattr(circuit, 'layout', None)
    if _layout is not None:
        try:
            _final = _layout.final_index_layout(filter_ancillas=True)
            _p2o = {ph: o for o, ph in enumerate(_final)}
            _measured = [_p2o.get(_c2p.get(c), _c2p.get(c)) for c in range(_ncl)]
        except Exception:
            _measured = None
    if _measured is None:
        _measured = [_c2p.get(c, c) for c in range(_ncl)]
    result = {'__counts__': _job.result()[0].join_data().get_counts(),
              '__job_id__': _jid,
              '__backend__': getattr(backend, 'name', None),
              '__measured_qubits__': _measured}
result
", env];
    Which[
        (* IBM Runtime SamplerV2: return an IBMJob carrying the job id + measurement. The raw
           counts are keyed by classical-bit order; reorder each outcome into ascending-qubit
           order via the captured classical-bit -> original-qubit map before building the
           measurement, so a permuted measurement or non-identity layout decodes correctly. *)
        AssociationQ[result] && KeyExistsQ[result, "__job_id__"],
        With[{bitCounts = KeyMap[
                qiskitReorderByQubit[#, Lookup[result, "__measured_qubits__", Automatic]] &,
                qiskitBitCounts[result["__counts__"]]
            ]},
            With[{qm = qiskitMeasurementFromBitCounts[bitCounts]},
                If[ env === "ibm" && StringQ[result["__job_id__"]],
                    iMethodIBMJob[qm, bitCounts, result["__job_id__"], result["__backend__"]],
                    qm
                ]
            ]
        ],
        AssociationQ[result],
        qiskitCountsToMeasurement[result],
        NumericArrayQ[result],
        QuantumState[Chop @ Normal[result]]["Reverse"],
        True,
        result
    ]
]

qc_QiskitCircuit["Decompose", n : _Integer ? Positive : 1] := QiskitCircuit[qc["EvalBytes", "decompose", "reps" -> n]]

qc_QiskitCircuit["QuantumOperator"] := Enclose @ Module[{mat = Chop @ Normal[ConfirmBy[qiskitMatrix[qc], NumericArrayQ]]},
    ConfirmAssert[SquareMatrixQ[mat]];
    QuantumOperator[mat]["Reverse"]["Sort"]
]

qc_QiskitCircuit["Matrix"] := qc["QuantumOperator"]["Matrix"]

qc_QiskitCircuit[opts : OptionsPattern[qiskitApply]] := qiskitApply[qc, Automatic, opts]

qc_QiskitCircuit[qs_QuantumState, opts : OptionsPattern[qiskitApply]] := qiskitApply[qc, qs, opts]


qiskitQASM[qc_, opts : OptionsPattern[Join[{"Version" -> 2, "OptimizationLevel" -> Automatic}, Options[qiskitInitBackend]]]] := Enclose @ Block[
    {$version = OptionValue["Version"], $optLevel = Replace[OptionValue["OptimizationLevel"], Automatic -> 1], env},
    env = Confirm @ qiskitInitBackend[qc, FilterRules[{opts}, Options[qiskitInitBackend]]];
    (* transpile against the live backend's Target, which carries per-instruction error +
       duration: layout and routing are error-aware. optimization_level is user-controllable. *)
    PythonEvaluate[Context[$version], "
from qiskit import transpile
circuit = transpile(qc, backend, optimization_level=<* $optLevel *>)
try:
    from qiskit_braket_provider import AWSBraketProvider
    if isinstance(provider, AWSBraketProvider):
        from qiskit_braket_provider.providers.adapter import convert_qiskit_to_braket_circuit
        from braket.circuits.serialization import IRType
        result = convert_qiskit_to_braket_circuit(circuit).to_ir(IRType.OPENQASM).source
except:
    if <* $version *> == 3:
        from qiskit import qasm3
        result = qasm3.dumps(circuit)
    else:
        from qiskit import qasm2
        result = qasm2.dumps(circuit)
result
", env]
]

(* QiskitTarget["FromBackend" -> name, "Provider" -> p] (the generic backend extractor): build an
   error-aware QiskitTarget by reading a live backend's own qiskit Target, which is already
   populated with per-instruction error + duration. ONE implementation for every qiskit-ecosystem
   backend (IBM, AWS Braket, and the local GenericV2 fake device) through the shared
   qiskitInitBackend -- the provider is a parameter, not hardcoded. The result is the same
   InstructionProperties-carrying spec the arbel adapter produces, so downstream transpilation is
   error-aware identically. *)
QiskitTarget["FromBackend" -> name_, opts : OptionsPattern[qiskitInitBackend]] := Enclose @ Block[{$session, qc, env, spec},
    qc  = ConfirmMatch[QuantumCircuitOperator[{"H" -> 1}]["Qiskit"], _QiskitCircuit];
    env = Confirm @ qiskitInitBackend[qc, "Backend" -> name,
        Sequence @@ FilterRules[{opts}, {"Provider", "FireOpal", "DensityMatrix"}]];
    spec = Confirm @ PythonEvaluate[Context[$session], "
from wolframclient.language import wl
t = backend.target
# Target.operation_names mixes real gates with control-flow ops (if_else, for_loop, ...) and
# directives (delay, barrier); only standard gate names are valid basis_gates, so filter to
# those. measure is kept for the instruction-properties walk (its error steers readout-aware
# layout) but not put in basis_gates (the carrier adds measure/reset universally).
try:
    from qiskit.circuit.library import get_standard_gate_name_mapping
    gateset = set(get_standard_gate_name_mapping())
except Exception:
    gateset = None
nonbasis = {'measure', 'reset', 'delay', 'barrier', 'if_else', 'for_loop', 'while_loop', 'switch_case'}
def _isgate(nm):
    return (nm in gateset) if gateset is not None else (nm not in nonbasis)
basis = sorted(nm for nm in t.operation_names if _isgate(nm))
nq = t.num_qubits
cm = t.build_coupling_map()
coupling = [list(e) for e in cm.get_edges()] if cm is not None else []
props = []
for nm in [n for n in t.operation_names if _isgate(n) or n == 'measure']:
    qmap = t[nm]
    if not hasattr(qmap, 'items'):
        continue
    for qargs, ip in qmap.items():
        if qargs is None:
            continue
        row = {}
        if ip is not None:
            if getattr(ip, 'error', None) is not None:
                row['Error'] = ip.error
            if getattr(ip, 'duration', None) is not None:
                row['Duration'] = ip.duration
        if row:
            props.append([nm, list(qargs), row])
result = wl.Association(
    wl.Rule('BasisGates', basis),
    wl.Rule('NumQubits', nq),
    wl.Rule('CouplingMap', coupling),
    wl.Rule('InstructionProperties', props))
result
", env];
    ConfirmMatch[QiskitTarget[spec], _QiskitTarget]
]

(* Async primitive submission through qiskit's own SamplerV2 / EstimatorV2 objects, so qiskit
   owns every wire format: circuit serialization, the estimator observable's layout +
   ObservablesArray, the PUB tuple, and the primitive options schema. The job is launched with
   run() but NOT awaited, so this returns immediately with the job id plus the classical-bit ->
   original-qubit decode map for sampler counts. obsTerms is a list of {pauliString, reCoeff,
   imCoeff} triples (estimator) or Null (sampler); primOpts is an already snake_cased option tree
   applied onto the primitive's own options object. Returns an Association or a Failure. *)
qiskitPrimitiveSubmit[
    qc_QiskitCircuit, primitive_String, obsTerms : _List | Null, shots_Integer,
    primOpts_Association, opts : OptionsPattern[Join[{"OptimizationLevel" -> Automatic, "InitialLayout" -> Automatic, "Transpile" -> True}, Options[qiskitInitBackend]]]
] := Enclose @ Block[{
    $primitive = primitive,
    $obsTerms = Replace[obsTerms, None -> Null],
    $shots = shots,
    $primOpts = primOpts,
    (* level 1 already runs the error-aware layout passes against backend.target; users may
       raise to 2 / 3 for more optimization. Default preserves prior submission behavior. *)
    $optLevel = Replace[OptionValue["OptimizationLevel"], Automatic -> 1],
    (* pin the transpiler's initial layout as device physical-qubit indices (0-indexed) that the
       logical qubits map onto; Automatic (or None) lets the error-aware layout pass choose.
       Both sentinels normalize to Null (Python None) once, here in the helper (like OptimizationLevel),
       so a direct caller may pass either. Ignored when $transpile is False. *)
    $initialLayout = Replace[OptionValue["InitialLayout"], (Automatic | None) -> Null],
    (* True: transpile to the backend's ISA (error-aware layout/routing/translation) before
       submitting. False: submit the circuit exactly as given: it must already be ISA for this
       backend, otherwise the V2 primitive rejects it server-side. *)
    $transpile = TrueQ[OptionValue["Transpile"]],
    env
},
    env = Confirm @ qiskitInitBackend[qc, FilterRules[{opts}, Options[qiskitInitBackend]]];
    With[{res = PythonEvaluate[Context[$primitive], "
from qiskit.transpiler import generate_preset_pass_manager
from qiskit.quantum_info import SparsePauliOp
from wolframclient.language import wl

primitive = <* $primitive *>
shots = <* $shots *>
opts = <* $primOpts *> or {}
transpile = <* $transpile *>

# the estimator measures observables, so any terminal measurement the circuit carries is
# stripped; the sampler needs a classical register, added only when the circuit has none
circ = qc
if primitive == 'estimator':
    circ = circ.remove_final_measurements(inplace=False)
elif circ.num_clbits == 0:
    circ = circ.copy()
    circ.measure_all()

# transpile to the backend's ISA (error-aware layout/routing), optionally pinning the initial
# layout to chosen physical qubits; or, when Transpile is False, submit the circuit exactly as
# given (it must already be ISA for this backend, or the V2 primitive rejects it server-side)
if transpile:
    isa = generate_preset_pass_manager(backend=backend, optimization_level=<* $optLevel *>, initial_layout=<* $initialLayout *>).run(circ)
else:
    isa = circ

# classical bit -> qubit map, so sampler counts can be reordered into ascending-qubit order
# regardless of the measure -> clbit layout. When the circuit was transpiled here, isa.layout maps
# each physical qubit back to the original (pre-transpile) logical qubit, so counts decode into the
# user's logical order even when a pinned layout put the circuit on other physical qubits. For a
# pass-through submission (Transpile False) the circuit is taken exactly as given, so counts decode
# by the physical qubit each clbit measures, with no layout re-interpretation.
clbit_to_phys = {}
for ins in isa.data:
    if ins.operation.name == 'measure':
        clbit_to_phys[isa.find_bit(ins.clbits[0]).index] = isa.find_bit(ins.qubits[0]).index
ncl = isa.num_clbits
measured = None
layout = getattr(isa, 'layout', None)
if transpile and layout is not None:
    try:
        final = layout.final_index_layout(filter_ancillas=True)
        phys_to_orig = {phys: orig for orig, phys in enumerate(final)}
        measured = [phys_to_orig.get(clbit_to_phys.get(c), clbit_to_phys.get(c)) for c in range(ncl)]
    except Exception:
        measured = None
if measured is None:
    measured = [clbit_to_phys.get(c, c) for c in range(ncl)]

# forward the user's option tree onto the primitive's own options object (qiskit owns the
# schema). The keys were validated against that schema before transpile by
# ibmValidatePrimitiveOptions, so a plain setattr here is safe.
def apply_opts(o, d):
    for k, v in (d or {}).items():
        if isinstance(v, dict):
            apply_opts(getattr(o, k), v)
        else:
            setattr(o, k, v)

if primitive == 'estimator':
    from qiskit_ibm_runtime import EstimatorV2
    terms = <* $obsTerms *>
    obs = SparsePauliOp.from_list([(t[0], complex(t[1], t[2])) for t in terms])
    # transpiled: the observable is on the original (logical) qubits, so map it onto the physical
    # register through the layout the pass manager produced. Pass-through (Transpile False): the
    # circuit is submitted as given, so the observable is taken on the circuit's own qubits and
    # only widened to the circuit's width, never re-laid-out.
    if transpile:
        obs = obs.apply_layout(isa.layout)
    else:
        obs = obs.apply_layout(None, num_qubits=isa.num_qubits)
    est = EstimatorV2(mode=backend)
    est.options.default_shots = shots
    est.options.dynamical_decoupling.enable = True
    apply_opts(est.options, opts)
    job = est.run([(isa, obs)])
else:
    from qiskit_ibm_runtime import SamplerV2
    smp = SamplerV2(mode=backend)
    smp.options.default_shots = shots
    smp.options.dynamical_decoupling.enable = True
    apply_opts(smp.options, opts)
    job = smp.run([(isa,)])

try:
    jid = job.job_id()
except Exception:
    jid = None

wl.Association(
    wl.Rule('JobID', jid),
    wl.Rule('MeasuredQubits', measured),
    wl.Rule('Primitive', primitive),
    wl.Rule('Backend', getattr(backend, 'name', None))
)
", env]},
        If[FailureQ[res], res, ConfirmBy[res, AssociationQ]]
    ]
]


qc_QiskitCircuit["QASM" | "QASM2", opts___] := QuantumQASM[qc, "Version" -> 2, opts]

qc_QiskitCircuit["QASM3", opts___] := QuantumQASM[qc, "Version" -> 3, opts]


qc_QiskitCircuit["QPY", opts : OptionsPattern[qiskitInitBackend]] := Enclose @ Block[{env},
    env = Confirm @ qiskitInitBackend[qc, opts];
    PythonEvaluate["
from qiskit import qpy
from io import BytesIO
from qiskit import transpile
import zlib

qc = transpile(qc, backend)
bytes = BytesIO()
qpy.dump(qc, bytes)
zlib.compress(bytes.getvalue())
", env]
]

ImportQPY[qpy_ByteArray] := Block[{$qpy = qpy}, PythonEvaluate[Context[$qpy], "
from qiskit import qpy
import zlib
from io import BytesIO
from wolframclient.language import wl
import pickle

qcs = qpy.load(BytesIO(zlib.decompress(<* $qpy *>)))
[wl.Wolfram.QuantumFramework.QiskitCircuit(pickle.dumps(qc)) for qc in qcs]
"]
]

Options[qiskitTranspile] = Join[{"InitialLayout" -> None, "OptimizationLevel" -> 0}, Options[qiskitInitBackend]]

qiskitTranspile[qc_QiskitCircuit, basisGates : {_String...} | None : None, opts : OptionsPattern[]] := Enclose @ Block[{
    $basisGates = Replace[basisGates, None -> Null],
    $initialLayout = Replace[OptionValue["InitialLayout"], None -> Null],
    $optimizationLevel = OptionValue["OptimizationLevel"],
    env = Confirm @ qiskitInitBackend[qc, FilterRules[{opts}, Options[qiskitInitBackend]]]
},
    PythonEvaluate[Context[$basisGates], "
import pickle
from qiskit import transpile
from wolframclient.language import wl

wl.Wolfram.QuantumFramework.QiskitCircuit(
    pickle.dumps(
        transpile(qc,
            backend=backend,
            basis_gates=<* $basisGates *>,
            initial_layout=<* $initialLayout *>,
            optimization_level=<* $optimizationLevel *>
        )
    )
)
", env]
]

qc_QiskitCircuit["Transpile", args___] := qiskitTranspile[qc, args]

qc_QiskitCircuit["Layout", opts : OptionsPattern[qiskitInitBackend]] := Enclose @ Block[{image, env},
    env = Confirm @ qiskitInitBackend[qc, opts];
    image = PythonEvaluate[Context[image], "
from qiskit.visualization import plot_circuit_layout
import PIL
import matplotlib

fig = plot_circuit_layout(qc, backend)
fig.canvas.draw()
image = PIL.Image.frombytes('RGBA', fig.canvas.get_width_height(), fig.canvas.buffer_rgba())
matplotlib.pyplot.close()
image
", env];
    image
]

qc_QiskitCircuit["Validate", opts : OptionsPattern[qiskitInitBackend]]:= Enclose[
    Confirm @ qiskitInitBackend[qc, opts, "FireOpal" -> True, "Provider" -> "IBMProvider"];
    PythonEvaluate["
import fireopal
from qiskit import transpile, qasm2
circuit = transpile(qc, backend)
qc_qasm = qasm2.dumps(qc)
circuit_qasm = qasm2.dumps(circuit)
fireopal.validate(
    circuits=[qc_qasm, circuit_qasm], credentials=fireopal_credentials, backend_name=backend.name
)
", "qctrl"]
]

ImportQASMCircuit::deprecated = "ImportQASMCircuit is deprecated. Use QuantumQASM[source] or QuantumCircuitOperator[source] to import OpenQASM directly.";

ImportQASMCircuit[file_String /; FileExistsQ[file], ___] := (
    Message[ImportQASMCircuit::deprecated];
    QuantumQASM[File[file]]
)

ImportQASMCircuit[str_String, ___] := (
    Message[ImportQASMCircuit::deprecated];
    QuantumQASM[str]
)


(* Formatting *)

QiskitCircuit /: MakeBoxes[qc : QiskitCircuit[_ByteArray], format_] :=
    BoxForm`ArrangeSummaryBox[
        "QiskitCircuit",
        qc,
        Show[Replace[$quantumServiceIcon, Except[_Image] -> QuantumCircuitOperator[{"Fourier"[3]}]["Icon", "GateBackgroundStyle" -> _ -> LightGray, "GateBoundaryStyle" -> _ -> Gray]], ImageSize -> 12],
        {
            {BoxForm`SummaryItem[{"Qubits: ", qc["Qubits"]}]},
            {BoxForm`SummaryItem[{"Depth: ", qc["Depth"]}]}
        },
        {
            {BoxForm`SummaryItem[{"Ops: ", qc["Ops"]}]}
        },
        format,
        "Interpretable" -> Automatic
    ]

