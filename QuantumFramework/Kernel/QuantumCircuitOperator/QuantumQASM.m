Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumQASM"]
PackageExport["QiskitTarget"]



(* Option-key normalization.

   The WL surface uses CamelCase option keys; qiskit's transpile / Target parameters
   are case-sensitive snake_case. A deterministic transform inserts "_" at each
   lower/digit -> upper boundary and lowercases, so "CouplingMap", "couplingMap" and
   "coupling_map" all map to qiskit's "coupling_map". This is NOT a curated dictionary,
   so a future qiskit parameter maps automatically; the normalized key is then validated
   against the live inspect.signature, so anything that does not resolve is a Failure. *)

camelToSnake[s_String] := ToLowerCase @ StringReplace[s, RegularExpression["(?<=[a-z0-9])([A-Z])"] -> "_$1"]

normalizeKeys[opts_Association] := KeyMap[camelToSnake, opts]


(* Transpile core: forward an Association of options straight to qiskit.transpile as
   **kwargs, validated at call time against the live signature (version-resilient).
   Returns a QiskitCircuit (pickled transpiled circuit) or a Failure. *)

qqasmTranspileBytes[bytes_ByteArray, opts_Association] := Enclose @ Block[{
    pythonBytes = bytes,
    pythonOpts = opts
},
    ConfirmBy[
        PythonEvaluate[Context[pythonBytes], "
import pickle, inspect, qiskit
from qiskit import transpile
from qiskit.transpiler import CouplingMap, Target
from qiskit.circuit import Measure, Reset
from wolframclient.language import wl

circuit = pickle.loads(<* pythonBytes *>)
kwargs = dict(<* pythonOpts *>)

accepted = sorted(inspect.signature(transpile).parameters)
unknown = sorted(k for k in kwargs if k not in accepted)
if unknown:
    result = wl.Failure('QuantumQASM', {
        'MessageTemplate': 'Unknown transpile option(s): `1`. Accepted (qiskit `2`): `3`.',
        'MessageParameters': [unknown, qiskit.__version__, accepted],
        'Unknown': unknown,
        'Accepted': accepted,
        'QiskitVersion': qiskit.__version__})
else:
    def build_target(spec):
        spec = dict(spec)
        if spec.get('coupling_map') is not None and not isinstance(spec['coupling_map'], CouplingMap):
            spec['coupling_map'] = CouplingMap([tuple(e) for e in spec['coupling_map']])
        if not hasattr(Target, 'from_configuration'):
            raise RuntimeError('qiskit Target.from_configuration is unavailable')
        t = Target.from_configuration(**spec)
        nq = t.num_qubits
        # measure / reset are universal; from_configuration omits them for a bare unitary
        # basis, which makes measured circuits fail at HighLevelSynthesis. Add idempotently.
        try:
            if 'measure' not in t.operation_names:
                t.add_instruction(Measure(), {(q,): None for q in range(nq)})
            if 'reset' not in t.operation_names:
                t.add_instruction(Reset(), {(q,): None for q in range(nq)})
        except Exception:
            pass
        return t
    try:
        if kwargs.get('coupling_map') is not None and not isinstance(kwargs['coupling_map'], CouplingMap):
            kwargs['coupling_map'] = CouplingMap([tuple(e) for e in kwargs['coupling_map']])
        if isinstance(kwargs.get('target'), dict):
            kwargs['target'] = build_target(kwargs['target'])
        out = transpile(circuit, **kwargs)
        result = wl.Wolfram.QuantumFramework.QiskitCircuit(pickle.dumps(out))
    except Exception as e:
        result = wl.Failure('QuantumQASM', {
            'MessageTemplate': 'qiskit transpile failed: `1`',
            'MessageParameters': [type(e).__name__ + ': ' + str(e)],
            'PythonError': str(e)})
result
"],
        MatchQ[_QiskitCircuit | _Failure]
    ]
]


(* Faithful OpenQASM dump: serialize the circuit as-is, no transpile. Guarded so a
   non-OPENQASM result is a Failure rather than a malformed string. *)

qqasmDumpBytes[bytes_ByteArray, version : (2 | 3)] := Enclose @ Block[{
    pythonBytes = bytes,
    pythonVersion = version
},
    ConfirmBy[
        PythonEvaluate[Context[pythonBytes], "
import pickle
from wolframclient.language import wl

circuit = pickle.loads(<* pythonBytes *>)
try:
    if <* pythonVersion *> == 3:
        from qiskit import qasm3
        s = qasm3.dumps(circuit)
    else:
        from qiskit import qasm2
        s = qasm2.dumps(circuit)
    if isinstance(s, str) and s.lstrip().startswith('OPENQASM'):
        result = s
    else:
        result = wl.Failure('QuantumQASM', {
            'MessageTemplate': 'OpenQASM export did not produce an OPENQASM string.'})
except Exception as e:
    result = wl.Failure('QuantumQASM', {
        'MessageTemplate': 'OpenQASM`1` export failed: `2`',
        'MessageParameters': [<* pythonVersion *>, type(e).__name__ + ': ' + str(e)],
        'PythonError': str(e)})
result
"],
        MatchQ[_String | _Failure]
    ]
]


(* QiskitTarget: a validated, normalized device-spec carrier (NOT a pickled Target).
   The live qiskit Target is rebuilt from the spec at point-of-use inside the transpile
   core, so there is no cross-version pickle to break. *)

validateTargetSpec[spec_Association] := Enclose @ Block[{pythonSpec = spec},
    ConfirmBy[
        PythonEvaluate[Context[pythonSpec], "
import inspect
from wolframclient.language import wl
try:
    from qiskit.transpiler import Target
except Exception as e:
    result = wl.Failure('QiskitTarget', {
        'MessageTemplate': 'qiskit Target is unavailable: `1`',
        'MessageParameters': [str(e)]})
else:
    if not hasattr(Target, 'from_configuration'):
        result = wl.Failure('QiskitTarget', {
            'MessageTemplate': 'qiskit Target.from_configuration is unavailable.'})
    else:
        accepted = sorted(inspect.signature(Target.from_configuration).parameters)
        spec = dict(<* pythonSpec *>)
        unknown = sorted(k for k in spec if k not in accepted)
        if unknown:
            result = wl.Failure('QiskitTarget', {
                'MessageTemplate': 'Unknown Target option(s): `1`. Accepted: `2`.',
                'MessageParameters': [unknown, accepted],
                'Unknown': unknown,
                'Accepted': accepted})
        else:
            result = True
result
"],
        MatchQ[True | _Failure]
    ]
]

QiskitTarget[spec_Association] /; ! KeyExistsQ[spec, "$QiskitTargetSpec"] := Module[{normSpec = normalizeKeys[spec], chk},
    chk = validateTargetSpec[normSpec];
    If[FailureQ[chk], chk, QiskitTarget[<|"$QiskitTargetSpec" -> normSpec|>]]
]

QiskitTarget[data_Association]["Spec"] /; KeyExistsQ[data, "$QiskitTargetSpec"] := data["$QiskitTargetSpec"]

QiskitTarget[data_Association]["NumQubits"] /; KeyExistsQ[data, "$QiskitTargetSpec"] :=
    Lookup[data["$QiskitTargetSpec"], "num_qubits", Missing["NotAvailable"]]

QiskitTarget[data_Association]["CouplingMap"] /; KeyExistsQ[data, "$QiskitTargetSpec"] :=
    Lookup[data["$QiskitTargetSpec"], "coupling_map", Missing["NotAvailable"]]

QiskitTarget[data_Association]["OperationNames"] /; KeyExistsQ[data, "$QiskitTargetSpec"] :=
    Union[Lookup[data["$QiskitTargetSpec"], "basis_gates", {}], {"measure", "reset"}]

(* Summary-box icon: the monochrome Qiskit globe logo (Wikimedia Commons,
   File:Qiskit-Logo.svg, public domain), rasterized to a 45 px transparent-background image and
   decompressed once at load. *)
$qiskitTargetIcon = Uncompress["1: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"];

QiskitTarget /: MakeBoxes[qt : QiskitTarget[data_Association] /; KeyExistsQ[data, "$QiskitTargetSpec"], fmt_] :=
    BoxForm`ArrangeSummaryBox[
        "QiskitTarget",
        qt,
        Show[$qiskitTargetIcon, ImageSize -> 18],
        {
            {BoxForm`SummaryItem[{"Qubits: ", qt["NumQubits"]}]},
            {BoxForm`SummaryItem[{"Gates: ", qt["OperationNames"]}]}
        },
        {
            {BoxForm`SummaryItem[{"Couplings: ", Replace[qt["CouplingMap"], {l_List :> Length[l], _ -> 0}]}]}
        },
        fmt
    ]


(* QuantumQASM: circuit -> OpenQASM string. Faithful dump by default; transpile to a
   native basis only when a basis / Target / coupling map / any transpile option is
   given. Returns an OpenQASM String or a Failure. *)

QuantumQASM[qco_QuantumCircuitOperator, args___] := With[{qk = qco["Qiskit"]},
    If[MatchQ[qk, _QiskitCircuit], QuantumQASM[qk, args], qk]
]

QuantumQASM[qc_QiskitCircuit, basisGates : {___String} | Automatic : Automatic, opts_Association] :=
    QuantumQASM[qc, basisGates, Sequence @@ Normal[opts]]

QuantumQASM[qc_QiskitCircuit, basisGates : {___String} | Automatic : Automatic, opts___Rule] := Module[{
    rawOpts = Association[opts],
    normOpts, version, bytes, tcirc
},
    normOpts = normalizeKeys[rawOpts];
    If[ Length[normOpts] =!= Length[rawOpts],
        Return @ Failure["QuantumQASM", <|"MessageTemplate" ->
            "Ambiguous duplicate option after key normalization (for example CouplingMap and coupling_map together)."|>]
    ];
    version = Lookup[normOpts, "version", Lookup[normOpts, "qasm_version", 3]];
    normOpts = KeyDrop[normOpts, {"version", "qasm_version"}];
    If[ basisGates =!= Automatic,
        If[ KeyExistsQ[normOpts, "basis_gates"],
            Return @ Failure["QuantumQASM", <|"MessageTemplate" ->
                "Native gates given both positionally and as the BasisGates option."|>]
        ];
        normOpts = Prepend[normOpts, "basis_gates" -> basisGates]
    ];
    If[ KeyExistsQ[normOpts, "target"],
        normOpts["target"] = Replace[normOpts["target"], {
            tgt_QiskitTarget :> tgt["Spec"],
            a_Association :> normalizeKeys[a]
        }]
    ];
    bytes = qc["Bytes"];
    If[ Length[normOpts] > 0,
        tcirc = qqasmTranspileBytes[bytes, normOpts];
        If[FailureQ[tcirc], Return @ tcirc];
        qqasmDumpBytes[tcirc["Bytes"], version],
        qqasmDumpBytes[bytes, version]
    ]
]
