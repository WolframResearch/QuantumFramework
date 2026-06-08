(* Prototype for a version-resilient Qiskit transpile + faithful QASM export layer.

   This file does NOT edit the kernel. It defines standalone helper functions in
   the Global context that demonstrate the proposed mechanism on top of the live
   QiskitCircuit object. Load the paclet first, then Get this file.

   Proposed home in the kernel: QuantumFramework/Kernel/QuantumCircuitOperator/Qiskit.m
   (the prototype functions map 1:1 onto the methods proposed in the design doc:
     qTranspile  -> qc_QiskitCircuit["Transpile", optsAssoc]
     qDumpQASM   -> qc_QiskitCircuit["QASM3" | "QASM2"]   (faithful, backend-free)
     qTarget     -> QiskitTarget[spec]                     (new compound object)). *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"];
Needs["Wolfram`QuantumFramework`"];

$pe = Wolfram`QuantumFramework`PackageScope`PythonEvaluate;

(* Option-key style: full CamelCase WL surface, mapped to qiskit snake_case by an
   ALGORITHMIC transform (not a hand-maintained dictionary) so that any future
   qiskit parameter auto-maps with no QF change. A key already in snake_case (or
   all-lowercase) passes through unchanged, so raw qiskit names also work. *)
camelToSnake[s_String] := ToLowerCase @ StringReplace[s, RegularExpression["(?<=.)([A-Z])"] -> "_$1"]
normalizeKeys[opts_Association] := KeyMap[camelToSnake, opts]

(* ------------------------------------------------------------------ *)
(* The version-resilience core: forward an Association of qiskit         *)
(* transpile parameters straight as **kwargs, validated against the     *)
(* live inspect.signature(transpile). Nothing about the option list is  *)
(* hardcoded in WL; a new/removed qiskit parameter needs no QF change.   *)
(* ------------------------------------------------------------------ *)

(* Keys whose values are rich Python objects needing coercion before the call.
   Everything else passes through verbatim. *)
qTranspile[qc_QiskitCircuit, opts_Association : <||>] := Enclose @ Block[{$b = qc["Bytes"], $opts = normalizeKeys[opts], result},
    result = $pe["Global`", "
import pickle, inspect
from qiskit import transpile
from qiskit.transpiler import CouplingMap, Target
from wolframclient.language import wl

kwargs = dict(<* $opts *>)

# --- version-resilience guard: accepted set read live from the installed qiskit ---
accepted = set(inspect.signature(transpile).parameters)
unknown = [k for k in kwargs if k not in accepted]
if unknown:
    result = wl.Failure('QiskitTranspileOption', wl.Association(
        wl.Rule('MessageTemplate', 'Unknown transpile option(s): `1`. Accepted (qiskit `2`): `3`.'),
        wl.Rule('MessageParameters', [unknown, __import__('qiskit').__version__, sorted(accepted)]),
        wl.Rule('Unknown', unknown)))
else:
    # --- coercion layer: only a handful of keys carry rich Python objects ---
    cm = kwargs.get('coupling_map')
    if cm is not None and not isinstance(cm, CouplingMap):
        kwargs['coupling_map'] = CouplingMap([tuple(e) for e in cm])
    tg = kwargs.get('target')
    if isinstance(tg, (bytes, bytearray)):          # a QiskitTarget byte-blob
        kwargs['target'] = pickle.loads(tg)
    try:
        out = transpile(pickle.loads(<* $b *>), **kwargs)
        result = wl.Wolfram.QuantumFramework.QiskitCircuit(pickle.dumps(out))
    except Exception as e:
        result = wl.Failure('QiskitTranspile', wl.Association(
            wl.Rule('MessageTemplate', 'qiskit transpile failed: `1`'),
            wl.Rule('MessageParameters', [type(e).__name__ + ': ' + str(e)]),
            wl.Rule('PythonError', str(e))))
result
"];
    ConfirmBy[result, MatchQ[#, _QiskitCircuit | _Failure] &]
]

(* convenience: accept a flat option sequence too *)
qTranspile[qc_QiskitCircuit, opts : (_Rule | _RuleDelayed) ..] := qTranspile[qc, Association[opts]]

(* ------------------------------------------------------------------ *)
(* Faithful QASM export: dump the circuit AS-IS, never transpile-first. *)
(* Guards the output so a non-OPENQASM result becomes a Failure rather   *)
(* than a fake string.                                                   *)
(* ------------------------------------------------------------------ *)

qDumpQASM[qc_QiskitCircuit, version : (2 | 3) : 3] := Enclose @ Block[{$b = qc["Bytes"], $v = version, result},
    result = $pe["Global`", "
import pickle
from wolframclient.language import wl
qc = pickle.loads(<* $b *>)
try:
    if <* $v *> == 3:
        from qiskit import qasm3
        s = qasm3.dumps(qc)
    else:
        from qiskit import qasm2
        s = qasm2.dumps(qc)
    result = s if isinstance(s, str) and s.lstrip().startswith('OPENQASM') else wl.Failure(
        'QiskitQASMExport', wl.Association(wl.Rule('MessageTemplate', 'export did not produce OPENQASM')))
except Exception as e:
    result = wl.Failure('QiskitQASMExport', wl.Association(
        wl.Rule('MessageTemplate', 'qiskit QASM`1` export failed: `2`'),
        wl.Rule('MessageParameters', [<* $v *>, type(e).__name__ + ': ' + str(e)])))
result
"];
    ConfirmBy[result, MatchQ[#, _String | _Failure] &]
]

(* ------------------------------------------------------------------ *)
(* QiskitTarget: the version-stable device descriptor.                  *)
(* Forwards a spec Association straight to Target.from_configuration via *)
(* the same validate-against-live-signature mechanism. Returns a pickled *)
(* Target blob usable as the "target" transpile option.                  *)
(* ------------------------------------------------------------------ *)

QiskitTarget[spec_Association] := Enclose @ Block[{$spec = normalizeKeys[spec], result},
    result = $pe["Global`", "
import pickle, inspect
from qiskit.transpiler import Target, CouplingMap
from wolframclient.language import wl

spec = dict(<* $spec *>)
accepted = set(inspect.signature(Target.from_configuration).parameters)
unknown = [k for k in spec if k not in accepted]
if unknown:
    result = wl.Failure('QiskitTargetOption', wl.Association(
        wl.Rule('MessageTemplate', 'Unknown Target option(s): `1`. Accepted: `2`.'),
        wl.Rule('MessageParameters', [unknown, sorted(accepted)])))
else:
    cm = spec.get('coupling_map')
    if cm is not None and not isinstance(cm, CouplingMap):
        spec['coupling_map'] = CouplingMap([tuple(e) for e in cm])
    try:
        t = Target.from_configuration(**spec)
        result = wl.Global.QiskitTargetObject(pickle.dumps(t))
    except Exception as e:
        result = wl.Failure('QiskitTarget', wl.Association(
            wl.Rule('MessageTemplate', 'Target.from_configuration failed: `1`'),
            wl.Rule('MessageParameters', [type(e).__name__ + ': ' + str(e)])))
result
"];
    ConfirmBy[result, MatchQ[#, _QiskitTargetObject | _Failure] &]
]

QiskitTargetObject[bytes_ByteArray]["Bytes"] := bytes
QiskitTargetObject /: MakeBoxes[t_QiskitTargetObject, fmt_] := ToBoxes["QiskitTarget[…]", fmt]
