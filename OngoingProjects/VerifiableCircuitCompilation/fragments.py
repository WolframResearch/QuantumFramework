"""Export small, real kickmix arithmetic fragments for independent QF checks."""
from fractions import Fraction
import io
import json
from pathlib import Path

import kickmix as km

HERE = Path(__file__).resolve().parent


def describe(circuit):
    buffer = io.StringIO()
    circuit.to_file(buffer, format="kmx")
    lines = [line for line in map(str.strip, buffer.getvalue().splitlines())
             if line and not line.startswith(("APPEND_TO_REGISTER", "REGISTER", "INFO_TAG"))]
    return {"NumQubits": circuit.num_qubits, "Lines": lines,
            "Measurements": sum(line.startswith("HMR ") for line in lines),
            "Toffolis": circuit.operation_statistics().num_toffoli_like}


def main():
    builder = km.CircuitBuilder()
    inputs = builder.create_quantum_register(2, name="input")
    scratch = builder.alloc_qubits(1)[0]
    builder.init_and(inputs[0], inputs[1], scratch)
    builder.z_pow(scratch, Fraction(3, 16))
    builder.del_and(inputs[0], inputs[1], scratch)
    phase = describe(builder.finish_circuit())

    builder = km.CircuitBuilder()
    lhs = builder.create_quantum_register(2, name="a")
    rhs = builder.create_quantum_register(2, name="b")
    builder.iadd(rhs, target=lhs)
    addition = describe(builder.finish_circuit())
    (HERE / "results/fragments.json").write_text(json.dumps({"AndPhase": phase,
                                                           "Addition": addition}, indent=2) + "\n")


if __name__ == "__main__":
    main()
