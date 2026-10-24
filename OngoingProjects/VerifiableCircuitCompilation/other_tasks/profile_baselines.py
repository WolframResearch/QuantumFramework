"""Run unchanged DQI/DFTHC baselines with default official problem settings."""

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time

import kickmix as km

HERE = Path(__file__).resolve().parent
SALT = "qf-other-tasks-feasibility-v1"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("dqi","dfthcsos"))
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    sys.path.insert(0,str(args.data/"verifiers"))
    paths = sorted((args.data/"verifiers"/args.task).rglob("*.py"))
    if args.task == "dqi":
        paths.append(args.data/"sample_submission/dqi.kmb")
    else:
        paths.extend(sorted((args.data/"sample_submission/dfthcsos").glob("*")))
        paths.extend((args.data/"verifiers/dfthcsos/FeMoCo76o113e.fixed_shifts.npz",
                      args.data/"verifiers/dfthcsos/FeMoCo76o113e.fcidump"))
    hashes = {str(p.relative_to(args.data)):digest(p) for p in paths if p.is_file()}
    metadata = {"task":args.task,"input_sha256":hashes,"official_code_modified":False,
                "generator_sha256":digest(Path(__file__)),"salt":SALT,
                "kickmix_version":km.__version__,"coverage":"DEFAULT_FULL_OFFICIAL_BASELINE"}
    output = HERE/(args.task+"-baseline.json")
    output.write_text(json.dumps({**metadata,"status":"RUNNING"},indent=2)+"\n")
    started = time.monotonic()
    print(json.dumps({"task":args.task,"status":"RUNNING"}),flush=True)
    try:
        if args.task == "dqi":
            verifier = importlib.import_module("dqi.verify_submission_dqi")
            circuit = km.Circuit.from_file(args.data/"sample_submission/dqi.kmb")
            static = circuit.operation_statistics()
            metadata["static"] = {"num_qubits":circuit.num_qubits,"num_ops":len(circuit),
                "toffoli_count":static.num_toffoli_like,"operation_counts":dict(static)}
            output.write_text(json.dumps({**metadata,"status":"RUNNING"},indent=2)+"\n")
            print(json.dumps({"static":metadata["static"]}),flush=True)
            result = verifier.verify(circuit,verifier.ProblemSpecDQI.benchmark_instance(),
                                     salt=SALT,max_workers=args.workers)
        else:
            verifier = importlib.import_module("dfthcsos.verify")
            result = verifier.verify(args.data/"sample_submission/dfthcsos",salt=SALT)
        if any(digest(args.data/name)!=sha for name,sha in hashes.items()):
            raise RuntimeError("An official verification input changed during profiling.")
        report = {**metadata,**result,"elapsed_seconds":time.monotonic()-started}
        output.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps(report,indent=2),flush=True)
    except Exception as error:
        output.write_text(json.dumps({**metadata,"status":"FAIL","failure_type":type(error).__name__,
                          "elapsed_seconds":time.monotonic()-started},indent=2)+"\n")
        raise


if __name__=="__main__":
    main()
