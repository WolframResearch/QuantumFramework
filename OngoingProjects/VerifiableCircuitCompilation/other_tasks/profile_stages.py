"""Extract native DQI flame-chart profiles and conditional substitution models.

SVG widths are rounded display values; stage Toffoli counts are approximate.
The displayed group hit counts are exact. Native degree-12 gf_mul is 63
Toffolis, verified by the separately compiled fragment.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import kickmix as km

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True)
    args=parser.parse_args()
    circuit=km.Circuit.from_file(args.data/"sample_submission/dqi.kmb")
    baseline=json.loads((HERE/"dqi-baseline.json").read_text())
    fields=json.loads((HERE/"field-results.json").read_text())
    if baseline["status"]!="PASS":raise ValueError("Full baseline verification is incomplete.")
    svg=circuit.flame_chart_svg();root=ET.fromstring(svg)
    last_rectangle=None;groups=[];gf_calls=0
    for element in root:
        tag=element.tag.split("}")[-1]
        if tag=="rect" and element.get("fill")=="red":last_rectangle=element
        if tag!="text":continue
        label="".join(element.itertext())
        match=re.fullmatch(r"gf_mul(?: \(x(\d+)\))?",label)
        if match:gf_calls+=int(match.group(1) or 1)
        if element.get("y")=="7":
            groups.append({"label":label,"approximate_static_toffoli":
                float(last_rectangle.get("width"))*baseline["max_toffoli_count"]/1000})
    fragments={ (r["name"],r["kind"]):r for r in fields["fragments"] }
    native=fragments[("native","multiply")]["toffoli_count"]
    scenarios=[]
    for name in ("subfield-6","subfield-4","subfield-3"):
        new=fragments[(name,"multiply")]["toffoli_count"]
        delta_q=fragments[(name,"inverse_phase")]["num_qubits"]-fragments[("native","inverse_phase")]["num_qubits"]
        modeled_t=baseline["max_toffoli_count"]-gf_calls*(native-new)
        modeled_q=baseline["num_qubits"]+delta_q
        cost=baseline["repetition_factor"]*modeled_q**2*modeled_t
        scenarios.append({"name":name,"static_toffoli_model":modeled_t,
            "assumed_qubits":modeled_q,"assumed_repetition_factor":baseline["repetition_factor"],
            "modeled_objective":cost,"modeled_improvement":baseline["score"]/cost})
    report={"status":"PROFILING_AND_CONDITIONAL_MODELS_ONLY","stages":groups,
            "gf_multiplications_static":gf_calls,"native_toffoli_per_multiplication":native,
            "multiplication_toffoli_fraction":gf_calls*native/baseline["max_toffoli_count"],
            "substitution_scenarios":scenarios,
            "limitations":["No modified decoder was built or verified.",
                "Static counts include classical conditions; baseline actual mean is separately recorded.",
                "Modeled peak qubits assume the fragment scratch increase adds to the baseline peak.",
                "The repetition factor assumes zero errors on the fixed full verification protocol."],
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "baseline_sha256":hashlib.sha256((HERE/"dqi-baseline.json").read_bytes()).hexdigest(),
            "field_results_sha256":hashlib.sha256((HERE/"field-results.json").read_bytes()).hexdigest()}
    (HERE/"dqi-profile.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":main()
