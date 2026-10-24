"""Compiled arithmetic fragments, not a candidate DQI decoder."""

import hashlib
import json
from pathlib import Path

import kickmix as km
import numpy as np

HERE=Path(__file__).resolve().parent
GF=km.GF2Field(12)


def frobenius(builder,power,target):
    """Compile the invertible binary Frobenius map using elementary row gates."""
    matrix=np.array([[(GF.frobenius(1<<j,power%12)>>i)&1 for j in range(12)]
                     for i in range(12)],dtype=np.uint8)
    operations=[]
    for col in range(12):
        pivot=next(row for row in range(col,12) if matrix[row,col])
        if pivot!=col:
            matrix[[col,pivot]]=matrix[[pivot,col]];operations.append(("swap",col,pivot))
        for row in range(12):
            if row!=col and matrix[row,col]:
                matrix[row]^=matrix[col];operations.append(("cx",col,row))
    assert np.array_equal(matrix,np.eye(12,dtype=np.uint8))
    for op,a,b in reversed(operations):getattr(builder,op)(target[a],target[b])


def xor_linear(builder,matrix,source,target):
    for out,row in enumerate(matrix):
        for bit,value in enumerate(row):
            if value:builder.cx(source[bit],target[out])


def six_product(builder,a,b,target):
    certificate=json.loads((HERE/"symbolic.json").read_text())["six_product"]
    temp=builder.alloc_qubits(1)[0]
    for term,form in enumerate(certificate["input_forms"]):
        indices=[j for j,v in enumerate(form) if v]
        first,*rest=indices
        for j in rest:builder.cx(a[j],a[first]);builder.cx(b[j],b[first])
        builder.init_and(a[first],b[first],temp)
        for j,row in enumerate(certificate["output_forms"]):
            if row[term]:builder.cx(temp,target[j])
        builder.del_and(a[first],b[first],temp)
        for j in reversed(rest):builder.cx(b[j],b[first]);builder.cx(a[j],a[first])
    builder.free(km.array([temp]))


def tower_multiply(builder,a,b,target,plan):
    degree=plan["subfield_degree"];field=km.GF2Field(degree)
    left,right,product=[builder.alloc_qubits(degree) for _ in range(3)]
    for evaluate,output in zip(plan["evaluation_maps"],plan["output_maps"]):
        xor_linear(builder,evaluate,a,left);xor_linear(builder,evaluate,b,right)
        for q in product:builder.reset(q)
        if degree==3:six_product(builder,left,right,product)
        else:builder.ixor_gf2_mul(left,right,target=product,field=field)
        xor_linear(builder,output,product,target)
        builder.del_gf2_mul(left,right,target=product,field=field)
        xor_linear(builder,evaluate,b,right);xor_linear(builder,evaluate,a,left)
    for reg in (left,right,product):builder.free(reg)


def multiply(builder,a,b,target,plan=None,control=None):
    if plan is None:
        builder.ixor_gf2_mul(a,b,target=target,field=GF,control=True if control is None else control)
    elif control is None:tower_multiply(builder,a,b,target,plan)
    else:
        product=builder.alloc_qubits(12)
        tower_multiply(builder,a,b,product,plan)
        builder.ccx(control,product,target)
        builder.del_gf2_mul(a,b,target=product,field=GF)
        builder.free(product)


def multiplication_circuit(plan=None,controlled=False):
    builder=km.CircuitBuilder()
    a,b,target=[builder.create_quantum_register(12,name=n) for n in ("a","b","target")]
    ctrl=builder.create_quantum_register(1,name="control")[0] if controlled else None
    multiply(builder,a,b,target,plan,ctrl)
    return builder.finish_circuit()


def inverse_phase_circuit(plan=None,mask=1):
    builder=km.CircuitBuilder();source=builder.create_quantum_register(12,name="input")
    if plan is None:
        output,scaffold=builder.init_gf2_inverse_with_scaffold(source,field=GF)
        for bit in range(12):
            if (mask>>bit)&1:builder.z(output[bit])
        builder.del_gf2_inverse_with_scaffold(source,target=output,scaffold=scaffold,field=GF)
        builder.free(output+scaffold)
        return builder.finish_circuit()
    output=builder.alloc_qubits(12)
    nodes={1:source,**{k:builder.alloc_qubits(12) for k in (2,4,8,10)}}
    for prev,new in ((1,2),(2,4),(4,8)):
        builder.cx(nodes[prev],output)
        frobenius(builder,prev,output)
        multiply(builder,nodes[prev],output,nodes[new],plan)
        frobenius(builder,-prev,output)
        builder.cx(nodes[prev],output)
    frobenius(builder,2,nodes[8])
    multiply(builder,nodes[8],nodes[2],nodes[10],plan)
    frobenius(builder,1,nodes[10])
    multiply(builder,nodes[10],source,output,plan)
    frobenius(builder,1,output)
    for bit in range(12):
        if (mask>>bit)&1:builder.z(output[bit])
    frobenius(builder,-1,output)
    builder.del_gf2_mul(nodes[10],source,target=output,field=GF)
    for q in output:builder.reset(q)
    frobenius(builder,-1,nodes[10])
    builder.del_gf2_mul(nodes[8],nodes[2],target=nodes[10],field=GF)
    frobenius(builder,-2,nodes[8])
    for prev,new in ((4,8),(2,4),(1,2)):
        builder.cx(nodes[prev],output)
        frobenius(builder,prev,output)
        builder.del_gf2_mul(nodes[prev],output,target=nodes[new],field=GF)
        frobenius(builder,-prev,output)
        builder.cx(nodes[prev],output)
    for reg in [output,*[nodes[k] for k in (2,4,8,10)]]:builder.free(reg)
    return builder.finish_circuit()


def packed(values,width=12):
    return ((np.asarray(values,dtype=np.uint64)[None,:] >> np.arange(width,dtype=np.uint64)[:,None])&1).astype(bool)


def unpack(values):
    return (values.astype(np.uint64)<<np.arange(len(values),dtype=np.uint64)[:,None]).sum(axis=0)


def check_circuit(circuit,inverse=False,mask=1):
    rng=np.random.default_rng(61283)
    # Every field element plus independent random partners and dirty targets.
    a=np.tile(np.arange(4096,dtype=np.uint64),2)
    sim=km.Simulator(batch_size=len(a),counting_operations=True,seed=931)
    public=set()
    if inverse:
        reg=circuit.register_data["input"];public.update(reg)
        sim.write_across_shots(reg,packed(a));sim.do(circuit)
        expected=np.array([bool((GF.invert(int(v))&mask).bit_count()%2) for v in a])
        assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
        assert np.array_equal(unpack(sim.read_across_shots(reg)),a)
    else:
        b=rng.integers(0,4096,len(a),dtype=np.uint64);target=rng.integers(0,4096,len(a),dtype=np.uint64)
        control=rng.integers(0,2,len(a),dtype=np.uint64)
        for name,values in (("a",a),("b",b),("target",target)):
            reg=circuit.register_data[name];public.update(reg);sim.write_across_shots(reg,packed(values))
        if "control" in circuit.register_data:
            reg=circuit.register_data["control"];public.update(reg);sim.write_across_shots(reg,packed(control,1))
        else:control[:]=1
        sim.do(circuit)
        if "control" in circuit.register_data:
            assert np.array_equal(unpack(sim.read_across_shots(circuit.register_data["control"])),control)
        expected=target^np.array([GF.mul(int(x),int(y)) for x,y in zip(a,b)],dtype=np.uint64)*control
        for name,values in (("a",a),("b",b),("target",expected)):
            assert np.array_equal(unpack(sim.read_across_shots(circuit.register_data[name])),values),name
        assert not sim.read_phase_flipped_across_shots().any()
    scratch=km.array([km.q(q) for q in range(circuit.num_qubits) if km.q(q) not in public])
    assert not sim.read_across_shots(scratch).any()
    return {"shots":len(a),"all_field_elements_covered":True,"scratch_clean":True,
            "input_restored":True,"phase_correct":True}


def main():
    symbolic=json.loads((HERE/"symbolic.json").read_text())
    assert symbolic["six_product"]["certificate"]
    records=[]
    for plan in [None,*symbolic["tower_plans"]]:
        if plan is not None:assert all(plan[k] for k in ("basis_certificate","interpolation_certificate","tensor_certificate"))
        name="native" if plan is None else "subfield-"+str(plan["subfield_degree"])
        for kind,controlled in (("multiply",False),("controlled_multiply",True),("inverse_phase",False)):
            circuit=inverse_phase_circuit(plan) if kind=="inverse_phase" else multiplication_circuit(plan,controlled)
            checks=check_circuit(circuit,kind=="inverse_phase");stats=circuit.operation_statistics()
            records.append({"name":name,"kind":kind,"num_qubits":circuit.num_qubits,
                "toffoli_count":stats.num_toffoli_like,"num_ops":len(circuit),"checks":checks})
    report={"status":"ARITHMETIC_FRAGMENTS_ONLY_NOT_DQI_VERIFICATION","fragments":records,
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "symbolic_sha256":hashlib.sha256((HERE/"symbolic.json").read_bytes()).hexdigest()}
    (HERE/"field-results.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
