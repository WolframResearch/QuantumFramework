"""Standalone four-address signed-Givens XOR loader.

SPDX-License-Identifier: Apache-2.0
Adapted from Dmitrii Gluzdov's public v8 competition notebook, gauge_qrom.py
and the attributed banked_lab.py. See NOTICE.md and LICENSE. This port
removes notebook/source-module dependencies and adds boundary validation.
"""
from functools import lru_cache
from math import ceil

import kickmix as km
import numpy as np

def batch_deltas(words, J=5):
    words = np.asarray(words)
    if (words.ndim != 2 or 0 in words.shape or
            not np.issubdtype(words.dtype, np.integer) or
            np.any(words < 0) or np.any(words >= 1024) or
            not isinstance(J, int) or isinstance(J, bool) or J < 1):
        raise ValueError("Expected a nonempty ten-bit integer table and positive batch size.")
    words = words.astype(np.int64)
    nv, na = words.shape
    padded = np.zeros((nv, ceil(na / J) * J), dtype=np.int64)
    padded[:, :na] = words
    batches = padded.reshape(nv, -1, J).transpose(1, 0, 2)
    result = batches.copy()
    result[1:] ^= batches[:-1]
    return result


def pack(row, brot):
    return sum(int(value) << (j * brot) for j, value in enumerate(row))


def mobius_transform(table):
    """Packed-output Boolean truth table <-> ANF coefficient table over GF(2)."""
    result = list(map(int, table))
    if not result or len(result) & (len(result) - 1):
        raise ValueError("truth-table length must be a power of two")
    n = (len(result) - 1).bit_length()
    for bit in range(n):
        step = 1 << bit
        for start in range(0, len(result), 2 * step):
            for offset in range(step):
                result[start + step + offset] ^= result[start + offset]
    return result


def del_lookup_anf(builder, *, table, address, target, low_bits=None):
    """Clear table[address] with exact conditional ANF phase correction.

    `target` must contain the given lookup value and remain allocated. Extra
    scratch is allocated only when the reset target is too narrow. Clean target
    wires are reused internally; all target and extra scratch finish at zero.
    """
    import kickmix as km
    address = list(address)
    target = list(target)
    n = len(address)
    if len(table) != 1 << n:
        raise ValueError("lookup table must cover the complete address space")
    if n < 2:
        raise ValueError("split ANF cleanup expects at least two address bits")
    low_bits = n // 2 if low_bits is None else low_bits
    if not 0 < low_bits < n:
        raise ValueError("both address halves must be nonempty")
    if any(value < 0 or value >= (1 << len(target)) for value in table):
        raise ValueError("lookup value does not fit the target")
    coefficients = mobius_transform(table)
    assert mobius_transform(coefficients) == list(table), "exact ANF identity failed"
    measurements = builder.alloc_bits(len(target))
    builder.hmr(km.array(target), measurements)
    parity = builder.alloc_bits(1)[0]

    low_count = 1 << low_bits
    high_bits = n - low_bits
    retained_masks = [mask for mask in range(1, low_count) if mask.bit_count() > 1]
    scratch_count = len(retained_masks) + high_bits - 1
    extra = builder.alloc_qubits(max(0, scratch_count - len(target)))
    scratch = target + list(extra)
    low = {1 << bit: address[bit] for bit in range(low_bits)}
    for mask, qubit in zip(retained_masks, scratch):
        bit = mask.bit_length() - 1
        parent = mask ^ (1 << bit)
        builder.init_and(low[parent], address[bit], qubit)
        low[mask] = qubit
    carries = scratch[len(retained_masks):scratch_count]

    def correct(high_mask, high_qubit):
        for low_mask in range(low_count):
            coefficient = coefficients[low_mask | (high_mask << low_bits)]
            if not coefficient:
                continue
            builder.bit_store(parity, False)
            while coefficient:
                bit = (coefficient & -coefficient).bit_length() - 1
                builder.cx(measurements[bit], parity)
                coefficient &= coefficient - 1
            if not low_mask and not high_mask:
                builder.neg_if(parity)
                continue
            builder.push_condition(parity)
            if not low_mask:
                builder.z(high_qubit)
            elif not high_mask:
                builder.z(low[low_mask])
            else:
                builder.cz(low[low_mask], high_qubit)
            builder.pop_condition()

    def walk(mask, highest_bit, control):
        correct(mask, control)
        for bit in range(highest_bit + 1, high_bits):
            child_mask = mask | (1 << bit)
            if not mask:
                walk(child_mask, bit, address[low_bits + bit])
            else:
                child = carries[mask.bit_count() - 1]
                builder.init_and(control, address[low_bits + bit], child)
                walk(child_mask, bit, child)
                builder.del_and(control, address[low_bits + bit], child)

    walk(0, -1, None)
    for mask in reversed(retained_masks):
        bit = mask.bit_length() - 1
        parent = mask ^ (1 << bit)
        builder.del_and(low[parent], address[bit], low[mask])
    builder.free(extra)
    return {"address_bits": n, "low_bits": low_bits, "high_bits": high_bits,
            "retained_products": len(retained_masks), "traversal_carries": high_bits - 1,
            "anf_toffoli": (1 << low_bits) - low_bits - 1 + (1 << high_bits) - high_bits - 1,
            "scratch_count": scratch_count, "anf_table_identity": "PASS"}



def xor_and(b,left,right,targets):
    targets=list(targets)
    if not targets:return
    q=targets[0]
    for t in targets[1:]:b.cx(q,t)
    b.ccx(left,right,q)
    for t in reversed(targets[1:]):b.cx(q,t)

def pebbled_init(b,table,address,outputs,scratch):
    address,outputs,scratch=map(list,(address,outputs,scratch))
    H,n=len(table),len(address)
    coefficients=mobius_transform(list(table)+[0]*((1<<n)-H))
    coefficients[H:]=[0]*((1<<n)-H)
    completed=mobius_transform(coefficients)
    assert completed[:H]==table
    stats={'monomial_ands':0,'recomputed_ands':0,'pebble_sites':0,'direct_leaf_ands':0,'scratch':len(scratch)}
    @lru_cache(None)
    def children(mask,last):return tuple((mask|(1<<j),j) for j in range(last+1,n) if (mask|(1<<j))<H)
    @lru_cache(None)
    def need(mask,last):
        return max([1+need(c,j) for c,j in children(mask,last) if children(c,j)]+[0])
    def targets(mask):
        v=coefficients[mask]
        return [outputs[j] for j in range(len(outputs)) if (v>>j)&1]
    def walk(mask,last,control,path,free):
        erased=[]
        if mask.bit_count()==5 and need(mask,last)>len(free):
            m=need(mask,last)-len(free)
            assert 0<m<=len(path)-1
            erased=path[:m]
            for a,c,t in reversed(erased):b.del_and(a,c,t);free.append(t)
            stats['pebble_sites']+=1
        ts=targets(mask)
        if ts:
            if control is None:b.x(km.array(ts))
            else:b.cx(control,km.array(ts))
        for child,j in children(mask,last):
            if not mask:
                walk(child,j,address[j],path,free)
            elif not children(child,j):
                ts=targets(child)
                if ts:xor_and(b,control,address[j],ts);stats['direct_leaf_ands']+=1
            else:
                assert free,('scratch exhausted',mask,j,len(path))
                t=free.pop();b.init_and(control,address[j],t);stats['monomial_ands']+=1
                walk(child,j,t,path+[(control,address[j],t)],free)
                b.del_and(control,address[j],t);free.append(t)
        for a,c,t in erased:
            free.remove(t);b.init_and(a,c,t);stats['recomputed_ands']+=1
    walk(0,-1,None,[],scratch.copy())
    return completed,stats

def build(delta,brot=10,return_report=False):
    raw=np.asarray(delta)
    if (raw.ndim!=2 or raw.shape[1]!=3 or raw.shape[0]<16 or
        raw.shape[0]%4 or raw.shape[0]>4096 or brot!=10 or
        not np.issubdtype(raw.dtype,np.integer) or
        np.any(raw<0) or np.any(raw>=1024)):
        raise ValueError('Expected 16..4096 address rows, a multiple of four, and three ten-bit integer words.')
    delta=raw.astype(np.int64);A,J=delta.shape
    if np.any(((delta>>9)^(delta>>8))&1):
        raise ValueError('Signed words require bit9=bit8.')
    W=10*J;BW=9*J;H=A//4;n=(A-1).bit_length()
    packed=[pack(row,10) for row in delta]
    ds=[[],[],[]];base=[]
    compact=lambda row:sum(((row>>(j*10))&511)<<(j*9) for j in range(J))
    for h in range(H):
        v=packed[4*h:4*h+4];base.append(v[0])
        ds[0].append(compact(v[0]^v[1]));ds[1].append(compact(v[0]^v[2]))
        ds[2].append(compact(v[0]^v[1]^v[2]^v[3]))
    def expand(v):
        return sum((((v>>(9*j))&511)|(((v>>(9*j+8))&1)<<9))<<(10*j) for j in range(J))
    restored=[base[h]^expand((ds[0][h] if a&1 else 0)^(ds[1][h] if a&2 else 0)^(ds[2][h] if (a&3)==3 else 0)) for h in range(H) for a in range(4)]
    assert restored==packed
    values=[base[h]|(ds[0][h]<<W)|(ds[1][h]<<(W+BW))|(ds[2][h]<<(W+2*BW)) for h in range(H)]
    b=km.CircuitBuilder();sel=b.create_quantum_register(n,'selection')
    tgt=[q for j in range(J) for q in b.create_quantum_register(10,f'word_{j}')]
    bank=list(b.alloc_qubits(3*BW));scratch=b.alloc_qubits(5)
    completed,stats=pebbled_init(b,values,sel[2:],tgt+bank,scratch);b.free(scratch)
    both=b.alloc_qubits(1)[0];b.init_and(sel[0],sel[1],both)
    for s,ctl in enumerate((sel[0],sel[1],both)):
        for j in range(J):
            for k in range(9):
                ts=[tgt[j*10+k]]+([tgt[j*10+9]] if k==8 else [])
                xor_and(b,ctl,bank[s*BW+j*9+k],ts)
    b.del_and(sel[0],sel[1],both);b.free(both)
    erase=del_lookup_anf(b,table=[v>>W for v in completed],address=sel[2:],target=bank,low_bits=min(5,n-3))
    b.free(km.array(bank));c=b.finish_circuit()
    report={'num_qubits':c.num_qubits,'static_toffoli':int(c.operation_statistics().num_toffoli_like),'table_identity':'PASS','init':stats,'cleanup':erase}
    assert c.num_qubits<=128
    return (c,report) if return_report else c
