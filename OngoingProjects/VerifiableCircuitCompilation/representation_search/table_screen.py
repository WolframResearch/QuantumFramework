"""Exact fixed-table representation screen. No circuits or factor changes.

Uses the published C30 generator's signed words (hash pinned), and only emits
aggregate support and rank information. Raw tables remain under local/.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def mobius(table):
    result = list(table)
    for bit in range((len(result)-1).bit_length()):
        for mask in range(len(result)):
            if mask & (1 << bit):
                result[mask] ^= result[mask ^ (1 << bit)]
    return result


def rank2(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def pack(row):
    return sum(int(value) << (10*j) for j, value in enumerate(row))


def compact(word):
    return sum(((word >> (10*j)) & 511) << (9*j) for j in range(3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--words', type=Path, required=True)
    args = parser.parse_args()
    word_hash = hashlib.sha256(args.words.read_bytes()).hexdigest()
    if word_hash != '8903db9b02b1c314cb3b6c8b760738791e42d4bbce7a7738170217b95ac84761':
        raise ValueError('Expected the verified C30 signed-word table.')
    words = np.fromfile(args.words, dtype='<u2').reshape(3976, 75).astype(np.int64)
    batches = words.reshape(3976,25,3).transpose(1,0,2).copy()
    batches[1:] ^= batches[:-1].copy()
    records = []
    for batch, delta in enumerate(batches):
        packed = [pack(row) for row in delta]
        values = []
        for start in range(0,len(packed),4):
            a,b,c,d = packed[start:start+4]
            values.append(a | compact(a^b)<<30 | compact(a^c)<<57 | compact(a^b^c^d)<<84)
        coefficients = mobius(values + [0]*30)
        coefficients[994:] = [0]*30
        if mobius(coefficients)[:994] != values:
            raise ArithmeticError('Prefix completion identity failed.')
        histogram = {str(degree):sum(bool(value) and mask.bit_count()==degree
                      for mask,value in enumerate(coefficients)) for degree in range(11)}
        records.append(dict(batch=batch,output_width=111,legal_rows=994,
            nonlinear_support=sum(int(histogram[str(k)]) for k in range(2,11)),
            support_by_degree=histogram,all_output_rank=rank2(values),
            difference_bank_rank=rank2([word>>30 for word in values])))
    report = dict(scope='EXACT_FIXED_TABLE_REPRESENTATION_SCREEN_NO_CIRCUIT',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        words_sha256=word_hash,batches=records,
        conclusion='Invalid suffix completion and invertible output basis changes cannot remove any legal-prefix ANF support. This is not a bound on general AND/XOR circuits.')
    (Path(__file__).parent/'table-screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:sorted({r[key] for r in records}) for key in
                      ['nonlinear_support','all_output_rank','difference_bank_rank']}))


if __name__ == '__main__':
    main()
