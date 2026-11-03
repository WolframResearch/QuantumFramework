import re, sys
sys.path.insert(0, '.')
from spec import S
base = __import__("os").path.normpath(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", ".."))
src = open(base + '/C-section.md').read().split('\n')
# entries
starts = [i for i, l in enumerate(src) if l.startswith('### ')]
nonf = next(i for i, l in enumerate(src) if l.startswith('## Non-formula records'))
ends = starts[1:] + [nonf]
# cell outputs
outs = {}
cur = None
for l in open('cells_out.txt').read().split('\n'):
    if l.startswith('=== '):
        cur = l[4:].strip()
    elif cur and l.strip():
        outs[cur] = l
def short(s, n=420):
    return s if len(s) <= n else s[:230] + ' ... ' + s[-150:]
over = {
 'C048': "{E^((GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2]*Tanh[r])/2)**Sech[r]^(1/2 + SuperDagger[\\[FormalA]1]**\\[FormalA]1)**E^(-1/2*(GeneralizedPower[NonCommutativeMultiply, \\[FormalA]1, 2]*Tanh[r])), <one-mode complex xi form with Tan[Sqrt[-(xi*Conjugate[xi])]]>, E^((xi*SuperDagger[\\[FormalA]1]**SuperDagger[\\[FormalA]2]*Tan[Sqrt[-(xi*Conjugate[xi])]])/Sqrt[-(xi*Conjugate[xi])])**Sec[Sqrt[-(xi*Conjugate[xi])]]^(SuperDagger[\\[FormalA]1]**\\[FormalA]1 + SuperDagger[\\[FormalA]2]**\\[FormalA]2)**E^((Sqrt[-(xi*Conjugate[xi])]*\\[FormalA]1**\\[FormalA]2*Tan[Sqrt[-(xi*Conjugate[xi])]])/xi)*Sec[Sqrt[-(xi*Conjugate[xi])]], {Tanh[r]/r, Sech[r]}}",
}
out = []
out.append('# C section: Combescure and Robert, Coherent States and Applications in Mathematical Physics (atlas lines 2981 to 4360)')
out.append('')
out.append('## Errors found')
out.append('')
out.append("No new error in the atlas text. The source misprints the atlas flags and corrects (C-018 missing 1/2, C-026 proof signs, C-028 amplitude, C-030 missing 2^n, C-047 sinh for tanh, C-050 endpoint, C-053 power, C-059 overlap parenthesis, C-067 sign, C-069 block-to-paired sign, C-076 argument sign, C-082 inverse label, C-085 Cayley sign, C-086 exponent sign and (-1)^n, C-091 Fresnel phase, C-094 global claim, C-098 Fourier of 1) were re-derived, and in each case the corrected form is the one that checks. Three corrections were tested against the uncorrected print: C-050 (the printed endpoint fails by an order-one residual in the 2-by-2 representation), C-085 (the one-mode hand check gives the minus sign) and C-098 (the exterior-algebra transform of 1 is the ordered top product, not alpha^* . alpha).")
out.append('')
out.append("Convention that affects QF users of this section: QF's `SqueezeOperator[xi]` is exp((conj(xi) a^2 - xi a^dag^2)/2), the atlas operator S(-xi) of C-025, C-027, C-048.")
out.append('')
out.append('Every QF cell below runs after these loading lines, each cell in a fresh session (scripts: `drafts/exact/C-work/cells/*.wl`, runner `runcells.wls`; non-QF exact checks: `exactchecks.wl`; outputs `cells_out.txt`, `exact_out.txt`):')
out.append('')
out.append('```wl')
out.append('')
out.append('Needs["Wolfram`QuantumFramework`"];')
out.append('Needs["Wolfram`QuantumFramework`SecondQuantization`"];')
out.append('```')
out.append('')
assert len(starts) == 90
for k, (s0, s1) in enumerate(zip(starts, ends)):
    block = src[s0:s1]
    title, topic, appears = block[0], block[1], block[2]
    assert topic.startswith('Topic:') and appears.startswith('Appears in:'), (title,)
    ai = next(i for i, l in enumerate(block) if l.startswith('Audit:'))
    seg = block[3:ai]
    while seg and not seg[0].strip(): seg = seg[1:]
    while seg and not seg[-1].strip(): seg = seg[:-1]
    sent, audit, cell, note = S[k]
    out += [title, topic, appears, '', sent, '']
    out += seg
    if not seg[-1].startswith('Conditions:'):
        out.append('')
    out.append('Audit: ' + audit)
    if cell:
        cnum = 'C-' + cell[1:]
        assert cnum in appears, (title, cell)
        code = open('cells/' + cell + '.wl').read().rstrip('\n')
        o = over.get(cell, short(outs[cell]))
        assert '`' not in o
        out.append('QF code (exact):')
        out.append('```wl')
        out.append(code)
        out.append('```')
        out.append('Output: `' + o + '`. ' + note)
    else:
        out.append('QF code (exact): none. ' + note)
    out.append('')
out += src[nonf:]
txt = '\n'.join(out).rstrip('\n') + '\n'
assert '—' not in txt and '–' not in txt
open(base + '/exact/C-section.md', 'w').write(txt)
print(len(starts), txt.count('QF code (exact):\n'), txt.count('QF code (exact): none'))
