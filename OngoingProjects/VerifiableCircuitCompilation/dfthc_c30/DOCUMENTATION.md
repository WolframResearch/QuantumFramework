# Loader mathematics and WL contract

The quantum operation is a reversible classical XOR oracle on the address
and angle-word registers, with measurement-assisted scratch cleanup. It is
used inside the DFTHC block encoding; no 113-electron state is simulated here.

For four packed rows f0,f1,f2,f3, define GF(2) differences d1=f0+f1,
d2=f0+f2 and d3=f0+f1+f2+f3. For low address bits a,b, the lookup polynomial is
`f0 + a d1 + b d2 + a b d3`. Substitution of the four Boolean addresses
gives each original row exactly. WL keeps the coefficients symbolic and
reduces only modulo 2, before testing any released numbers.

The signed ten-bit interval is [-256,255] modulo 1024. Bits 9 and 8 agree
on all 512 values. Nine compact bits therefore suffice for every difference;
expansion duplicates bit 8 into bit 9. Three words use 81 bank wires. The
FeMoCo allocation is 12 address + 30 target + 81 bank + 5 scratch = 128.
The one-wire low-bit AND and measurement cleanup reuse those scratch wires.
Q128 is the ANGLE limit; the full block encoding has additional persistent
registers and shares scratch between its stages.

The Givens gauge decodes the public words, eliminates leaves in reverse tree
order and rounds back to the signed ten-bit grid. The supplied branch-zero
threshold is 1e-12; it is retained from the released implementation. This
numerical gauge is accepted only if decoded vectors differ from the public
ones by row signs to within 1e-11 and pass the official metric gates. WL proves
`(-u)(-u)^T = u u^T`; the generator checks the actual vectors and norms.
The gauge does not optimize weights, offsets, shifts, tree or rank.

Scratch holding an integer Boolean lookup F(a) can be reset in the X basis.
A measured bit m contributes phase `(-1)^(m F(a))`; applying the same phase
correction cancels it since `(-1)^(2 m F(a))=1`. The split-ANF eraser evaluates
this phase polynomial from the measured target bits and address monomials.
Its low and high products use measurement-based AND cleanup. The fanout
sandwich around CCX XORs the same AND into every target even when those
targets start dirty. WL checks the phase identity for arbitrary integers,
ANF linearity, fanout cancellation, and XOR involution exactly. Circuit
simulation checks the emitted implementation, including phases and scratch.

## Documented WL symbols

All symbols are core WL. Installed Wolfram 15.0 reference pages for
[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html)
and [BinaryReadList](https://reference.wolfram.com/language/ref/BinaryReadList.html)
were converted with nb-reader and read. Integer PolynomialMod reduces
coefficients modulo 2; polynomial moduli also impose Boolean idempotence.
BinaryReadList opens and closes a file passed by name, and uses
`"UnsignedInteger16", ByteOrdering->-1` for the generator's little-endian
word cache. Live usage and execution confirm that BitGet is zero-indexed and
the bit operations thread over the arrays. No QuantumFramework API is used.

FourAddress accepts a nonempty four-row matrix of signed ten-bit integers;
other inputs return Failure. The exact checks include both interval endpoints,
a nonconstant three-word table and malformed inputs. Full table checks
reconstruct every released row and validate batch accumulation and reversal.

## Required regimes

| Regime | Executed evidence |
|---|---|
| General symbolic | Arbitrary f0,f1,f2,f3 GF(2) interpolation; integer measurement-phase cancellation; symbolic row-sign projector and Givens norm. |
| Exactly solvable | All 512 signed words; a nonconstant three-word example; XOR involution and dirty-target fanout. |
| Limiting | Zero and constant lookup tables, so differences vanish; tested in WL and emitted circuits at four address sizes. |
| Numerical reference | All 298200 angle words; fixed-factor metric gates; full official 3976-address, 64-shot, 25-circuit forward/reverse verification. |
| Failure and edge | Invalid shape/type/word/batch rejection; signed endpoints; Q128; randomized dirty targets; phase, selection and scratch checks. |

The combined local estimate is the geometric mean of the baseline-to-current
cost ratios across inverse square root, DQI and DFTHC. It is not a reported
Kaggle leaderboard score. The public factors' two-body error is close to the
20 mHa gate; the exact factors and words are hash-bound to the evidence.
