# Range | [SpanFromLeft]

> [Range](https://reference.wolfram.com/language/ref/Range.html)[*i*_*max*] — generates the list `{1,2,…,*i*_*max*}`.
> [Range](https://reference.wolfram.com/language/ref/Range.html)[*i*_*min*,*i*_*max*] — generates the list `{*i*_*min*,…,*i*_*max*}`.
> [Range](https://reference.wolfram.com/language/ref/Range.html)[*i*_*min*,*i*_*max*,*di*]  — uses step `*di*`.

## Details

The arguments to [Range](https://reference.wolfram.com/language/ref/Range.html) need not be integers.

[Range](https://reference.wolfram.com/language/ref/Range.html) starts from `*i*_*min*` and successively adds increments of `*di*` until the result is greater than `*i*_*max*`.

[Range](https://reference.wolfram.com/language/ref/Range.html) uses the standard Wolfram Language iteration specification, as applied to a single variable.

[Range](https://reference.wolfram.com/language/ref/Range.html) has attribute [Listable](https://reference.wolfram.com/language/ref/Listable.html).

## Examples

### Basic Examples

```wolfram
Range[4]
(* Output *)
{1,2,3,4}
```

```wolfram
Range[1.2,2.2,0.15]
(* Output *)
{1.2,1.3499999999999999,1.5,1.65,1.7999999999999998,1.95,2.1}
```

```wolfram
Range[x,x+4]
(* Output *)
{x,1+x,2+x,3+x,4+x}
```

### Scope

Use a step of 2:

```wolfram
Range[1,10,2]
(* Output *)
{1,3,5,7,9}
```

Use a negative step:

```wolfram
Range[10,1,-1]
(* Output *)
{10,9,8,7,6,5,4,3,2,1}
```

Use an exact numeric-valued step:

```wolfram
Range[0,10,π]
(* Output *)
{0,π,2 π,3 π}
```

Use a machine-number step:

```wolfram
Range[0,10,N[π]]
(* Output *)
{0.,3.141592653589793,6.283185307179586,9.42477796076938}
```

Use a precision-24 step:

```wolfram
Range[0,10, N[π, 24]]
(* Output *)
{0,3.1415926535897932384626433832795028842,6.2831853071795864769252867665590057684,9.42477796076937971538793014983850865259}
```

Range of very large numbers:

```wolfram
Range[2^225,2^225 + 5]
(* Output *)
{53919893334301279589334030174039261347274288845081144962207220498432,53919893334301279589334030174039261347274288845081144962207220498433,53919893334301279589334030174039261347274288845081144962207220498434,53919893334301279589334030174039261347274288845081144962207220498435,53919893334301279589334030174039261347274288845081144962207220498436,53919893334301279589334030174039261347274288845081144962207220498437}
```

### Generalizations & Extensions

Use a symbolic step:

```wolfram
Range[a,b,(b-a)/4]
(* Output *)
{a,a+(1)/(4) (-a+b),a+(1)/(2) (-a+b),a+(3)/(4) (-a+b),b}
```

Use a list of range specifications:

```wolfram
Range[{5,2,6,3}]
(* Output *)
{{1,2,3,4,5},{1,2},{1,2,3,4,5,6},{1,2,3}}
```

### Applications

Produce a geometric sequence:

```wolfram
q^Range[5]
(* Output *)
{q,q^2,q^3,q^4,q^5}
```

Form a polynomial from coefficients:

```wolfram
coeff = {-2,9,5,3,-3,-6,-7,-4,8,3};
```

```wolfram
poly = coeff.x^Range[0,Length[coeff]-1]
(* Output *)
-2+9 x+5 x^2+3 x^3-3 x^4-6 x^5-7 x^6-4 x^7+8 x^8+3 x^9
```

Form a random permutation:

```wolfram
RandomSample[Range[10]]
(* Output *)
{10,6,8,1,9,7,3,5,4,2}
```

Find an inverse permutation:

```wolfram
perm = RandomSample[Range[10]]
(* Output *)
{4,2,8,1,10,5,3,6,9,7}
```

```wolfram
inverse = perm;
inverse[[perm]]=Range[Length[perm]];
inverse
(* Output *)
{4,2,7,1,6,8,10,3,9,5}
```

### Properties & Relations

[Range](https://reference.wolfram.com/language/ref/Range.html)[*i*_*min*,*i*_*max*,*di*] is equivalent to [Table](https://reference.wolfram.com/language/ref/Table.html)[*i*,{*i*_*min*,*i*_*max*,*di*}]:

```wolfram
Range[-4,9,3]
(* Output *)
{-4,-1,2,5,8}
```

```wolfram
Table[i,{i,-4,9,3}]
(* Output *)
{-4,-1,2,5,8}
```

Use [Range](https://reference.wolfram.com/language/ref/Range.html) or [Span](https://reference.wolfram.com/language/ref/Span.html) (;;) as [Part](https://reference.wolfram.com/language/ref/Part.html) specification:

```wolfram
list={a,b,c,d,e};
```

```wolfram
list[[Range[1,5,2]]]
(* Output *)
{a,c,e}
```

```wolfram
list[[1;;5;;2]]
(* Output *)
{a,c,e}
```

### Possible Issues

For some step sizes, [Range](https://reference.wolfram.com/language/ref/Range.html) may not include the upper limit given:

```wolfram
Range[0,10,3]
(* Output *)
{0,3,6,9}
```

Even though the lower limit was exact, the inexact step makes the first element inexact:

```wolfram
Range[0,1,.1]
(* Output *)
{0.,0.1,0.2,0.30000000000000004,0.4,0.5,0.6000000000000001,0.7000000000000001,0.8,0.9,1.}
```

[Range](https://reference.wolfram.com/language/ref/Range.html)  accepts [Quantity](https://reference.wolfram.com/language/ref/Quantity.html) expressions as limits and steps:

```wolfram
Range[1,9,2]
(* Output *)
{1,3,5,7,9}
```

For [Quantity](https://reference.wolfram.com/language/ref/Quantity.html) expressions, [Precision](https://reference.wolfram.com/language/ref/Precision.html) is taken into account when determining whether elements are within the bounds of the limits:

```wolfram
Precision[UnitConvert[1]]
(* Output *)
9.181924092215384
```

```wolfram
Range[40]
(* Output *)
{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39}
```

```wolfram
Length@%
(* Output *)
39
```

### Neat Examples

Make nested ranges:

```wolfram
Range[Range[5]]
(* Output *)
{{1},{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}}
```

```wolfram
Range[Range[Range[3]]]
(* Output *)
{{{1}},{{1},{1,2}},{{1},{1,2},{1,2,3}}}
```

```wolfram
Nest[Range,3,6]
(* Output *)
{{{{{{1}}}}},{{{{{1}}}},{{{{1}}},{{{1}},{{1},{1,2}}}}},{{{{{1}}}},{{{{1}}},{{{1}},{{1},{1,2}}}},{{{{1}}},{{{1}},{{1},{1,2}}},{{{1}},{{1},{1,2}},{{1},{1,2},{1,2,3}}}}}}
```

Show it in tree form:

```wolfram
TreeForm[%]
(* Output *)
![image](img/image_001.png)
```

## Tech Notes ▪Vectors and Matrices ▪Making Tables of Values

## Related Guides ▪Constructing Lists ▪Regular & Coordinate Arrays ▪List Manipulation ▪Incrementals ▪Graphics Annotation & Appearance

## Related Links [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/03-first-look-at-lists.html)[: First Look at Lists](https://www.wolfram.com/language/elementary-introduction/03-first-look-at-lists.html) [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/20-options.html)[: Options](https://www.wolfram.com/language/elementary-introduction/20-options.html) [NKS|Online](http://www.wolframscience.com/nks/search/?q=Range) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0)
