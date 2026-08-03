# RotateLeft | [SpanFromLeft]

> [RotateLeft](https://reference.wolfram.com/language/ref/RotateLeft.html)[*expr*,*n*] — cycles the elements in `*expr*` `*n*` positions to the left.
> [RotateLeft](https://reference.wolfram.com/language/ref/RotateLeft.html)[*expr*] — cycles one position to the left.
> [RotateLeft](https://reference.wolfram.com/language/ref/RotateLeft.html)[*expr*,{*n*_1,*n*_2,…}] — cycles elements at successive levels `*n*_*i*` positions to the left.

## Details

[RotateLeft](https://reference.wolfram.com/language/ref/RotateLeft.html)[*expr*,-*n*] rotates `*n*` positions to the right.

[RotateLeft](https://reference.wolfram.com/language/ref/RotateLeft.html) can be used on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) and [Association](https://reference.wolfram.com/language/ref/Association.html) objects.

## Examples

### Basic Examples

Rotate two positions to the left:

```wolfram
RotateLeft[{a,b,c,d,e},2]
(* Output *)
{c,d,e,a,b}
```

Rotate one position to the left:

```wolfram
RotateLeft[{a,b,c,d,e}]
(* Output *)
{b,c,d,e,a}
```

Rotate [Association](https://reference.wolfram.com/language/ref/Association.html) one position to the left:

```wolfram
RotateLeft[<|1->a,2->b,3->c|>]
(* Output *)
<|2->b,3->c,1->a|>
```

Rotate [Association](https://reference.wolfram.com/language/ref/Association.html) on the first and second levels:

```wolfram
RotateRight[<|1->{a,b},2->{b,c},3->{c,d}|>,{1,1}]
(* Output *)
<|3->{d,c},1->{b,a},2->{c,b}|>
```

### Generalizations & Extensions

Rotate one position left at the first level, and right at the second level:

```wolfram
RotateLeft[{{a,b,c},{d,e,f},{g,h,i}},{1,-1}]
(* Output *)
{{f,d,e},{i,g,h},{c,a,b}}
```

Rotate an expression with any head:

```wolfram
RotateLeft[f[x,y,z]]
(* Output *)
f[y,z,x]
```

### Applications

Successively rotate a list left:

```wolfram
NestList[RotateLeft,{a,b,c,d,e},4]
(* Output *)
{{a,b,c,d,e},{b,c,d,e,a},{c,d,e,a,b},{d,e,a,b,c},{e,a,b,c,d}}
```

Rotate successive rows of a matrix by their row number:

```wolfram
MapIndexed[RotateLeft,Table[{a,b,c,d},{5}]]//TableForm
(* Output *)
{{b, c, d, a}, {c, d, a, b}, {d, a, b, c}, {a, b, c, d}, {b, c, d, a}}
```

Rotate a 2D image:

```wolfram
ArrayPlot[RotateLeft[Table[If[x^2+y^2<1,1,0],{x,-1,1,1/100},{y,-1,1,1/100}],{40,20}]]
```

*([Graphics])*

Rotate operands:

```wolfram
RotateLeft[a.b.c.d.e]
(* Output *)
b.c.d.e.a
```

## Tech Notes ▪Rearranging Lists ▪Rearranging Nested Lists

## Related Guides ▪Rearranging & Restructuring Lists ▪Parts of Matrices

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=RotateLeft) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0)
