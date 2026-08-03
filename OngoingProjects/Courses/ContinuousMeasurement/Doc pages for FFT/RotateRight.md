# RotateRight | [SpanFromLeft]

> [RotateRight](https://reference.wolfram.com/language/ref/RotateRight.html)[*expr*,*n*] — cycles the elements in `*expr*` `*n*` positions to the right.
> [RotateRight](https://reference.wolfram.com/language/ref/RotateRight.html)[*expr*] — cycles one position to the right.
> [RotateRight](https://reference.wolfram.com/language/ref/RotateRight.html)[*expr*,{*n*_1,*n*_2,…}] — cycles elements at successive levels `*n*_*i*` positions to the right.

## Details

[RotateRight](https://reference.wolfram.com/language/ref/RotateRight.html)[*expr*,-*n*] rotates `*n*` positions to the left.

[RotateRight](https://reference.wolfram.com/language/ref/RotateRight.html) can be used on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) and [Association](https://reference.wolfram.com/language/ref/Association.html) objects.

## Examples

### Basic Examples

Rotate two positions to the right:

```wolfram
RotateRight[{a,b,c,d,e},2]
(* Output *)
{d,e,a,b,c}
```

Rotate one position to the right:

```wolfram
RotateRight[{a,b,c,d,e}]
(* Output *)
{e,a,b,c,d}
```

Rotate [Association](https://reference.wolfram.com/language/ref/Association.html) one position to the right:

```wolfram
RotateRight[<|1->a,2->b,3->c|>]
(* Output *)
<|3->c,1->a,2->b|>
```

Rotate [Association](https://reference.wolfram.com/language/ref/Association.html) on the first and second levels:

```wolfram
RotateRight[<|1->{a,b},2->{b,c},3->{c,d}|>,{1,1}]
(* Output *)
<|3->{d,c},1->{b,a},2->{c,b}|>
```

### Generalizations & Extensions

Rotate one position right at the first level, and left at the second level:

```wolfram
RotateRight[{{a,b,c},{d,e,f},{g,h,i}},{1,-1}]
(* Output *)
{{h,i,g},{b,c,a},{e,f,d}}
```

Rotate an expression with any head:

```wolfram
RotateRight[f[x,y,z]]
(* Output *)
f[z,x,y]
```

### Applications

Successively rotate a list right:

```wolfram
NestList[RotateRight,{a,b,c,d,e},4]
(* Output *)
{{a,b,c,d,e},{e,a,b,c,d},{d,e,a,b,c},{c,d,e,a,b},{b,c,d,e,a}}
```

Rotate successive rows of a matrix by their row number:

```wolfram
MapIndexed[RotateRight,Table[{a,b,c,d},{5}]]//TableForm
(* Output *)
{{d, a, b, c}, {c, d, a, b}, {b, c, d, a}, {a, b, c, d}, {d, a, b, c}}
```

Rotate a 2D image:

```wolfram
ArrayPlot[RotateRight[Table[If[x^2+y^2<1,1,0],{x,-1,1,1/100},{y,-1,1,1/100}],{40,20}]]
```

*([Graphics])*

Rotate operands:

```wolfram
RotateRight[a.b.c.d.e]
(* Output *)
e.a.b.c.d
```

## Tech Notes ▪Rearranging Lists ▪Rearranging Nested Lists

## Related Guides ▪Rearranging & Restructuring Lists ▪Parts of Matrices

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0)
