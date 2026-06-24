# Exponent | [SpanFromLeft]

> [Exponent](https://reference.wolfram.com/language/ref/Exponent.html)[*expr*,*form*] — gives the maximum power with which `*form*` appears in the expanded form of `*expr*`.
> [Exponent](https://reference.wolfram.com/language/ref/Exponent.html)[*expr*,*form*,*h*] — applies `*h*` to the set of exponents with which `*form*` appears in `*expr*`.

## Details and Options

The default taken for `*h*` is [Max](https://reference.wolfram.com/language/ref/Max.html).

`*form*` can be a product of terms.

[Exponent](https://reference.wolfram.com/language/ref/Exponent.html) works whether or not `*expr*` is explicitly given in expanded form.

[Exponent](https://reference.wolfram.com/language/ref/Exponent.html)[0,*x*] is `-[Infinity](https://reference.wolfram.com/language/ref/Infinity.html)`.

[Exponent](https://reference.wolfram.com/language/ref/Exponent.html)[*expr*,{*form*_1,*form*_2,…}] gives the list of exponents for each of the `*form*_*i*`.

[Exponent](https://reference.wolfram.com/language/ref/Exponent.html) takes the following options:

| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| --- | --- | --- |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations |

## Examples

### Basic Examples

Find the highest exponent of $x$:

```wolfram
Exponent[1+x^2+a x^3,x]
(* Output *)
3
```

### Scope

The degree of a polynomial:

```wolfram
Exponent[(x^2+1)^3+1,x]
(* Output *)
6
```

Exponents may be rational numbers or symbolic expressions:

```wolfram
Exponent[x^(n+1)+2 Sqrt[x]+1,x]
(* Output *)
Max[(1)/(2),1+n]
```

The lowest exponent in a polynomial:

```wolfram
Exponent[(x^2+1)^3-1,x,Min]
(* Output *)
2
```

The list of all exponents with which $x$ appears:

```wolfram
Exponent[1+x^2+a x^3,x,List]
(* Output *)
{0,2,3}
```

### Options

#### Modulus

The degree of a polynomial over the integers modulo 2:

```wolfram
Exponent[2x^2-x+1,x,Modulus->2]
(* Output *)
1
```

#### Trig

With [Trig](https://reference.wolfram.com/language/ref/Trig.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Exponent](https://reference.wolfram.com/language/ref/Exponent.html) recognizes dependencies between trigonometric functions:

```wolfram
Exponent[Tan[x],Cos[x],Trig->True]
(* Output *)
-1
```

### Applications

Compute the leading coefficient:

```wolfram
LeadingCoefficient[poly_,x_]:=Coefficient[poly,x,Exponent[poly,x]]
```

```wolfram
LeadingCoefficient[2+3x+17x^5,x]
(* Output *)
17
```

Compute the leading term:

```wolfram
LeadingTerm[poly_,x_]:=
With[{n=Exponent[poly,x]},Coefficient[poly,x,n]x^n]
```

```wolfram
LeadingTerm[2+3x+17x^5,x]
(* Output *)
17 x^5
```

### Properties & Relations

The number of complex roots of a polynomial is equal to its degree:

```wolfram
f=(x+1)^5-2x+3;
```

```wolfram
Exponent[f,x]
(* Output *)
5
```

Use [Solve](https://reference.wolfram.com/language/ref/Solve.html) to find the roots:

```wolfram
Length[x/.Solve[f==0,x]]
(* Output *)
5
```

Length of the [CoefficientList](https://reference.wolfram.com/language/ref/CoefficientList.html) of a polynomial is one more than its degree:

```wolfram
f=(x^2+2x-1)^7-3;
```

```wolfram
Exponent[f,x]
(* Output *)
14
```

```wolfram
Length[CoefficientList[f,x]]
(* Output *)
15
```

### Possible Issues

[Exponent](https://reference.wolfram.com/language/ref/Exponent.html) is purely syntactical; it does not attempt to recognize zero coefficients:

```wolfram
zero=Sqrt[2]+Sqrt[3]-Sqrt[5+2Sqrt[6]];
f=zero x^2+x+1;
```

```wolfram
Exponent[f,x]
(* Output *)
2
```

```wolfram
Exponent[RootReduce[f],x]
(* Output *)
1
```

## Tech Notes ▪Picking Out Pieces of Algebraic Expressions ▪Finding the Structure of a Polynomial

## Related Guides ▪Polynomial Algebra ▪Formula Manipulation

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2003 (5.0)
