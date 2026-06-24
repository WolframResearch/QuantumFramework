# FactorSquareFree | [SpanFromLeft]

> [FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html)[*poly*] — pulls out any multiple factors in a polynomial.

## Details and Options

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | coefficient field to be used |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations  |

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html)[*poly*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] extends the coefficient field to include algebraic numbers that appear in the coefficients of `*poly*`.

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) automatically threads over lists, as well as equations, inequalities and logic functions.

## Examples

### Basic Examples

Pull out multiple factors:

```wolfram
FactorSquareFree[x^5-x^3-x^2+1]
(* Output *)
(-1+x)^2 (1+2 x+2 x^2+x^3)
```

A complete factorization:

```wolfram
Factor[x^5-x^3-x^2+1]
(* Output *)
(-1+x)^2 (1+x) (1+x+x^2)
```

### Scope

A univariate polynomial:

```wolfram
FactorSquareFree[x^4-9x^3+29x^2-39x+18]
(* Output *)
(-3+x)^2 (2-3 x+x^2)
```

A multivariate polynomial:

```wolfram
FactorSquareFree[x^5-x^3y^2-x^2y^3+y^5]
(* Output *)
(x-y)^2 (x^3+2 x^2 y+2 x y^2+y^3)
```

A rational function:

```wolfram
FactorSquareFree[(x^3+x^2)/(x^2-4y^2)-(x+1)/(x^2-4y^2)]
(* Output *)
((-1+x) (1+x)^2)/(x^2-4 y^2)
```

A polynomial with complex coefficients:

```wolfram
FactorSquareFree[x^4-2I x^3-2x^2+2 I x+1]
(* Output *)
(-ⅈ+x)^2 (-1+x^2)
```

Non-polynomial expressions:

```wolfram
FactorSquareFree[E^(3x)-3E^(2x)+3E^x-1]
(* Output *)
(-1+ℯ^x)^3
```

```wolfram
FactorSquareFree[x+2Sqrt[x]+1]
(* Output *)
(1+Sqrt[x])^2
```

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) threads over lists:

```wolfram
FactorSquareFree[{(x^2-1)(x-1),(x^4-1)(x^2-1)}]
(* Output *)
{(-1+x)^2 (1+x),(-1+x^2)^2 (1+x^2)}
```

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) threads over equations and inequalities:

```wolfram
FactorSquareFree[1<x^4+2x^2+1<2]
(* Output *)
1<(1+x^2)^2<2
```

Square-free factorization of a polynomial over the integers modulo 3:

```wolfram
FactorSquareFree[x^6+1,Modulus ->3]
(* Output *)
(1+x^2)^3
```

Square-free factorization of a polynomial over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
FactorSquareFree[ℱ[1]x^4+ℱ[246]x^3+ℱ[4875]x^2+ℱ[4608]x+ℱ[304]]
(* Output *)
![image](img/image_001.png)
```

Compute the square-free factorization of a polynomial of degree $10000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[2000]; q=rpoly[2000];
r=Expand[p^2 q^3];
```

```wolfram
FactorSquareFree[r]//Short//AbsoluteTiming
(* Output *)
{0.2815757,-(-912-892 x+<<2976>>+695 x^1999+86 x^2000)^2 (460+<<2982>>+222 x^<<4>>)^3}
```

### Options

#### Extension

By default, algebraic number coefficients are treated as independent variables:

```wolfram
FactorSquareFree[x^2+2Sqrt[2]x+2]
(* Output *)
2+2 Sqrt[2] x+x^2
```

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html), algebraic dependencies between coefficients are recognized:

```wolfram
FactorSquareFree[x^2+2Sqrt[2]x+2,Extension->Automatic]
(* Output *)
(Sqrt[2]+x)^2
```

Square-free factorization over a finite field:

```wolfram
ℱ=FiniteField[2,3];
```

```wolfram
FactorSquareFree[x^4+1,Extension->ℱ]
(* Output *)
![image](img/image_003.png)
```

#### Modulus

Pull out multiple factors over the integers modulo 2:

```wolfram
FactorSquareFree[1+x^4,Modulus->2]
(* Output *)
(1+x)^4
```

#### Trig

Pull out multiple factors in a trigonometric expression:

```wolfram
FactorSquareFree[1-Cos[2x],Trig->True]
(* Output *)
2 Sin[x]^2
```

### Properties & Relations

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) only pulls out multiple factors:

```wolfram
f=x^9+9x^8+21x^7-27x^6-153x^5-81x^4+239x^3+207x^2-108x-108;
```

```wolfram
FactorSquareFree[f]
(* Output *)
(3+x)^3 (-4+x^2) (-1+x^2)^2
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) gives a complete factorization:

```wolfram
Factor[f]
(* Output *)
(-2+x) (-1+x)^2 (1+x)^2 (2+x) (3+x)^3
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) is effectively the inverse of [FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html):

```wolfram
FactorSquareFree[x^5-x^3-x^2+1]
(* Output *)
(-1+x)^2 (1+2 x+2 x^2+x^3)
```

```wolfram
Expand[%]
(* Output *)
1-x^2-x^3+x^5
```

[FactorSquareFreeList](https://reference.wolfram.com/language/ref/FactorSquareFreeList.html) gives a list of factors:

```wolfram
FactorSquareFreeList[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
{{1,1},{2+x,1},{3+x,3},{-1+x^2,2}}
```

A univariate polynomial has multiple factors if and only if its [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) is zero:

```wolfram
Discriminant[x^5-x^3-x^2+1,x]
(* Output *)
0
```

```wolfram
FactorSquareFree[x^5-x^3-x^2+1]
(* Output *)
(-1+x)^2 (1+2 x+2 x^2+x^3)
```

```wolfram
Discriminant[x^5-x^3-x^2-1,x]
(* Output *)
7684
```

```wolfram
FactorSquareFree[x^5-x^3-x^2-1]
(* Output *)
-1-x^2-x^3+x^5
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Factoring & Decomposition ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2007 (6.0) ▪ 2022 (13.2) ▪ 2023 (13.3)
