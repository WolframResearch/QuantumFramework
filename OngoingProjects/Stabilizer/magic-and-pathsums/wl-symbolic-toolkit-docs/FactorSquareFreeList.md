# FactorSquareFreeList | [SpanFromLeft]

> [FactorSquareFreeList](https://reference.wolfram.com/language/ref/FactorSquareFreeList.html)[*poly*] — gives a list of square[Hyphen]free factors of a polynomial, together with their exponents.

## Details and Options

[FactorSquareFreeList](https://reference.wolfram.com/language/ref/FactorSquareFreeList.html) takes the following options:

| [Extension](https://reference.wolfram.com/language/ref/Extension.html) | [None](https://reference.wolfram.com/language/ref/None.html) | coefficient field to be used |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |
| [Trig](https://reference.wolfram.com/language/ref/Trig.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to do trigonometric as well as algebraic transformations  |

## Examples

### Basic Examples

List square-free factors of a polynomial:

```wolfram
FactorSquareFreeList[x^5-x^3-x^2+1]
(* Output *)
{{1,1},{-1+x,2},{1+2 x+2 x^2+x^3,1}}
```

List irreducible factors of a polynomial:

```wolfram
FactorList[x^5-x^3-x^2+1]
(* Output *)
{{1,1},{-1+x,2},{1+x,1},{1+x+x^2,1}}
```

### Scope

A univariate polynomial:

```wolfram
FactorSquareFreeList[x^4-9x^3+29x^2-39x+18]
(* Output *)
{{1,1},{-3+x,2},{2-3 x+x^2,1}}
```

A multivariate polynomial:

```wolfram
FactorSquareFreeList[x^5-x^3y^2-x^2y^3+y^5]
(* Output *)
{{1,1},{x-y,2},{x^3+2 x^2 y+2 x y^2+y^3,1}}
```

A rational function:

```wolfram
FactorSquareFreeList[(x^3+x^2)/(x^2-4y^2)-(x+1)/(x^2-4y^2)]
(* Output *)
{{1,1},{-1+x,1},{1+x,2},{x^2-4 y^2,-1}}
```

A polynomial with complex coefficients:

```wolfram
FactorSquareFreeList[x^4-2I x^3-2x^2+2 I x+1]
(* Output *)
{{1,1},{-ⅈ+x,2},{-1+x^2,1}}
```

Non-polynomial expressions:

```wolfram
FactorSquareFreeList[E^(3x)-3E^(2x)+3E^x-1]
(* Output *)
{{1,1},{-1+ℯ^x,3}}
```

```wolfram
FactorSquareFreeList[x+2Sqrt[x]+1]
(* Output *)
{{1,1},{1+Sqrt[x],2}}
```

Square-free factors of a polynomial over the integers modulo 3:

```wolfram
FactorSquareFreeList[x^6+1,Modulus ->3]
(* Output *)
{{1,1},{1+x^2,3}}
```

Square-free factors of a polynomial over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
FactorSquareFreeList[ℱ[1]x^4+ℱ[246]x^3+ℱ[4875]x^2+ℱ[4608]x+ℱ[304]]
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
(facs=FactorSquareFreeList[r]);//AbsoluteTiming
(* Output *)
{0.3672414,Null}
```

```wolfram
{Exponent[#[[1]],x],#[[2]]}&/@facs
(* Output *)
{{0,1},{2000,2},{2000,3}}
```

### Options

#### Extension

By default algebraic number coefficients are treated as independent variables:

```wolfram
FactorSquareFreeList[x^2+2Sqrt[2]x+2]
(* Output *)
{{1,1},{2+2 Sqrt[2] x+x^2,1}}
```

With [Extension](https://reference.wolfram.com/language/ref/Extension.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) algebraic dependencies between coefficients are recognized:

```wolfram
FactorSquareFreeList[x^2+2Sqrt[2]x+2,Extension->Automatic]
(* Output *)
{{1,1},{Sqrt[2]+x,2}}
```

Square-free factorization over a finite field:

```wolfram
ℱ=FiniteField[2,3];
```

```wolfram
FactorSquareFreeList[x^4+1,Extension->ℱ]
(* Output *)
![image](img/image_003.png)
```

#### Modulus

Find square-free factors over the integers modulo 2:

```wolfram
FactorSquareFreeList[1+x^4,Modulus->2]
(* Output *)
{{1,1},{1+x,4}}
```

#### Trig

List multiple factors in a trigonometric expression:

```wolfram
FactorSquareFreeList[1-Cos[2x],Trig->True]
(* Output *)
{{2,1},{Sin[x],2}}
```

### Applications

Remove root multiplicities:

```wolfram
f=x^9+9x^8+21x^7-27x^6-153x^5-81x^4+239x^3+207x^2-108x-108;
```

```wolfram
rtsf=x/.Solve[f==0,x]
(* Output *)
{-3,-3,-3,-2,-1,-1,1,1,2}
```

```wolfram
g=Times@@First/@FactorSquareFreeList[f]
(* Output *)
(3+x) (-4+x^2) (-1+x^2)
```

```wolfram
rtsg=x/.Solve[g==0,x]
(* Output *)
{-3,-2,-1,1,2}
```

```wolfram
rtsg===Union[rtsf]
(* Output *)
True
```

### Properties & Relations

[FactorSquareFreeList](https://reference.wolfram.com/language/ref/FactorSquareFreeList.html) gives a list of square-free factors:

```wolfram
f=x^9+9x^8+21x^7-27x^6-153x^5-81x^4+239x^3+207x^2-108x-108;
```

```wolfram
FactorSquareFreeList[f]
(* Output *)
{{1,1},{3+x,3},{-4+x^2,1},{-1+x^2,2}}
```

This multiplies the factors together:

```wolfram
Times@@Power@@@%
(* Output *)
(3+x)^3 (-4+x^2) (-1+x^2)^2
```

[FactorSquareFree](https://reference.wolfram.com/language/ref/FactorSquareFree.html) gives a product of factors:

```wolfram
FactorSquareFree[x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]
(* Output *)
(2+x) (3+x)^3 (-1+x^2)^2
```

[Expand](https://reference.wolfram.com/language/ref/Expand.html) combines all the factors back together:

```wolfram
Expand[%]
(* Output *)
54+81 x-63 x^2-151 x^3-35 x^4+59 x^5+43 x^6+11 x^7+x^8
```

[FactorList](https://reference.wolfram.com/language/ref/FactorList.html) gives a list of irreducible factors:

```wolfram
FactorList[f]
(* Output *)
{{1,1},{-2+x,1},{-1+x,2},{1+x,2},{2+x,1},{3+x,3}}
```

A univariate polynomial has multiple factors if and only if its [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) is zero:

```wolfram
Discriminant[x^5-x^3-x^2+1,x]
(* Output *)
0
```

```wolfram
FactorSquareFreeList[x^5-x^3-x^2+1]
(* Output *)
{{1,1},{-1+x,2},{1+2 x+2 x^2+x^3,1}}
```

```wolfram
Discriminant[x^5-x^3-x^2-1,x]
(* Output *)
7684
```

```wolfram
FactorSquareFreeList[x^5-x^3-x^2-1]
(* Output *)
{{1,1},{-1-x^2-x^3+x^5,1}}
```

[FactorInteger](https://reference.wolfram.com/language/ref/FactorInteger.html) gives a list of prime factors of an integer:

```wolfram
FactorInteger[1234567]
(* Output *)
{{127,1},{9721,1}}
```

```wolfram
Times@@Power@@@%
(* Output *)
1234567
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Factoring & Decomposition ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 2022 (13.2) ▪ 2023 (13.3)
