# PolynomialQuotientRemainder | [SpanFromLeft]

> [PolynomialQuotientRemainder](https://reference.wolfram.com/language/ref/PolynomialQuotientRemainder.html)[*p*,*q*,*x*]  —  gives a list of the quotient and remainder of `*p*` and `*q*`, treated as polynomials in `*x*`.

## Details and Options

The remainder will always have a degree not greater than `*q*`.

## Examples

### Basic Examples

Find the quotient and remainder after dividing one polynomial by another:

```wolfram
PolynomialQuotientRemainder[x^4+2x+1,x^2+1,x]
(* Output *)
{-1+x^2,2+2 x}
```

The dividend is equal to the product of the quotient and the divisor plus the remainder:

```wolfram
x^4+2x+1==%[[1]](x^2+1)+%[[2]]//Expand
(* Output *)
True
```

Find the quotient and remainder for polynomials with symbolic coefficients:

```wolfram
PolynomialQuotientRemainder[x^3,a x+b,x]
(* Output *)
{(b^2)/(a^3)-(b x)/(a^2)+(x^2)/(a),-(b^3)/(a^3)}
```

### Scope

The resulting polynomials will have coefficients that are rational expressions of input coefficients:

```wolfram
PolynomialQuotientRemainder[x^2+x+1,2x+1,x]
(* Output *)
{(1)/(4)+(x)/(2),(3)/(4)}
```

```wolfram
PolynomialQuotientRemainder[x^2+ b x+1,a x+1,x]
(* Output *)
{-(1)/(a^2)+(b)/(a)+(x)/(a),1+(1)/(a^2)-(b)/(a)}
```

```wolfram
PolynomialQuotientRemainder[x^2+x+1,Pi x+1,x]
(* Output *)
{-(1)/(π^2)+(1)/(π)+(x)/(π),1+(1)/(π^2)-(1)/(π)}
```

Polynomial quotient and remainder over the integers modulo $5$:

```wolfram
PolynomialQuotientRemainder[x^2+4 x+1,2x+1,x,Modulus->5]
(* Output *)
{3+3 x,3}
```

Polynomial quotient and remainder over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialQuotientRemainder[ℱ[1]x^5+ℱ[123]x+ℱ[456],ℱ[789]x^3+ℱ[987]x+ℱ[654],x]
(* Output *)
![image](img/image_001.png)
```

[PolynomialQuotientRemainder](https://reference.wolfram.com/language/ref/PolynomialQuotientRemainder.html) also works for rational functions:

```wolfram
{quo,rem}=PolynomialQuotientRemainder[(x^2+1)/((x+2)(x-1)),(x^2-1)/(x+3),x]
(* Output *)
{(-9+x^2)/((-1+x) (-2+x+x^2)),(2)/(-1+x)}
```

The quotient and remainder of division of $\frac{f}{h}$ by $\frac{g}{k}$ are $\frac{q k}{d h}$ and $\frac{r}{d}$, where $d=gcd(g,h)$:

```wolfram
(x^2+1)/((x+2)(x-1))==quo((x^2-1)/(x+3))+rem//Together
(* Output *)
True
```

$q$ and $r$ are uniquely determined by the condition that the degree of $r$ is less than the degree of $\frac{g}{d}$:

```wolfram
Exponent[Numerator[rem],x]<Exponent[Cancel[(x^2-1)/Denominator[rem]],x]
(* Output *)
True
```

### Options

#### Modulus

Use a prime modulus:

```wolfram
PolynomialQuotientRemainder[x^2+4 x+1,2x+1,x,Modulus->2]
(* Output *)
{1+x^2,0}
```

```wolfram
PolynomialQuotientRemainder[x^2+4 x+1,2x+1,x,Modulus->3]
(* Output *)
{1+2 x,0}
```

### Applications

Express the rational function as a polynomial and simple fraction:

```wolfram
{f,g}={x^2+2x+1,x+2};
```

```wolfram
{q,r}=PolynomialQuotientRemainder[f,g,x]
(* Output *)
{x,1}
```

The transformed rational function:

```wolfram
f/g==q+r/g
(* Output *)
(1+2 x+x^2)/(2+x)==x+(1)/(2+x)
```

```wolfram
Simplify[%]
(* Output *)
True
```

### Properties & Relations

For a polynomial $f$, $f=g q+r$:

```wolfram
{f,g}={x^2+4x+1,x+2};
```

```wolfram
{q,r}=PolynomialQuotientRemainder[f,g,x]
(* Output *)
{2+x,-3}
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to verify identity:

```wolfram
Expand[q g+r]==f
(* Output *)
True
```

[PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) and [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html):

```wolfram
{PolynomialQuotient[f,g,x],PolynomialRemainder[f,g,x]}
(* Output *)
{2+x,-3}
```

[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html) generalizes [PolynomialQuotientRemainder](https://reference.wolfram.com/language/ref/PolynomialQuotientRemainder.html) for multivariate polynomials:

```wolfram
PolynomialReduce[x^2+4x+1,{x+2},{x}]
(* Output *)
{{2+x},-3}
```

```wolfram
PolynomialQuotientRemainder[x^2+4x+1,x+2,x]
(* Output *)
{2+x,-3}
```

## Related Guides ▪Polynomial Division ▪Finite Fields

## History Introduced in 2007 (6.0) | Updated in 2023 (13.3)
