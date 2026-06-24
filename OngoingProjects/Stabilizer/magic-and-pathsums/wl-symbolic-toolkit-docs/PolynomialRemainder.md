# PolynomialRemainder | [SpanFromLeft]

> [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html)[*p*,*q*,*x*] — gives the remainder from dividing `*p*` by `*q*`, treated as polynomials in `*x*`.

## Details and Options

The degree of the result in `*x*` is guaranteed to be smaller than the degree of `*q*`.

Unlike [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html), [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) performs divisions in generating its results.

With the option [Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n*, the remainder is computed modulo `*n*`.

## Examples

### Basic Examples

Find the remainder after dividing one polynomial by another:

```wolfram
PolynomialRemainder[x^4+2x+1,x^2+1,x]
(* Output *)
2+2 x
```

The difference of the dividend and the remainder is a polynomial multiple of the divisor:

```wolfram
(x^4+2x+1-%)/(x^2+1)//Cancel
(* Output *)
-1+x^2
```

If the dividend is a multiple of the divisor, then the remainder is zero:

```wolfram
PolynomialRemainder[x^2+2x+1,x+1,x]
(* Output *)
0
```

Find the remainder of division for polynomials with symbolic coefficients:

```wolfram
PolynomialRemainder[x^3,a x+b,x]
(* Output *)
-(b^3)/(a^3)
```

Coefficients of the quotient are rational functions of the input coefficients:

```wolfram
(x^3-%)/(a x+b)//Cancel
(* Output *)
(b^2-a b x+a^2 x^2)/(a^3)
```

### Scope

The resulting polynomial will have coefficients that are rational numbers:

```wolfram
PolynomialRemainder[x^2+x+1,2x+1,x]
(* Output *)
(3)/(4)
```

The resulting polynomials will have coefficients that are rational functions:

```wolfram
{PolynomialRemainder[x^2+ b x+1,a x+1,x],PolynomialRemainder[x^2+x+1,Pi x+1,x]}
(* Output *)
{1+(1)/(a^2)-(b)/(a),1+(1)/(π^2)-(1)/(π)}
```

Polynomial remainder over the integers modulo $5$:

```wolfram
PolynomialRemainder[x^2+4 x+1,2x+1,x,Modulus->5]
(* Output *)
3
```

Polynomial remainder over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialRemainder[ℱ[1]x^5+ℱ[123]x+ℱ[456],ℱ[789]x^3+ℱ[987]x+ℱ[654],x]
(* Output *)
![image](img/image_001.png)
```

[PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) also works for rational functions:

```wolfram
PolynomialRemainder[(x^2+1)/((x+2)(x-1)),(x^2-1)/(x+3),x]
(* Output *)
(2)/(-1+x)
```

The quotient and remainder of division of $\frac{f}{h}$ by $\frac{g}{k}$ are $\frac{q k}{d h}$ and $\frac{r}{d}$, where $d=gcd(g,h)$:

```wolfram
((x^2+1)/((x+2)(x-1))-%)/((x^2-1)/(x+3))//Together
(* Output *)
((-3+x) (3+x))/((-1+x)^2 (2+x))
```

$q$ and $r$ are uniquely determined by the condition that the degree of $r$ is less than the degree of $\frac{g}{d}$:

```wolfram
Exponent[Numerator[%%],x]<Exponent[Cancel[(x^2-1)/Denominator[%%]],x]
(* Output *)
True
```

### Options

#### Modulus

Use 2 as a modulus:

```wolfram
PolynomialRemainder[x^2+4 x+1,2x+1,x,Modulus->2]
(* Output *)
0
```

This time use 5:

```wolfram
PolynomialRemainder[x^2+4 x+1,2x+1,x,Modulus->5]
(* Output *)
3
```

### Applications

Euclid's algorithm for the greatest common divisor:

```wolfram
Euclid[f_,g_,x_]/;Exponent[f,x]<Exponent[g,x] :=Euclid[g,f,x]
Euclid[f_,g_?NumericQ,x_]:=If[PossibleZeroQ[g],f,1]
Euclid[f_,g_,x_]:=Euclid[g,PolynomialRemainder[f,g,x],x]
```

Use it to compute a GCD:

```wolfram
Euclid[(x+1)(x+2),(x+2)(x+3),x]
(* Output *)
-4-2 x
```

Divide by the leading coefficient:

```wolfram
Expand[%/Coefficient[%,x,1]]
(* Output *)
2+x
```

```wolfram
PolynomialGCD[(x+1)(x+2),(x+2)(x+3)]
(* Output *)
2+x
```

### Properties & Relations

For a polynomial $f$, $f=\mathit{qg}+r$, where $q$ is given by [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html):

```wolfram
{f,g}={x^2+4x+1,x+2};
q=PolynomialQuotient[f,g,x]
(* Output *)
2+x
```

And `*r*` is given by [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html):

```wolfram
r=PolynomialRemainder[f,g,x]
(* Output *)
-3
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to verify the identity:

```wolfram
Expand[q g+r]==f
(* Output *)
True
```

To get both quotient and remainder use [PolynomialQuotientRemainder](https://reference.wolfram.com/language/ref/PolynomialQuotientRemainder.html):

```wolfram
PolynomialQuotientRemainder[f,g,x]
(* Output *)
{2+x,-3}
```

[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html) generalizes [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) for multivariate polynomials:

```wolfram
PolynomialReduce[x^2+4x+1,{x+2},{x}]
(* Output *)
{{2+x},-3}
```

[PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) requires knowing the variable of interest:

```wolfram
PolynomialRemainder[x^2,x+a,x]
(* Output *)
a^2
```

[PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) is not provided that information and can choose when more than one variable is present:

```wolfram
PolynomialMod[x^2,x+a]
(* Output *)
x^2
```

### Possible Issues

The variable assumed for the polynomials matters:

```wolfram
{PolynomialRemainder[x+y,x-y,x],PolynomialRemainder[x+y,x-y,y]}
(* Output *)
{2 y,2 x}
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Division ▪Rational Functions ▪Finite Fields

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=PolynomialRemainder) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2023 (13.3)
