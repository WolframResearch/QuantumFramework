# PolynomialQuotient | [SpanFromLeft]

> [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html)[*p*,*q*,*x*] — gives the quotient of `*p*` and `*q*`, treated as polynomials in `*x*`, with any remainder dropped.

## Details and Options

With the option [Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n*, the quotient is computed modulo `*n*`.

## Examples

### Basic Examples

The quotient of two polynomials:

```wolfram
PolynomialQuotient[x^4+2x+1,x^2+1,x]
(* Output *)
-1+x^2
```

The degree of the remainder is less than the degree of the divisor:

```wolfram
x^4+2x+1-%(x^2+1)//Expand
(* Output *)
2+2 x
```

The quotient of $x^{2}$ by $x+a$, with the remainder dropped:

```wolfram
PolynomialQuotient[x^2,x+a,x]
(* Output *)
-a+x
```

```wolfram
x^2-% (x+a)//Expand
(* Output *)
a^2
```

If the degree of the dividend is less than the degree of the divisor, then the quotient is zero:

```wolfram
PolynomialQuotient[x^2+2x+1,x^3,x]
(* Output *)
0
```

### Scope

The resulting polynomial will have coefficients that are rational expressions of input coefficients:

```wolfram
PolynomialQuotient[x^2+x+1,2x+1,x]
(* Output *)
(1)/(4)+(x)/(2)
```

```wolfram
PolynomialQuotient[x^2+ b x+1,a x+1,x]
(* Output *)
-(1)/(a^2)+(b)/(a)+(x)/(a)
```

```wolfram
PolynomialQuotient[x^2+x+1,Pi x+1,x]
(* Output *)
-(1)/(π^2)+(1)/(π)+(x)/(π)
```

Polynomial quotient over the integers modulo $5$:

```wolfram
PolynomialQuotient[x^2+4 x+1,2x+1,x,Modulus->5]
(* Output *)
3+3 x
```

Polynomial quotient over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialQuotient[ℱ[1]x^5+ℱ[123]x+ℱ[456],ℱ[789]x^3+ℱ[987]x+ℱ[654],x]
(* Output *)
![image](img/image_001.png)
```

[PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) also works for rational functions:

```wolfram
PolynomialQuotient[(x^2+1)/((x+2)(x-1)),(x^2-1)/(x+3),x]
(* Output *)
(-9+x^2)/((-1+x) (-2+x+x^2))
```

The quotient and remainder of division of $\frac{f}{h}$ by $\frac{g}{k}$ are $\frac{q k}{d h}$ and  $\frac{r}{d}$, where $d=gcd(g,h)$:

```wolfram
(x^2+1)/((x+2)(x-1))-%((x^2-1)/(x+3))//Together
(* Output *)
(2)/(-1+x)
```

$q$ and $r$ are uniquely determined by the condition that the degree of $r$ is less than the degree of $\frac{g}{d}$:

```wolfram
Exponent[Numerator[%],x]<Exponent[Cancel[(x^2-1)/Denominator[%]],x]
(* Output *)
True
```

### Options

#### Modulus

Use a prime modulus:

```wolfram
PolynomialQuotient[x^2+4 x+1,2x+1,x,Modulus->2]
(* Output *)
1+x^2
```

```wolfram
PolynomialQuotient[x^2+4 x+1,2x+1,x,Modulus->3]
(* Output *)
1+2 x
```

### Applications

When the divisor $g$ divides the dividend $f$, then the quotient $q$ of $f$ by $g$ satisfies $f=q g$:

```wolfram
f=x^4-4;
g=x^2-2;
```

Use [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) to check that $g$ divides $f$:

```wolfram
PolynomialGCD[f,g]==g
(* Output *)
True
```

Verify that $f=q g$:

```wolfram
q=PolynomialQuotient[f,g,x]
(* Output *)
2+x^2
```

```wolfram
Expand[f-q g]
(* Output *)
0
```

In general, the quotient $q$ of $f$ by $g$ satisfies $f=q g+r$:

```wolfram
f=x^4-5x+4;
g=x^2-2;
```

```wolfram
q=PolynomialQuotient[f,g,x]
(* Output *)
2+x^2
```

The degree of the remainder $r$ is less than the degree of $g$:

```wolfram
r=Expand[f-q g]
(* Output *)
8-5 x
```

```wolfram
Exponent[r,x]<Exponent[g,x]
(* Output *)
True
```

Factor a polynomial by finding one root at a time:

```wolfram
f=x^3+6x^2+11x+6;
```

```wolfram
FindRoot[f,{x,0}]
(* Output *)
{x->-1.0000000000000004}
```

Take a quotient by the first factor:

```wolfram
PolynomialQuotient[f,x+1,x]
(* Output *)
6+5 x+x^2
```

Find another root and compute the quotient:

```wolfram
FindRoot[%,{x,0}]
(* Output *)
{x->-2.0000000000000004}
```

```wolfram
PolynomialQuotient[%%,x+2,x]
(* Output *)
3+x
```

Verify the obtained factorization:

```wolfram
f-(x+1)(x+2)(x+3)//Expand
(* Output *)
0
```

### Properties & Relations

For a polynomial `*f*`, `*f*==*gq*+*r*`, where `*r*` is given by [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html):

```wolfram
{f,g}={x^2+4x+1,x+2};
```

```wolfram
q=PolynomialQuotient[f,g,x]
(* Output *)
2+x
```

```wolfram
r=PolynomialRemainder[f,g,x]
(* Output *)
-3
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to verify identity:

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

[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html) generalizes [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) for multivariate polynomials:

```wolfram
PolynomialReduce[x^2+4x+1,{x+2},{x}]
(* Output *)
{{2+x},-3}
```

Use [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) to find a common divisor:

```wolfram
{f,g}={x^2+3x+2,x^2+5x+6};
```

```wolfram
h=PolynomialGCD[f,g]
(* Output *)
2+x
```

Use [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) to see the resulting factorization:

```wolfram
f==PolynomialQuotient[f,h,x]h
(* Output *)
2+3 x+x^2==(1+x) (2+x)
```

```wolfram
g==PolynomialQuotient[g,h,x]h
(* Output *)
6+5 x+x^2==(2+x) (3+x)
```

For rational functions common divisors are not automatically canceled:

```wolfram
f/g
(* Output *)
(2+3 x+x^2)/(6+5 x+x^2)
```

[Cancel](https://reference.wolfram.com/language/ref/Cancel.html) effectively uses [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) to cancel common divisors:

```wolfram
Cancel[%]
(* Output *)
(1+x)/(3+x)
```

The [Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html) polynomials are defined as quotients:

```wolfram
PolynomialQuotient[x^5-1,x-1,x]
(* Output *)
1+x+x^2+x^3+x^4
```

```wolfram
Cyclotomic[5,x]
(* Output *)
1+x+x^2+x^3+x^4
```

### Possible Issues

The result depends on what is assumed to be a variable:

```wolfram
{PolynomialQuotient[x^3+y^2,x-y,x],PolynomialQuotient[x^3+y^2,x-y,y]}
(* Output *)
{x^2+x y+y^2,-x-y}
```

The result from [PolynomialQuotient](https://reference.wolfram.com/language/ref/PolynomialQuotient.html) depends on recognizing zeros:

```wolfram
PolynomialQuotient[x^3+x+1,zero x^2+x+1,x]
(* Output *)
-(1)/(zero^2)+(x)/(zero)
```

```wolfram
PolynomialQuotient[x^3+x+1,x+1,x]
(* Output *)
2-x+x^2
```

This is a hidden zero:

```wolfram
zero=Sin[Sqrt[2]+Sqrt[3]]-Sin[Sqrt[5+2Sqrt[6]]];
```

```wolfram
FullSimplify[zero]
(* Output *)
0
```

The result is as if the hidden zero was not zero:

```wolfram
PolynomialQuotient[x^3+x+1,zero x^2+x+1,x]//N
(* Output *)
-5.070602400912918×10^30-2.251799813685248×10^15 x
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Division ▪Polynomial Algebra ▪Rational Functions ▪Finite Fields

## History Introduced in 1988 (1.0) | Updated in 2023 (13.3)
