# PolynomialExtendedGCD | [SpanFromLeft]

> [PolynomialExtendedGCD](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html)[*poly*_1,*poly*_2,*x*]  — gives the extended GCD of `*poly*_1 and `*poly*_2 treated as univariate polynomials in `*x*`.
> [PolynomialExtendedGCD](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html)[*poly*_1,*poly*_2,*x*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*]  — gives the extended GCD over the integers modulo the prime `*p*`.

## Details and Options

[PolynomialExtendedGCD](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html) takes the following options:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |

## Examples

### Basic Examples

Compute the extended GCD:

```wolfram
{f,g}={2x^5-2x,(x^2-1)^2};
```

```wolfram
{d,{a,b}}=PolynomialExtendedGCD[f,g,x]
(* Output *)
{-1+x^2,{(x)/(4),(1)/(4) (-4-2 x^2)}}
```

The second part gives coefficients of a linear combination of polynomials that yields the GCD:

```wolfram
a f+b g==d//Expand
(* Output *)
True
```

Compute the extended GCD of polynomials with coefficients involving symbolic parameters:

```wolfram
PolynomialExtendedGCD[a(x+b)^2,(x+a)(x+b),x]
(* Output *)
{b+x,{-(1)/(a (a-b)),(1)/(a-b)}}
```

### Scope

Polynomials with numeric coefficients:

```wolfram
PolynomialExtendedGCD[(x-1)(x-2)^2,(x-1)(x^2-3),x]
(* Output *)
{-1+x,{7+4 x,9-4 x}}
```

Polynomials with symbolic coefficients:

```wolfram
PolynomialExtendedGCD[(x-a)(b x-c)^2,(x-a)(x^2-b c),x]
(* Output *)
{-a+x,{(b^3 c+c^2+2 b c x)/(b^6 c^2-2 b^3 c^3+c^4),(-b^5 c+3 b^2 c^2-2 b^3 c x)/(b^6 c^2-2 b^3 c^3+c^4)}}
```

Relatively prime polynomials:

```wolfram
PolynomialExtendedGCD[a x^2+b x+c,x-r,x]
(* Output *)
{1,{(1)/(c+b r+a r^2),(-b-a r-a x)/(c+b r+a r^2)}}
```

Polynomials with complex coefficients:

```wolfram
PolynomialExtendedGCD[(x-I)(x^2+1),(x-1)(x+I),x]
(* Output *)
{ⅈ+x,{(ⅈ)/(2),(1)/(2) ⅈ ((-1+2 ⅈ)-x)}}
```

Compute the extended GCD of polynomials over the integers modulo 3:

```wolfram
PolynomialExtendedGCD[(x-1)^2(x-2)^3,x^3-1,x,Modulus->3]
(* Output *)
{1+x+x^2,{2,1+x+x^2}}
```

Compute the extended GCD of polynomials over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
PolynomialExtendedGCD[ℱ[1]x^2+ℱ[246]x+ℱ[4436],ℱ[3]x^3+ℱ[1771]x,x]
(* Output *)
![image](img/image_001.png)
```

### Options

#### Modulus

Extended GCD over the integers:

```wolfram
PolynomialExtendedGCD[(x-1)^2(x-2)^2,(x-1)(x^2-3),x]
(* Output *)
{-1+x,{(1)/(2) (19+11 x),(1)/(2) (-26+36 x-11 x^2)}}
```

Extended GCD over the integers modulo 2:

```wolfram
PolynomialExtendedGCD[(x-1)^2(x-2)^2,(x-1)(x^2-3),x,Modulus->2]
(* Output *)
{1+x^2,{1,1+x}}
```

### Applications

Given polynomials $a$, $b$, and $c$, find polynomials $f$ and $g$ such that $a f+b g=c$:

```wolfram
a=(x-1)(x^2+7x-9);
b=(x-1)^2(x^3-5x+3);
c=(x-1)(x^2-7);
```

```wolfram
{d,{r,s}}=PolynomialExtendedGCD[a,b,x]
(* Output *)
{-1+x,{(1)/(579) (1244-2360 x+53 x^2+484 x^3),(1)/(579) (-3925-484 x)}}
```

A solution exists if and only if $c$ is divisible by $d$:

```wolfram
h=Cancel[c/d]
(* Output *)
-7+x^2
```

```wolfram
{f,g}=h{r,s}
(* Output *)
{(1)/(579) (-7+x^2) (1244-2360 x+53 x^2+484 x^3),(1)/(579) (-3925-484 x) (-7+x^2)}
```

```wolfram
a f+b g==c//Expand
(* Output *)
True
```

### Properties & Relations

The extended GCD of $f$ and $g$ is `{*d*,{*r*,*s*}}`, such that $d=GCD(f,g)$ and $f r+g s=d$:

```wolfram
f=2a(x^2-1)^2(x^3-2);
g=4a(x-1)^3(x^2-3);
```

```wolfram
{d,{r,s}}=PolynomialExtendedGCD[f,g,x]
(* Output *)
{1-2 x+x^2,{-(166-32 x-42 x^2)/(736 a),-(-12+188 x+224 x^2+158 x^3+42 x^4)/(1472 a)}}
```

```wolfram
r f+s g==d//Expand
(* Output *)
True
```

`*d*` is equal to [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html)[*f*,*g*] up to a factor not containing `*x*`:

```wolfram
PolynomialGCD[f,g]
(* Output *)
2 a (-1+x)^2
```

```wolfram
Cancel[d/%]
(* Output *)
(1)/(2 a)
```

$\mathit{r}$ and $\mathit{s}$ are uniquely determined by the following [Exponent](https://reference.wolfram.com/language/ref/Exponent.html) conditions:

```wolfram
Exponent[r,x]<Exponent[g,x]-Exponent[d,x]
(* Output *)
True
```

```wolfram
Exponent[s,x]<Exponent[f,x]-Exponent[d,x]
(* Output *)
True
```

Use [Cancel](https://reference.wolfram.com/language/ref/Cancel.html) or [PolynomialRemainder](https://reference.wolfram.com/language/ref/PolynomialRemainder.html) to prove that `*d*` divides `*f*` and `*g*`:

```wolfram
Cancel[f/d]
(* Output *)
2 a (1+x)^2 (-2+x^3)
```

```wolfram
PolynomialRemainder[g,d,x]
(* Output *)
0
```

## Related Guides ▪Polynomial Division ▪Finite Fields

## History Introduced in 2007 (6.0) | Updated in 2023 (13.3)
