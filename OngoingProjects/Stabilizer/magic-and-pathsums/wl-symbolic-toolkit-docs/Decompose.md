# Decompose | [SpanFromLeft]

> [Decompose](https://reference.wolfram.com/language/ref/Decompose.html)[*poly*,*x*] — decomposes a polynomial, if possible, into a composition of simpler polynomials.

## Details and Options

[Decompose](https://reference.wolfram.com/language/ref/Decompose.html) gives a list of the polynomials `*P*_*i*` which can be composed as $P_{1} (P_{2} (\ldots x \ldots))$ to give the original polynomial.

The set of polynomials `*P*_*i*` is not necessarily unique.

Decomposition is an operation which is independent of polynomial factorization.

## Examples

### Basic Examples

Represent a polynomial as a composition of polynomials:

```wolfram
Decompose[x^2+1,x]
(* Output *)
{1+x,x^2}
```

```wolfram
(1+x)/.(x->x^2)
(* Output *)
1+x^2
```

Compositions of the same polynomials in different orders:

```wolfram
Decompose[x^4+x^2,x]
(* Output *)
{x+x^2,x^2}
```

```wolfram
Decompose[(x+x^2)^2,x]
(* Output *)
{x^2,x+x^2}
```

A decomposition with three polynomials:

```wolfram
Decompose[(x^2+x)^4+1, x]
(* Output *)
{1+x,x^4,x+x^2}
```

### Scope

A composition of more than two polynomials:

```wolfram
Decompose[(x^2+1)^4+3, x]
(* Output *)
{4+2 x+x^2,2 x+x^2,x^2}
```

No decomposition:

```wolfram
Decompose[2 x+1,x]
(* Output *)
{1+2 x}
```

A polynomial with symbolic coefficients:

```wolfram
Decompose[(a x^3+1)^2+b,x]
(* Output *)
{1+b+2 a x+a^2 x^2,x^3}
```

A polynomial with complex coefficients:

```wolfram
Decompose[(x^2+I)^2,x]
(* Output *)
{-1+2 ⅈ x+x^2,x^2}
```

Decompose a polynomial over the integers modulo 3:

```wolfram
Decompose[x^4+4x^3+4x^2+6x+24,x,Modulus->3]
(* Output *)
{x^2,2 x+x^2}
```

### Options

#### Modulus

Decompose a polynomial over integers modulo 2:

```wolfram
Decompose[x^3+2x^2+4x+7,x,Modulus->2]
(* Output *)
{1+x,x^3}
```

### Applications

Solve some polynomial equations of degrees higher than 4 in terms of radicals:

```wolfram
f=x^8+4x^7+2x^6-8x^5-5x^4+8x^3+2x^2-4x+8;
```

```wolfram
{a,b,c}=Decompose[f,x]
(* Output *)
{8+2 x+x^2,-2 x+x^2,x+x^2}
```

Solve $a(b(c(x)))=0$ by solving $a(y)=0$ and then $b(z)=y_{i}$ etc:

```wolfram
s1=x/.Solve[a==0,x]
(* Output *)
{-1-ⅈ Sqrt[7],-1+ⅈ Sqrt[7]}
```

```wolfram
s2=Join@@(x/.Solve[b==#,x]&/@s1)
(* Output *)
{1-(-1)^(3/4) 7^(1/4),1+(-1)^(3/4) 7^(1/4),1-(-7)^(1/4),1+(-7)^(1/4)}
```

```wolfram
s3=Join@@(x/.Solve[c==#,x]&/@s2)
(* Output *)
{(1)/(2) (-1-Sqrt[5-4 (-1)^(3/4) 7^(1/4)]),(1)/(2) (-1+Sqrt[5-4 (-1)^(3/4) 7^(1/4)]),(1)/(2) (-1-Sqrt[5+4 (-1)^(3/4) 7^(1/4)]),(1)/(2) (-1+Sqrt[5+4 (-1)^(3/4) 7^(1/4)]),(1)/(2) (-1-Sqrt[5-4 (-7)^(1/4)]),(1)/(2) (-1+Sqrt[5-4 (-7)^(1/4)]),(1)/(2) (-1-Sqrt[5+4 (-7)^(1/4)]),(1)/(2) (-1+Sqrt[5+4 (-7)^(1/4)])}
```

Check that these indeed are the roots of `f`:

```wolfram
Expand[f/.x->#&/@s3]
(* Output *)
{0,0,0,0,0,0,0,0}
```

Wolfram Language solvers use [Decompose](https://reference.wolfram.com/language/ref/Decompose.html) automatically:

```wolfram
Solve[f==0,x]
(* Output *)
{{x->(1)/(2) (-1-Sqrt[5-4 (-7)^(1/4)])},{x->(1)/(2) (-1+Sqrt[5-4 (-7)^(1/4)])},{x->(1)/(2) (-1-Sqrt[5+4 (-7)^(1/4)])},{x->(1)/(2) (-1+Sqrt[5+4 (-7)^(1/4)])},{x->(1)/(2) (-1-Sqrt[5-4 (-1)^(3/4) 7^(1/4)])},{x->(1)/(2) (-1+Sqrt[5-4 (-1)^(3/4) 7^(1/4)])},{x->(1)/(2) (-1-Sqrt[5+4 (-1)^(3/4) 7^(1/4)])},{x->(1)/(2) (-1+Sqrt[5+4 (-1)^(3/4) 7^(1/4)])}}
```

### Properties & Relations

Composition of polynomials given by [Decompose](https://reference.wolfram.com/language/ref/Decompose.html) gives the original polynomial:

```wolfram
f=(x^2+1)^4+3;
```

```wolfram
Decompose[f, x]
(* Output *)
{4+2 x+x^2,2 x+x^2,x^2}
```

Use [Fold](https://reference.wolfram.com/language/ref/Fold.html) to compose the polynomials:

```wolfram
Fold[#2/.x->#1&,x,Reverse[%]]
(* Output *)
4+2 (2 x^2+x^4)+(2 x^2+x^4)^2
```

Use [Expand](https://reference.wolfram.com/language/ref/Expand.html) to show that the result is equal to `f`:

```wolfram
Expand[%==f]
(* Output *)
True
```

Use [Factor](https://reference.wolfram.com/language/ref/Factor.html) to represent a polynomial as a product of irreducible factors:

```wolfram
f=x^2+3x+2;
g=x^4-x^2+17;
```

```wolfram
Factor[{f,g}]
(* Output *)
{(1+x) (2+x),17-x^2+x^4}
```

`f` can be factored but not decomposed; `g` can be decomposed but not factored:

```wolfram
Decompose[f,x]
(* Output *)
{2+3 x+x^2}
```

```wolfram
Decompose[g, x]
(* Output *)
{17-x+x^2,x^2}
```

### Possible Issues

[Decompose](https://reference.wolfram.com/language/ref/Decompose.html) ignores possible decompositions with inner polynomials that are linear:

```wolfram
Decompose[(x+1)^3+1,x]
(* Output *)
{2+3 x+3 x^2+x^3}
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Factoring & Decomposition ▪Polynomial Algebra

## History Introduced in 1988 (1.0)
