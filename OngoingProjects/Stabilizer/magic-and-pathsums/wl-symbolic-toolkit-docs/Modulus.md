# Modulus | [SpanFromLeft]

> [Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*n* — is an option that can be given in certain algebraic functions to specify that integers should be treated modulo `*n*`.

## Details

[Modulus](https://reference.wolfram.com/language/ref/Modulus.html) appears as an option in [Solve](https://reference.wolfram.com/language/ref/Solve.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Factor](https://reference.wolfram.com/language/ref/Factor.html), [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html), and [PolynomialLCM](https://reference.wolfram.com/language/ref/PolynomialLCM.html), as well as in linear algebra functions such as [Inverse](https://reference.wolfram.com/language/ref/Inverse.html), [LinearSolve](https://reference.wolfram.com/language/ref/LinearSolve.html), and [Det](https://reference.wolfram.com/language/ref/Det.html).

Arithmetic is usually done over the full ring $\mathbb{Z}$ of integers; setting the option [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) specifies that arithmetic should instead be done in the finite ring $\mathbb{Z}_{n}$.

The setting [Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->0 specifies the full ring $\mathbb{Z}$ of integers.

Some functions require that [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) be set to a prime, or a power of a prime. $\mathbb{Z}_{n}$ is a finite field when $n$ is prime.

Equations for [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) can be given in [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) and related functions.

## Examples

### Basic Examples

Solve equations:

```wolfram
Reduce[x^5==y^4+xy+1,{x,y},Modulus->4]
(* Output *)
(x==1&&y==0)||(x==1&&y==3)||(x==2&&y==1)||(x==2&&y==3)||(x==3&&y==2)||(x==3&&y==3)
```

Factor polynomials:

```wolfram
Table[Factor[x^6+3x^4+x,Modulus->Prime[k]],{k,1,3}]
(* Output *)
{x (1+x^3+x^5),x (1+x) (1+2 x+x^2+2 x^3+x^4),x (2+x) (4+x) (2+x+4 x^2+x^3)}
```

Compute inverse:

```wolfram
Inverse[{{1,2},{3,4}},Modulus->3]
(* Output *)
{{1,1},{0,1}}
```

### Scope

Compute [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html) over the integers modulo 2:

```wolfram
PolynomialGCD[x^3+x^2+x+1,(x+1)^4,Modulus->2]
(* Output *)
(1+x)^3
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) a polynomial over the integers modulo 3:

```wolfram
Factor[x^3+2x^2+3x+5,Modulus->3]
(* Output *)
(1+x) (2+x+x^2)
```

Find a [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) over the integers modulo 5:

```wolfram
GroebnerBasis[{x^3-2x y+y^2,y^3-x^2y+3},{x,y},Modulus->5]
(* Output *)
{2+4 y^2+2 y^3+2 y^4+y^5+4 y^6+3 y^7+y^8+y^9,4+x+y+4 y^3+4 y^6+y^7+4 y^8}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) equations over the integers modulo 7:

```wolfram
Reduce[x^3-3x y+2y^2==0&&y^3-x^2y+3==0,{x,y},Modulus->7]
(* Output *)
(x==1&&y==4)||(x==5&&y==1)
```

Compute the determinant of a matrix modulo 8:

```wolfram
Det[{{1,2,3},{21,210,2100},{2221,3332,4443}},Modulus->8]
(* Output *)
6
```

Find a modulus for which a system of equations has a solution:

```wolfram
Solve[x^2+y^2==1&&2x+3y==5&&x^3-x y==7,{x,y},Modulus->Automatic]
(* Output *)
{{Modulus->191721,y->96231,x->47377}}
```

### Properties & Relations

[Factor](https://reference.wolfram.com/language/ref/Factor.html) a polynomial over a finite field:

```wolfram
Factor[x^4-2,Modulus->3]
(* Output *)
(2+x+x^2) (2+2 x+x^2)
```

[Factor](https://reference.wolfram.com/language/ref/Factor.html) a polynomial over a finite [Extension](https://reference.wolfram.com/language/ref/Extension.html) of rationals:

```wolfram
Factor[x^4-2,Extension->Sqrt[2]]
(* Output *)
-(Sqrt[2]-x^2) (Sqrt[2]+x^2)
```

## Tech Notes ▪Polynomials Modulo Primes

## Related Guides ▪Polynomial Algebra ▪Polynomial Factoring & Decomposition ▪Polynomial Division

## History Introduced in 1988 (1.0)
