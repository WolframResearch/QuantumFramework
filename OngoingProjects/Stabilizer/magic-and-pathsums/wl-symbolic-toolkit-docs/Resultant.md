# Resultant | [SpanFromLeft]

> [Resultant](https://reference.wolfram.com/language/ref/Resultant.html)[*poly*_1,*poly*_2,*var*] — computes the resultant of the polynomials `*poly*_1 and `*poly*_2 with respect to the variable `*var*`.
> [Resultant](https://reference.wolfram.com/language/ref/Resultant.html)[*poly*_1,*poly*_2,*var*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*] — computes the resultant modulo the prime `*p*`.

## Details and Options

The resultant of two polynomials `*p*` and `*q*`, both with leading coefficient 1, is the product of all the differences `*p*_*i*-*q*_*j*` between roots of the polynomials. The resultant is always a number or a polynomial.

[Resultant](https://reference.wolfram.com/language/ref/Resultant.html) takes the following options:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |

## Examples

### Basic Examples

The resultant vanishes exactly when the polynomials have roots in common:

```wolfram
Resultant[(x-a)(x-b),(x-c)(x-d)(x-e),x]
(* Output *)
-(b-c) (-a+c) (b-d) (-a+d) (b-e) (-a+e)
```

### Scope

Resultant of polynomials with numeric coefficients:

```wolfram
Resultant[x^2-2x+7,x^3-x+5,x]
(* Output *)
265
```

Resultant of polynomials with parametric coefficients:

```wolfram
Resultant[a x^2+b x+c,(a-b)x^3+(a-c)x^2+(b-c)x+9,x]
(* Output *)
81 a^3-9 a b^3+9 b^4-18 a^3 c+36 a^2 b c-36 a b^2 c+a b^3 c-b^4 c+18 a^2 c^2+a^3 c^2-4 a^2 b c^2+3 a b^2 c^2+b^3 c^2+2 a^2 c^3-4 a b c^3+a c^4
```

Resultant over integers modulo 3:

```wolfram
Resultant[x^2-2x+7,x^3-x+5,x,Modulus->3]
(* Output *)
1
```

Resultant over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
Resultant[ℱ[1]x^5+ℱ[123]x y+ℱ[456],ℱ[789]x^3+ℱ[987]x+ℱ[654],x]
(* Output *)
![image](img/image_001.png)
```

The resultant reflects the multiplicities of roots:

```wolfram
Resultant[(x-a)^2(x-b),(x-c)^2(x-d)(x-e),x]
(* Output *)
(a-c)^4 (-b+c)^2 (b-d) (-a+d)^2 (b-e) (-a+e)^2
```

Compute the resultant of two polynomials of degree $1000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[1000]; q=rpoly[1000];
```

```wolfram
Resultant[p,q,x]//Short//AbsoluteTiming
(* Output *)
{2.9830545,149005787833207904470061883781050202587<<7863>>1745054253417528757002023825934056331228}
```

### Generalizations & Extensions

The resultant of rational functions is defined using the multiplicative property:

```wolfram
f_1=2 x+3;f_2=3x^2-7;g_1=x^2-x+5;g_2=2x^3-3x+8;
```

```wolfram
r1=Resultant[f_1f_2^2,g_1g_2^3,x]
(* Output *)
10245478724269744968742734862760
```

```wolfram
r2=Resultant[f_1,g_1,x]Resultant[f_2,g_1,x]^2Resultant[f_1,g_2,x]^3Resultant[f_2,g_2,x]^6
(* Output *)
10245478724269744968742734862760
```

```wolfram
r1==r2
(* Output *)
True
```

```wolfram
r1=Resultant[f_1/f_2,g_1/g_2^2,x]
(* Output *)
(84413315)/(979708)
```

```wolfram
r2=Resultant[f_1,g_1,x]Resultant[f_2,g_1,x]^-1Resultant[f_1,g_2,x]^-2Resultant[f_2,g_2,x]^2
(* Output *)
(84413315)/(979708)
```

```wolfram
r1==r2
(* Output *)
True
```

### Options

#### Method

This compares timings of the available methods of resultant computation:

```wolfram
Timing[Resultant[x^100-2x^77+3x+4,x^97+3x^55-11x^3+9,x,Method->#]//Short]&/@{Automatic,"Modular" ,"Subresultants", "BezoutMatrix","SylvesterMatrix"}//Column
(* Output *)
{{{0.,231215389619608480952934812126954<<52>>5924418056229593467896831334882576}}, {{0.,231215389619608480952934812126954<<52>>5924418056229593467896831334882576}}, {{0.,231215389619608480952934812126954<<52>>5924418056229593467896831334882576}}, {{0.0156001,231215389619608480952934812126954<<52>>5924418056229593467896831334882576}}, {{0.0156001,231215389619608480952934812126954<<52>>5924418056229593467896831334882576}}}
```

```wolfram
Timing[Resultant[a x^10+b x^5+(a+b)x+a b c,c x^5-a x^3+(a-b)x+7,x,Method->#]//Short]&/@{Automatic,"Modular" ,"Subresultants", "BezoutMatrix","SylvesterMatrix"}//Column
(* Output *)
![image](img/image_003.png)
```

#### Modulus

By default the resultant is computed over the rational numbers:

```wolfram
Resultant[a x^3+b x^2+c x+d,d x^3+c x^2+b x+a,x]
(* Output *)
a^6-a^4 b^2+2 a^3 b^2 c-2 a^4 c^2-a^2 b^2 c^2+a^2 c^4-2 a^2 b^3 d+6 a^3 b c d+2 a b^3 c d-4 a^2 b c^2 d-2 a b c^3 d-3 a^4 d^2-a^2 b^2 d^2-b^4 d^2+4 a b^2 c d^2+a^2 c^2 d^2+b^2 c^2 d^2+2 a c^3 d^2-6 a b c d^3-2 b c^2 d^3+3 a^2 d^4+2 b^2 d^4+c^2 d^4-d^6
```

Compute the resultant of the same polynomials over the integers modulo 2:

```wolfram
Resultant[a x^3+b x^2+c x+d,d x^3+c x^2+b x+a,x,Modulus->2]
(* Output *)
a^6+a^4 b^2+a^2 b^2 c^2+a^2 c^4+a^4 d^2+a^2 b^2 d^2+b^4 d^2+a^2 c^2 d^2+b^2 c^2 d^2+a^2 d^4+c^2 d^4+d^6
```

Compute the resultant of the same polynomials over the integers modulo 3:

```wolfram
Resultant[a x^3+b x^2+c x+d,d x^3+c x^2+b x+a,x,Modulus->3]
(* Output *)
a^6+2 a^4 b^2+2 a^3 b^2 c+a^4 c^2+2 a^2 b^2 c^2+a^2 c^4+a^2 b^3 d+2 a b^3 c d+2 a^2 b c^2 d+a b c^3 d+2 a^2 b^2 d^2+2 b^4 d^2+a b^2 c d^2+a^2 c^2 d^2+b^2 c^2 d^2+2 a c^3 d^2+b c^2 d^3+2 b^2 d^4+c^2 d^4+2 d^6
```

### Applications

Decide whether two polynomials have common roots:

```wolfram
Resultant[x^3-5x^2-7x+3,x^3-8x^2+9x-11,x]
(* Output *)
-10321
```

```wolfram
PolynomialGCD[x^3-5x^2-7x+3,x^3-8x^2+9x-11]
(* Output *)
1
```

```wolfram
Resultant[x^3-5x^2-7x+14,x^3-8x^2+9x+58,x]
(* Output *)
0
```

```wolfram
PolynomialGCD[x^3-5x^2-7x+14,x^3-8x^2+9x+58]
(* Output *)
2+x
```

Find conditions for two polynomials to have common roots:

```wolfram
Reduce[Resultant[x^3-2a x^2+a^2x-1,x^2-2a x+3,x]==0,a]
(* Output *)
a==2||a==Root||a==Root||a==Root
```

```wolfram
Apply[PolynomialGCD[##,Extension->Automatic]&,{x^3-2a x^2+a^2x-1,x^2-2a x+3}/.{ToRules[%]},{1}]
(* Output *)
{-1+x,-9+x-2 Root+3 Root^2,-9+x-2 Root+3 Root^2,-9+x-2 Root+3 Root^2}
```

### Properties & Relations

The resultant is zero if and only if the polynomials have a common root:

```wolfram
Resultant[(x-1)(x-2)(x-3),(x-4)(x-5)(x-6),x]
(* Output *)
-8640
```

```wolfram
Resultant[(x-1)(x-2)(x-3),(x-4)(x-5)(x-1),x]
(* Output *)
0
```

The polynomials have a zero resultant if and only if they have a nonconstant [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html):

```wolfram
Resultant[x^3-1,x^3+2x^2+2x-1,x]
(* Output *)
16
```

```wolfram
PolynomialGCD[x^3-1,x^3+2x^2+2x-1]
(* Output *)
1
```

```wolfram
Resultant[x^3-1,x^3+2x^2+2x+1,x]
(* Output *)
0
```

```wolfram
PolynomialGCD[x^3-1,x^3+2x^2+2x+1]
(* Output *)
1+x+x^2
```

The resultant can be represented in terms of roots as $R(f,g)=a_{n}^{m}b_{m}^{n}\prod_{1 \leq i \leq n,1 \leq j \leq m}(r_{i}-s_{j})$:

```wolfram
Resultant[a_3(x-r_1)(x-r_2)(x-r_3),b_2(x-s_1)(x-s_2),x]
(* Output *)
a_3^2 b_2^3 (r_1-s_1) (r_2-s_1) (r_3-s_1) (r_1-s_2) (r_2-s_2) (r_3-s_2)
```

Equation $(-1)^{n(n-1)/2}a_{n}D(f)=R(f,f')$ relates [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) and [Resultant](https://reference.wolfram.com/language/ref/Resultant.html):

```wolfram
f=3x^7-5x+4;
n=Exponent[f,x];
```

```wolfram
Resultant[f,D[f,x],x]
(* Output *)
4720053663936
```

```wolfram
(-1)^(n (n-1)/2)Coefficient[f,x,n]Discriminant[f,x]
(* Output *)
4720053663936
```

[GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) can also be used to find conditions for common roots:

```wolfram
Resultant[(x-a)(x-b),(x-c)(x-d)(x-e),x]
(* Output *)
-(b-c) (-a+c) (b-d) (-a+d) (b-e) (-a+e)
```

```wolfram
Select[GroebnerBasis[{(x-a)(x-b),(x-c)(x-d)(x-e)},x],FreeQ[#,x]&]//Factor
(* Output *)
{(a-c) (b-c) (a-d) (b-d) (a-e) (b-e)}
```

The same problem can also be solved using [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Resolve](https://reference.wolfram.com/language/ref/Resolve.html), and [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html):

```wolfram
Reduce[Exists[x,(x-a)(x-b)==0&&(x-c)(x-d)(x-e)==0],{a,b,c,d,e}]
(* Output *)
c==a||c==b||d==a||d==b||e==a||e==b
```

```wolfram
Resolve[Exists[x,(x-a)(x-b)==0&&(x-c)(x-d)(x-e)==0]]
(* Output *)
a==c||a==d||a==e||b==c||b==d||b==e
```

```wolfram
Eliminate[{(x-a)(x-b)==0,(x-c)(x-d)(x-e)==0},x]//Simplify
(* Output *)
(a-c) (-b+c) (a-d) (b-d) (a-e) (b-e)==0
```

### Possible Issues

The following two polynomials have no common root:

```wolfram
Resultant[x-1+10^(-17),x-1,x]
(* Output *)
-(1)/(100000000000000000)
```

Using approximate coefficients they will appear to have a common root:

```wolfram
Resultant[x-1.+10^(-17),x-1,x]
(* Output *)
0.
```

Using higher precision shows they have no common root:

```wolfram
Resultant[x-1+10^(-17),x-1,x]
(* Output *)
-1.×10^-17
```

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Systems ▪Polynomial Algebra ▪Finite Fields

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=Resultant) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2022 (13.2) ▪ 2023 (13.3)
