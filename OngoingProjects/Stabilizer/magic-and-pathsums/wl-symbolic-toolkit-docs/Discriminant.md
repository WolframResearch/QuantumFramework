# Discriminant | [SpanFromLeft]

> [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html)[*poly*,*var*]  — computes the discriminant of the polynomial `*poly*` with respect to the variable `*var*`.
> [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html)[*poly*,*var*,[Modulus](https://reference.wolfram.com/language/ref/Modulus.html)->*p*]  — computes the discriminant modulo $p$.

## Details and Options

The discriminant of a polynomial with leading coefficient one is the product over all pairs of roots $x_{i}$, $x_{j}$ of $(x_{i}-x_{j})^{2}$.

[Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) takes the following options:

| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| --- | --- | --- |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | modulus to assume for integers |

Typical possible values for [Method](https://reference.wolfram.com/language/ref/Method.html) are [Automatic](https://reference.wolfram.com/language/ref/Automatic.html), `"SylvesterMatrix"`, `"BezoutMatrix"`, `"Subresultants"` and `"Modular"`.

## Examples

### Basic Examples

Discriminant of a quadratic:

```wolfram
Discriminant[a x^2+b x+c,x]
(* Output *)
b^2-4 a c
```

### Scope

Discriminant of a polynomial with numeric coefficients:

```wolfram
Discriminant[x^10-5x^7-3 x+9,x]
(* Output *)
177945374758153510836
```

Discriminant of a general cubic:

```wolfram
Discriminant[a x^3+b x^2+c x +d,x]
(* Output *)
b^2 c^2-4 a c^3-4 b^3 d+18 a b c d-27 a^2 d^2
```

Discriminant of a general quintic:

```wolfram
Discriminant[Sum[a_ix^i,{i,0,5}],x]
(* Output *)
a_1^2 a_2^2 a_3^2 a_4^2-4 a_0 a_2^3 a_3^2 a_4^2-4 a_1^3 a_3^3 a_4^2+18 a_0 a_1 a_2 a_3^3 a_4^2-27 a_0^2 a_3^4 a_4^2-4 a_1^2 a_2^3 a_4^3+16 a_0 a_2^4 a_4^3+18 a_1^3 a_2 a_3 a_4^3-80 a_0 a_1 a_2^2 a_3 a_4^3-6 a_0 a_1^2 a_3^2 a_4^3+144 a_0^2 a_2 a_3^2 a_4^3-27 a_1^4 a_4^4+144 a_0 a_1^2 a_2 a_4^4-128 a_0^2 a_2^2 a_4^4-192 a_0^2 a_1 a_3 a_4^4+256 a_0^3 a_4^5-4 a_1^2 a_2^2 a_3^3 a_5+16 a_0 a_2^3 a_3^3 a_5+16 a_1^3 a_3^4 a_5-72 a_0 a_1 a_2 a_3^4 a_5+108 a_0^2 a_3^5 a_5+18 a_1^2 a_2^3 a_3 a_4 a_5-72 a_0 a_2^4 a_3 a_4 a_5-80 a_1^3 a_2 a_3^2 a_4 a_5+356 a_0 a_1 a_2^2 a_3^2 a_4 a_5+24 a_0 a_1^2 a_3^3 a_4 a_5-630 a_0^2 a_2 a_3^3 a_4 a_5-6 a_1^3 a_2^2 a_4^2 a_5+24 a_0 a_1 a_2^3 a_4^2 a_5+144 a_1^4 a_3 a_4^2 a_5-746 a_0 a_1^2 a_2 a_3 a_4^2 a_5+560 a_0^2 a_2^2 a_3 a_4^2 a_5+1020 a_0^2 a_1 a_3^2 a_4^2 a_5-36 a_0 a_1^3 a_4^3 a_5+160 a_0^2 a_1 a_2 a_4^3 a_5-1600 a_0^3 a_3 a_4^3 a_5-27 a_1^2 a_2^4 a_5^2+108 a_0 a_2^5 a_5^2+144 a_1^3 a_2^2 a_3 a_5^2-630 a_0 a_1 a_2^3 a_3 a_5^2-128 a_1^4 a_3^2 a_5^2+560 a_0 a_1^2 a_2 a_3^2 a_5^2+825 a_0^2 a_2^2 a_3^2 a_5^2-900 a_0^2 a_1 a_3^3 a_5^2-192 a_1^4 a_2 a_4 a_5^2+1020 a_0 a_1^2 a_2^2 a_4 a_5^2-900 a_0^2 a_2^3 a_4 a_5^2+160 a_0 a_1^3 a_3 a_4 a_5^2-2050 a_0^2 a_1 a_2 a_3 a_4 a_5^2+2250 a_0^3 a_3^2 a_4 a_5^2-50 a_0^2 a_1^2 a_4^2 a_5^2+2000 a_0^3 a_2 a_4^2 a_5^2+256 a_1^5 a_5^3-1600 a_0 a_1^3 a_2 a_5^3+2250 a_0^2 a_1 a_2^2 a_5^3+2000 a_0^2 a_1^2 a_3 a_5^3-3750 a_0^3 a_2 a_3 a_5^3-2500 a_0^3 a_1 a_4 a_5^3+3125 a_0^4 a_5^4
```

Discriminants are squares of differences of roots:

```wolfram
Discriminant[(x-a)(x-b)(x-c)(x-d),x]
(* Output *)
(a-b)^2 (a-c)^2 (b-c)^2 (a-d)^2 (b-d)^2 (c-d)^2
```

Discriminant over integers modulo 3:

```wolfram
Discriminant[x^5-x y+y^2-1,x,Modulus->3]
(* Output *)
2+y^2+2 y^5+y^6+2 y^8
```

Discriminant over a finite field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
Discriminant[ℱ[1]x^5+ℱ[123]x y+ℱ[456],x]
(* Output *)
![image](img/image_001.png)
```

Compute the discriminant of a polynomial of degree $1000$:

```wolfram
rpoly[n_]:=RandomInteger[{-2^10,2^10},{n+1}].x^Range[0,n]
SeedRandom[1234];
p=rpoly[1000];
```

```wolfram
Discriminant[p,x]//Short//AbsoluteTiming
(* Output *)
{2.2045553,-77921107669926779964331155238639496107952<<10256>>692072454910452578130689479008673384862896}
```

### Options

#### Method

This compares timings of the available methods of discriminant computation:

```wolfram
Timing[Discriminant[x^100-2x^77+3x+4,x,Method->#]//Short]&/@{Automatic,"Modular" ,"Subresultants", "BezoutMatrix","SylvesterMatrix"}//Column
(* Output *)
{{{0.,542284250071622148740342384159<<200>>551931558053046829820246441557}}, {{0.,542284250071622148740342384159<<200>>551931558053046829820246441557}}, {{0.,542284250071622148740342384159<<200>>551931558053046829820246441557}}, {{0.,542284250071622148740342384159<<200>>551931558053046829820246441557}}, {{0.0312002,542284250071622148740342384159<<200>>551931558053046829820246441557}}}
```

```wolfram
Timing[Discriminant[a x^10+b x^5+(a+b)x+a b c,x,Method->#]//Short]&/@{Automatic,"Modular" ,"Subresultants", "BezoutMatrix","SylvesterMatrix"}//Column
(* Output *)
{{{0.0156001,387420489 a^18+3874204890 a^17 b+<<55>>+<<1>>-10000000000 a^18 b^9 c^9}}, {{0.0312002,387420489 a^18+3874204890 a^17 b+<<55>>+<<1>>-10000000000 a^18 b^9 c^9}}, {{0.0156001,387420489 a^18+3874204890 a^17 b+<<55>>+<<1>>-10000000000 a^18 b^9 c^9}}, {{0.0312002,387420489 a^18+3874204890 a^17 b+<<55>>+<<1>>-10000000000 a^18 b^9 c^9}}, {{0.0156001,387420489 a^18+3874204890 a^17 b+<<55>>+<<1>>-10000000000 a^18 b^9 c^9}}}
```

#### Modulus

By default the discriminant is computed over the rational numbers:

```wolfram
Discriminant[a x^3+b x^2+c x+d,x]
(* Output *)
b^2 c^2-4 a c^3-4 b^3 d+18 a b c d-27 a^2 d^2
```

Compute the discriminant of the same polynomial over the integers modulo 2:

```wolfram
Discriminant[a x^3+b x^2+c x+d,x,Modulus->2]
(* Output *)
b^2 c^2+a^2 d^2
```

Compute the discriminant of the same polynomial over the integers modulo 3:

```wolfram
Discriminant[a x^3+b x^2+c x+d,x,Modulus->3]
(* Output *)
b^2 c^2+2 a c^3+2 b^3 d
```

### Applications

Decide whether a polynomial has multiple roots:

```wolfram
Discriminant[x^11-2x^10+x^9-2x^3+11x^2-16x+7,x]
(* Output *)
0
```

```wolfram
FactorSquareFree[x^11-2x^10+x^9-2x^3+11x^2-16x+7]
(* Output *)
(-1+x)^2 (7-2 x+x^9)
```

```wolfram
Discriminant[x^11-2x^10+x^9-2x^3+11x^2-16x-7,x]
(* Output *)
813583724618463631662004
```

```wolfram
FactorSquareFree[x^11-2x^10+x^9-2x^3+11x^2-16x-7]
(* Output *)
-7-16 x+11 x^2-2 x^3+x^9-2 x^10+x^11
```

Find the condition for a cubic to have multiple roots:

```wolfram
Solve[Discriminant[x^3+x+ c ,x]==0,c]
(* Output *)
{{c->-(2 ⅈ)/(3 Sqrt[3])},{c->(2 ⅈ)/(3 Sqrt[3])}}
```

```wolfram
FactorSquareFree[x^3+x+ c/.%,Extension->Automatic]
(* Output *)
{-(1)/(27) ⅈ (2 Sqrt[3]-3 ⅈ x) (Sqrt[3]+3 ⅈ x)^2,(1)/(27) ⅈ (Sqrt[3]-3 ⅈ x)^2 (2 Sqrt[3]+3 ⅈ x)}
```

### Properties & Relations

The discriminant is zero if and only if the polynomial has multiple roots:

```wolfram
Discriminant[(x-1)(x-2)(x-3),x]
(* Output *)
4
```

```wolfram
Discriminant[(x-1)(x-2)(x-1),x]
(* Output *)
0
```

The discriminant can be represented in terms of roots as $D(f)==a_{n}^{2 n-2}\prod_{1 \leq i<j \leq n}(r_{i}-r_{j})^{2}$:

```wolfram
x/.Solve[a x^2+b x+c==0,x]
(* Output *)
{(-b-Sqrt[b^2-4 a c])/(2 a),(-b+Sqrt[b^2-4 a c])/(2 a)}
```

```wolfram
Expand[a^(2 2-2)(%[[1]]-%[[2]])^2]
(* Output *)
b^2-4 a c
```

```wolfram
Discriminant[a x^2+b x+c,x]
(* Output *)
b^2-4 a c
```

Equation $(-1)^{n(n-1)/2}a_{n}D(f)==R(f,f')$ relates [Discriminant](https://reference.wolfram.com/language/ref/Discriminant.html) and [Resultant](https://reference.wolfram.com/language/ref/Resultant.html):

```wolfram
f=3x^7-5x+4;
```

```wolfram
Resultant[f,D[f,x],x]
(* Output *)
4720053663936
```

```wolfram
(-1)^(7 6/2)Coefficient[f,x,7]Discriminant[f,x]
(* Output *)
4720053663936
```

### Possible Issues

Using exact coefficients, this indicates no common root:

```wolfram
Discriminant[(x-1+10^(-17))(x-1),x]
(* Output *)
(1)/(10000000000000000000000000000000000)
```

With approximate coefficients, this does indicate a common root:

```wolfram
Discriminant[(x-1.+10^(-17))(x-1),x]
(* Output *)
0
```

in this case, using higher precision resolves the problem:

```wolfram
Discriminant[(x-1+10^(-17))(x-1),x]
(* Output *)
1.×10^-34
```

## Related Guides ▪Polynomial Systems ▪Polynomial Algebra ▪Finite Fields

## History Introduced in 2007 (6.0) | Updated in 2022 (13.2) ▪ 2023 (13.3)
