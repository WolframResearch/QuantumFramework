# Cyclotomic | [SpanFromLeft]

> [Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html)[*n*,*x*] — gives the `*n*`$^{th}$ cyclotomic polynomial in `*x*`.

## Details

The cyclotomic polynomial $C_{n}(x)$ of order $n$ is defined to be $\prod_{k}(x-e^{2 \pi i k/n})$, where the product runs over integers $k$ less than $n$ that are relatively prime to $n$.

## Examples

### Basic Examples

Give the fifth cyclotomic polynomial:

```wolfram
Cyclotomic[5,x]
(* Output *)
1+x+x^2+x^3+x^4
```

The roots are the primitive fifth roots of $1$:

```wolfram
Solve[Cyclotomic[5,x]==0,x]
(* Output *)
{{x->-(-1)^(1/5)},{x->(-1)^(2/5)},{x->-(-1)^(3/5)},{x->(-1)^(4/5)}}
```

### Scope

[](https://reference.wolfram.com/language/ref/.html) formatting:

```wolfram
Cyclotomic[n,x]//
(* Output *)
n
```

### Applications

Compare a factored polynomial with [Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html):

```wolfram
Factor[x^11-1]
(* Output *)
(-1+x) (1+x+x^2+x^3+x^4+x^5+x^6+x^7+x^8+x^9+x^10)
```

```wolfram
Cyclotomic[11,x]
(* Output *)
1+x+x^2+x^3+x^4+x^5+x^6+x^7+x^8+x^9+x^10
```

Values of successive cyclotomic polynomials at 1:

```wolfram
Table[Cyclotomic[n,1],{n,2,20}]
(* Output *)
{2,3,2,5,1,7,2,3,1,11,1,13,1,1,2,17,1,19,1}
```

Calculate [unique primes](http://mathworld.wolfram.com/UniquePrime.html) for which the decimal expansion of $1/p$ has a unique period:

```wolfram
up[n_]:=If[Length[#]===1,#[[1,1]]]&[FactorInteger[Cyclotomic[n,10]/GCD[Cyclotomic[n,10],n]]]
```

```wolfram
Select[Table[up[n],{n,30}],NumberQ]
(* Output *)
{3,11,37,101,333667,9091,9901,909091,1111111111111111111,11111111111111111111111,99990001}
```

```wolfram
RealDigits[1/%]//Column
(* Output *)
{{{{{3}},0}}, {{{{9,0}},-1}}, {{{{2,7,0}},-1}}, {{{{9,9,0,0}},-2}}, {{{{2,9,9,7,0,0,0,0,0}},-5}}, {{{{1,0,9,9,9,8,9,0,0,0}},-3}}, {{{{1,0,0,9,9,9,8,9,9,0,0,0}},-3}}, {{{{1,0,9,9,9,9,9,8,9,0,0,0,0,0}},-5}}, {{{{9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},-18}}, {{{{9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},-22}}, {{{{1,0,0,0,0,9,9,9,9,9,9,9,8,9,9,9,9,0,0,0,0,0,0,0}},-7}}}
```

Plot cyclotomic polynomials:

```wolfram
Plot[Evaluate[Table[Cyclotomic[k,z],{k,0,10}]],{z,-1,1}]
```

*([Graphics])*

Plot arguments of roots of [Cyclotomic](https://reference.wolfram.com/language/ref/Cyclotomic.html):

```wolfram
ListPlot[Table[{Arg[#],n}&/@(z/.NSolve[Cyclotomic[n,z]==0,z]),{n,50}],PlotMarkers->None]
```

*([Graphics])*

Plot the degree and number of terms of cyclotomic polynomials:

```wolfram
{ListPlot[Table[Exponent[Cyclotomic[n,x],x],{n,500}]],
ListPlot[Table[Length[Cyclotomic[n,x]],{n,500}]]}
(* Output *)
{[Graphics],[Graphics]}
```

### Properties & Relations

Integrate a cyclotomic polynomial:

```wolfram
Integrate[1/Cyclotomic[3,x],x]
(* Output *)
(2 ArcTan[(1+2 x)/(Sqrt[3])])/(Sqrt[3])
```

Factor a cyclotomic polynomial over an extension field:

```wolfram
Factor[Cyclotomic[5,x]]
(* Output *)
1+x+x^2+x^3+x^4
```

```wolfram
Factor[Cyclotomic[5,x],Extension->{(-1)^(1/5)}]
(* Output *)
((-1)^(2/5)-x) (-1+(-1)^(1/5)-(-1)^(2/5)+(-1)^(3/5)-x) ((-1)^(1/5)+x) ((-1)^(3/5)+x)
```

Generate cyclotomic polynomials from a definition:

```wolfram
c[n_,z_]:=Collect[Product[z-ℯ^((2 π ⅈ k)/(n)),{k,Select[Range[n],CoprimeQ[#,n]&]}],z,FullSimplify]
```

```wolfram
Table[c[n,z],{n,10}]
(* Output *)
{-1+z,1+z,1+z+z^2,1+z^2,1+z+z^2+z^3+z^4,1-z+z^2,1+z+z^2+z^3+z^4+z^5+z^6,1+z^4,1+z^3+z^6,1-z+z^2-z^3+z^4}
```

Use an alternative definition, valid for $n>2$:

```wolfram
c[n_,z_]:=Cancel[Times@@((1-z^#)^MoebiusMu[n/#]&/@Divisors[n])]
```

```wolfram
Table[c[n,z],{n,6}]
(* Output *)
{1-z,1+z,1+z+z^2,1+z^2,1+z+z^2+z^3+z^4,1-z+z^2}
```

```wolfram
Table[Cyclotomic[n,z],{n,6}]
(* Output *)
{-1+z,1+z,1+z+z^2,1+z^2,1+z+z^2+z^3+z^4,1-z+z^2}
```

Form products of cyclotomic polynomials:

```wolfram
f[n_, z_]:=Expand[Times@@(Cyclotomic[#,z]&/@Divisors[n])]
```

```wolfram
Table[f[n,z],{n,1,8}]
(* Output *)
{-1+z,-1+z^2,-1+z^3,-1+z^4,-1+z^5,-1+z^6,-1+z^7,-1+z^8}
```

Plot the Riemann surface of an inverse of a cyclotomic polynomial over the complex plane:

```wolfram
ParametricPlot3D[Evaluate[{Re[Cyclotomic[10,x+I y]],Im[Cyclotomic[10,x+I y]],x}],{x,-2,2},{y,-2,2},Mesh->False,BoxRatios->1]
```

*([Graphics3D])*

Plot the complex roots of successive derivatives of the 50$^{th}$ cyclotomic polynomial:

```wolfram
Graphics[Point[Flatten[Map[{Re[#],Im[#]}&,Table[z/.NSolve[D[Cyclotomic[50,z],{z,k}]==0,z],{k,EulerPhi[50]-2}],{-1}],1]]]
```

*([Graphics])*

### Neat Examples

The first cyclotomic polynomial with a coefficient other than 0, ±1:

```wolfram
Cyclotomic[105,x]
(* Output *)
1+x+x^2-x^5-x^6-2 x^7-x^8-x^9+x^12+x^13+x^14+x^15+x^16+x^17-x^20-x^22-x^24-x^26-x^28+x^31+x^32+x^33+x^34+x^35+x^36-x^39-x^40-2 x^41-x^42-x^43+x^46+x^47+x^48
```

Nonzero coefficients of successive cyclotomic polynomials:

```wolfram
ArrayPlot[CoefficientList[Array[Cyclotomic[#,x]&,40],x]]
```

![image](img/image_001.png)

## Tech Notes ▪Algebraic Operations on Polynomials

## Related Guides ▪Polynomial Algebra

## Related Links [MathWorld](http://mathworld.wolfram.com/CyclotomicPolynomial.html) [The Wolfram Functions Site](http://functions.wolfram.com/Polynomials/Cyclotomic/) [NKS|Online](http://www.wolframscience.com/nks/search/?q=Cyclotomic) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0)
