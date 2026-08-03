# DiscreteConvolve | [SpanFromLeft]

> [DiscreteConvolve](https://reference.wolfram.com/language/ref/DiscreteConvolve.html)[*f*,*g*,*n*,*m*]  — gives the convolution with respect to `*n*` of the expressions `*f*` and `*g*`.
> [DiscreteConvolve](https://reference.wolfram.com/language/ref/DiscreteConvolve.html)[*f*,*g*,{*n*_1,*n*_2,…},{*m*_1,*m*_2,…}] — gives the multidimensional convolution.

## Details and Options

The convolution $(\mathit{f}★\mathit{g})(\mathit{m})$ of two sequences $f(n)$ and $g(n)$ is given by $\sum_{n=-\infty}^{\infty}f(n) g(m-n)$.

The multidimensional convolution is given by $\sum_{n_{1}=-\infty}^{\infty}\sum_{n_{2}=-\infty}^{\infty}\cdots \mathit{f}(n_{1},n_{2},\ldots) g(m_{1}-n_{1},m_{2}-n_{2},\ldots)$.

The following options can be given:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions to make about parameters |
| --- | --- | --- |
| [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to generate conditions on parameters |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | method to use |
| [VerifyConvergence](https://reference.wolfram.com/language/ref/VerifyConvergence.html) | [True](https://reference.wolfram.com/language/ref/True.html) | whether to verify convergence |

## Examples

### Basic Examples

Convolve a sequence with [DiscreteDelta](https://reference.wolfram.com/language/ref/DiscreteDelta.html):

```wolfram
DiscreteConvolve[DiscreteDelta[n],Sin[n],n,m]
(* Output *)
Sin[m]
```

Convolve two exponential sequences:

```wolfram
DiscreteConvolve[(1/2)^n UnitStep[n],(1/2)^n UnitStep[n],n,m]
(* Output *)
{, {{2^-m (1+m), m>=0}, {0, True}}}
```

Convolve two [UnitBox](https://reference.wolfram.com/language/ref/UnitBox.html) sequences and plot the result:

```wolfram
DiscreteConvolve[UnitBox[n/5],UnitBox[n/5],n,m];
```

```wolfram
DiscretePlot[%,{m,-8,8},PlotRange->All]
```

*([Graphics])*

### Scope

#### Univariate Convolution

Convolution sums the product of translates:

```wolfram
f[n_]:=UnitBox[n/7]
```

```wolfram
g[n_]:=2 UnitBox[n/3+1]
```

```wolfram
GraphicsRow@Table[DiscretePlot[{f[x],g[x-i]},{x,-10,10},ImageSize->110,Ticks->None,
PlotRange -> All], {i, -1, 2}]
(* Output *)
![image](img/image_001.png)
```

```wolfram
DiscreteConvolve[f[n], g[n],n,m]
(* Output *)
2 ({, {{1, m==-7||m==1}, {2, m==-6||m==0}, {3, -6<m<0}, {0, True}}})
```

```wolfram
DiscretePlot[%,{m,-10,10}]
```

*([Graphics])*

Convolution of elementary sequences:

```wolfram
DiscreteConvolve[1/(2n+1)^2, 1/(3n-1)^2, n, m]
(* Output *)
(π (π+6 m π+(1+6 m) π Sec[((1)/(6)+m) π]^2-12 Tan[((1)/(6)+m) π]))/((1+6 m)^3)
```

Convolution of piecewise sequences:

```wolfram
DiscreteConvolve[4^(-n) UnitStep[n],7^n UnitStep[n] , n, m]
(* Output *)
{, {{(1)/(27) 4^-m (-1+28^(1+m)), m>=0}, {0, True}}}
```

#### Multivariate Convolution

Perform discrete convolution with a unit box:

```wolfram
DiscreteConvolve[UnitBox[m/4,n/7+1],m^2 +n,{m,n},{r,s}]
(* Output *)
35 (2+r^2)+35 (7+s)
```

### Generalizations & Extensions

Multiplication by [UnitStep](https://reference.wolfram.com/language/ref/UnitStep.html) effectively gives the convolution over a finite interval:

```wolfram
DiscreteConvolve[2^n  UnitStep[n],n UnitStep[n],n,m]
(* Output *)
{, {{-2+2^(1+m)-m, m>=0}, {0, True}}}
```

```wolfram
Sum[2^n (m-n),{n,0,m}]
(* Output *)
-2+2^(1+m)-m
```

### Options

#### Assumptions

Specify assumptions on a variable or parameter:

```wolfram
DiscreteConvolve[UnitStep[n]/3^n, 2^n UnitStep[-n], n, m,
Assumptions -> m>0]
(* Output *)
(2 3^(1-m))/(5)
```

```wolfram
DiscreteConvolve[UnitStep[n]/3^n, 2^n UnitStep[-n], n, m,
Assumptions -> m<0]
(* Output *)
(3 2^(1+m))/(5)
```

#### GenerateConditions

Generate conditions for the range of a parameter:

```wolfram
DiscreteConvolve[E^(-a Abs[n]), 1, n, m, GenerateConditions -> True]
(* Output *)
(1+ℯ^a)/(-1+ℯ^a)
```

### Applications

Obtain a particular solution for a linear difference equation:

```wolfram
y[n]/.RSolve[{y[n]-2y[n-1]==n,y[0]==0}, y,n][[1]]
(* Output *)
-2+2^(1+n)-n
```

```wolfram
DiscreteConvolve[2^k UnitStep[k],k UnitStep[k],k,n,Assumptions->n>0]
(* Output *)
-2+2^(1+n)-n
```

Obtain the step response of a linear, time-invariant system given its impulse response `*h*`:

```wolfram
h=(4/5)^n UnitStep[n];
```

```wolfram
DiscretePlot[h,{n,0,20}]
```

*([Graphics])*

The step response corresponding to this system:

```wolfram
DiscreteConvolve[h,UnitStep[n],n,m]
(* Output *)
{, {{5-4^(1+m) 5^-m, m>=0}, {0, True}}}
```

```wolfram
DiscretePlot[%,{m,-5,20}]
```

*([Graphics])*

### Properties & Relations

[DiscreteConvolve](https://reference.wolfram.com/language/ref/DiscreteConvolve.html) computes a sum over the set of integers:

```wolfram
DiscreteConvolve[1/(3n+1)^2,UnitStep[n],n,m,Assumptions -> m >0]
(* Output *)
(1)/(27) (4 π^2-3 PolyGamma[1,(4)/(3)+m])
```

```wolfram
Sum[UnitStep[m-n]/(3n+1)^2,{n,-Infinity,Infinity},
  Assumptions->Element[m,Integers]&&m > 0]
(* Output *)
(1)/(27) (4 π^2-3 PolyGamma[1,(4)/(3)+m])
```

Convolution with [DiscreteDelta](https://reference.wolfram.com/language/ref/DiscreteDelta.html) gives the value of a sequence at `*m*`:

```wolfram
DiscreteConvolve[f[n],DiscreteDelta[n],n,m]
(* Output *)
f[m]
```

Scaling:

```wolfram
DiscreteConvolve[a r[n],s[n],n,m]
(* Output *)
a DiscreteConvolve[r[n],s[n],n,m]
```

```wolfram
DiscreteConvolve[r[n],a s[n],n,m]
(* Output *)
a DiscreteConvolve[r[n],s[n],n,m]
```

Commutativity:

```wolfram
f[n_]:=2^n UnitStep[n]
```

```wolfram
g[n_]:=UnitStep[n]
```

```wolfram
DiscreteConvolve[f[n],g[n],n,m]
(* Output *)
{, {{-1+2^(1+m), m>=0}, {0, True}}}
```

```wolfram
DiscreteConvolve[g[n],f[n],n,m]
(* Output *)
{, {{-1+2^(1+m), m>=0}, {0, True}}}
```

Distributivity:

```wolfram
DiscreteConvolve[h[n]+i[n],j[n],n,m]
(* Output *)
DiscreteConvolve[h[n],j[n],n,m]+DiscreteConvolve[i[n],j[n],n,m]
```

```wolfram
DiscreteConvolve[h[n],i[n]+j[n],n,m]
(* Output *)
DiscreteConvolve[h[n],i[n],n,m]+DiscreteConvolve[h[n],j[n],n,m]
```

The [ZTransform](https://reference.wolfram.com/language/ref/ZTransform.html) of a causal convolution is the product of the individual transforms:

```wolfram
f[n_]:=2^n UnitStep[n]
```

```wolfram
g[n_]:=UnitStep[n]
```

```wolfram
ZTransform[DiscreteConvolve[f[n],g[n],n,m],m,z]
(* Output *)
(z^2)/(2-3 z+z^2)
```

```wolfram
Simplify[ZTransform[f[n],n,z]ZTransform[g[n],n,z]]
(* Output *)
(z^2)/(2-3 z+z^2)
```

Similarly for [GeneratingFunction](https://reference.wolfram.com/language/ref/GeneratingFunction.html):

```wolfram
GeneratingFunction[DiscreteConvolve[f[n],g[n],n,m],m,z]
(* Output *)
(1)/(1-3 z+2 z^2)
```

```wolfram
GeneratingFunction[f[n],n,z]GeneratingFunction[g[n],n,z]
(* Output *)
(1)/((1-2 z) (1-z))
```

```wolfram
%%-%//Simplify
(* Output *)
0
```

The [FourierSequenceTransform](https://reference.wolfram.com/language/ref/FourierSequenceTransform.html) of a convolution is the product of the individual transforms:

```wolfram
f[n_]:=3^(-n) UnitStep[n]
```

```wolfram
g[n_]:=2^(n)UnitStep[-n]
```

```wolfram
DiscreteConvolve[f[n],g[n],n,m]
(* Output *)
2^m ({, {{(6)/(5), m<=0}, {(6^(1-m))/(5), True}}})
```

```wolfram
Together[FourierSequenceTransform[%, m, w]]
(* Output *)
-(6 ℯ^(ⅈ w))/(2-7 ℯ^(ⅈ w)+3 ℯ^(2 ⅈ w))
```

```wolfram
Simplify[FourierSequenceTransform[f[n], n, w]
 FourierSequenceTransform[g[n], n, w]]
(* Output *)
-(6 ℯ^(ⅈ w))/(2-7 ℯ^(ⅈ w)+3 ℯ^(2 ⅈ w))
```

### Interactive Examples

This demonstrates the discrete-time convolution operation $(f* g) (y)$:

```wolfram
DynamicModule[{f,g,h},
f[n_]:=(3/4)^n UnitStep[n];
g[n_]:=UnitBox[n/5];
h[n_]=DiscreteConvolve[f[m],g[m],m,n];
Manipulate[GraphicsDiscretePlot[{f[m],g[n-m], PlotLabel -> Text[Style[f[m], g[n-m], Italic]], Epilog -> Arrowheads[Medium]Arrow[n-1n-0.1]Text[Style[n, Italic], n+1.5-0.9]],
DiscretePlot[f[m]g[n-m],{m,-5,12},PlotRange -> -11.5, Ticks -> Automatic01, AxesLabel -> Automatic, PlotStyle -> PointSize[0.04], PlotLabel -> Text[Style[f[m]·g[n-m], Italic]]],
DiscretePlot[h[m],{m,-5,n},PlotRange -> -512-14, Ticks -> Automatic024, AxesLabel -> Text[Style[n, Italic]]None, PlotStyle -> PointSize[0.04], PlotLabel -> Text[Style[(f*g)[n], Italic]]]},ImageSize->480],{{n,1},-4,11,1}]]
```

## Related Guides ▪Summation Transforms ▪Fourier Analysis ▪Additive Number Theory

## Related Links [MathWorld](https://mathworld.wolfram.com/Convolution.html)

## History Introduced in 2008 (7.0)
