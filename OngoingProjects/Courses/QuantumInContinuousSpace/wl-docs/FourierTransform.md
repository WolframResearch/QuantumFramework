# FourierTransform | [SpanFromLeft]

> [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html)[*f[t]*,*t*,*ω*] — gives the symbolic Fourier transform of `*f[t]*` in the variable `*t*` as `*F[ω]*` in the variable `*ω*`.
> [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html)[*f[t]*,*t*,*ω*^{^}] — gives the numeric Fourier transform at the numerical value `*ω*^{^}`.
> [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html)[*f*[*t*_1,…,*t*_n],{*t*_1,…,*t*_n},{*ω*_1,…,*ω*_n}] — gives the multidimensional Fourier transform of `*f*[*t*_1,…,*t*_n]`.

## Details and Options

The Fourier transform and its inverse are a way to transform between the time domain and the frequency domain.

Fourier transforms are typically used to reduce ordinary and partial differential equations to algebraic or ordinary differential equations, respectively. They are also used extensively in control theory and signal processing. Finally, they have applications in studying quantum mechanical phenomena, noise filtering, etc.

The Fourier transform of the time domain function $f(t)$ is the frequency domain function $F(\omega)$:

$$
\begin{pmatrix}
\text{[Graphics]} & \text{[Graphics]} & \text{[Graphics]} \\
"time domain" &  & "frequency domain"
\end{pmatrix}
$$

The Fourier transform of a function $f(t)$ is by default defined to be $\frac{1}{\sqrt{2 \pi}}\int_{-\infty}^{\infty}f(t)e^{i \omega t}d t$.

The multidimensional Fourier transform of a function $f(t_{1},t_{2},\ldots,t_{n})$ is by default defined to be $\frac{1}{(2 \pi)^{n/2}}\int_{-\infty}^{\infty}\cdots \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} f(t_{1},t_{2},\ldots,t_{n})e^{i (t_{1}\omega_{1}+t_{2}\omega_{2}+\ldots+t_{n}\omega_{n})}\mathrm{d}t_{1}\mathrm{d}t_{2}\ldots \mathrm{d}t_{n}$ or when using vector notation $\frac{1}{(2 \pi)^{n/2}}\int_{t \in \mathbb{R}^{n}}f(t) e^{i \omega.t}\mathrm{d}t$.

Different choices of definitions can be specified using the option [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html).

The integral is computed using numerical methods if the third argument, $\omega$, is given a numerical value.

The asymptotic Fourier transform can be computed using [Asymptotic](https://reference.wolfram.com/language/ref/Asymptotic.html).

There are several related Fourier transformations:

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) | infinite continuous-time functions (FT)
[FourierSequenceTransform](https://reference.wolfram.com/language/ref/FourierSequenceTransform.html) | infinite discrete-time functions (DTFT)
[FourierCoefficient](https://reference.wolfram.com/language/ref/FourierCoefficient.html) | finite continuous-time functions (FS)
[Fourier](https://reference.wolfram.com/language/ref/Fourier.html) | finite discrete-time functions (DFT)

The Fourier transform is an automorphism in the Schwartz vector space of functions whose derivatives are rapidly decreasing and thus induces an automorphism in its dual: the space of tempered distributions. These include absolutely integrable functions, well-behaved functions of polynomial growth and compactly supported distributions.

Hence, [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) not only works with absolutely integrable functions, but it can also handle a variety of tempered distributions such as [DiracDelta](https://reference.wolfram.com/language/ref/DiracDelta.html) to enlarge the pool of functions or generalized functions it can effectively transform.

The following options can be given:

| [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | digits of absolute accuracy sought |
| --- | --- | --- |
| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions to make about parameters |
| [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) | {0,1} | parameters to define the Fourier transform |
| [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to generate answers that involve conditions on parameters |
| [PerformanceGoal](https://reference.wolfram.com/language/ref/PerformanceGoal.html) | [$PerformanceGoal](https://reference.wolfram.com/language/ref/$PerformanceGoal.html) | aspects of performance to optimize |
| [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | digits of precision sought |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the precision used in internal computations |

Common settings for [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) include:

| `{0,1}` | $\frac{1}{\sqrt{2 \pi}}\int_{-\infty}^{\infty}f(t)e^{i \omega t}d t$ | default setting/physics |
| --- | --- | --- |
| `{1,-1}` | $\int_{-\infty}^{\infty}f(t)e^{-i \omega t}d t$ | systems engineering/mathematics |
| `{-1,1}` | $\frac{1}{2 \pi}\int_{-\infty}^{\infty}f(t)e^{i \omega t}d t$ | classical physics |
| `{0,-2[Pi](https://reference.wolfram.com/language/ref/Pi.html)}` | $\int_{-\infty}^{\infty}f(t)e^{-i 2 \pi \omega t}d t$ | ordinary frequency |
| `{*a*,*b*}` | $\sqrt{\frac{|b|}{(2 \pi)^{1-a}}}\int_{-\infty}^{\infty}f(t) e^{i b \omega t}d t$ | general setting |

In [](https://reference.wolfram.com/language/ref/.html), [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) is output using `ℱ`.

## Examples

### Basic Examples

Compute the Fourier transform of a function:

```wolfram
FourierTransform[UnitBox[t/2],t,ω]
(* Output *)
Sqrt[(2)/(π)] Sinc[ω]
```

Plot the function and its Fourier transform:

```wolfram
{Plot[UnitBox[t/2],{t,-2,2},Exclusions->None],Plot[%,{ω,-10,10}]}
(* Output *)
{[Graphics],[Graphics]}
```

Fourier transform of $\delta(t)$:

```wolfram
FourierTransform[DiracDelta[t],t,ω]
(* Output *)
(1)/(Sqrt[2 π])
```

For the systems engineering convention, change the parameters:

```wolfram
FourierTransform[DiracDelta[t],t,ω,FourierParameters->{1,-1}]
(* Output *)
1
```

The Fourier transform of a Gaussian is another Gaussian:

```wolfram
FourierTransform[E^(-t^2),t,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

Plot both Gaussians:

```wolfram
{Plot[E^(-t^2),{t,-5,5},PlotRange->Full],Plot[%,{ω,-5,5},PlotRange->Full]}
(* Output *)
{[Graphics],[Graphics]}
```

Compute the Fourier transform of a multivariate function:

```wolfram
FourierTransform[1/Sqrt[x^2+y^2],{x,y},{u,v}]
(* Output *)
(1)/(Sqrt[u^2+v^2])
```

Plot the result:

```wolfram
Plot3D[%,{u,-2,2},{v,-2,2},PlotRange->{0,10},Mesh->None]
```

*([Graphics3D])*

Compute the transform at a single point:

```wolfram
FourierTransform[Exp[-t^6] Sin[t^2],t,0.3]
(* Output *)
0.21026103935895005-5.984887656489812×10^-15 ⅈ
```

### Scope

#### Basic Uses

Fourier transform of a function for a symbolic parameter $\omega$:

```wolfram
FourierTransform[HeavisideTheta[t],t,ω]
(* Output *)
(ⅈ)/(Sqrt[2 π] ω)+Sqrt[(π)/(2)] DiracDelta[ω]
```

Fourier transforms of trigonometric functions:

```wolfram
FourierTransform[Sin[α t],t,ω]
(* Output *)
ⅈ Sqrt[(π)/(2)] DiracDelta[-α+ω]-ⅈ Sqrt[(π)/(2)] DiracDelta[α+ω]
```

```wolfram
FourierTransform[Cos[α t],t,ω]
(* Output *)
Sqrt[(π)/(2)] DiracDelta[-α+ω]+Sqrt[(π)/(2)] DiracDelta[α+ω]
```

Evaluate the Fourier transform for a numerical value of the parameter $\omega$:

```wolfram
FourierTransform[Exp[-t^2+t],t,0.3]
(* Output *)
0.8777740781424105+0.13266257670596573 ⅈ
```

[](https://reference.wolfram.com/language/ref/.html) formatting:

```wolfram
FourierTransform[f[t],t,ω]//
(* Output *)
ℱ_t[f(t)](ω)
```

#### Elementary Functions

Fourier transform of a power function:

```wolfram
FourierTransform[t^2,t,ω]
(* Output *)
-Sqrt[2 π] DiracDelta^′′[ω]
```

Polynomial:

```wolfram
FourierTransform[3t^2+5t-7,t,ω]
(* Output *)
-7 Sqrt[2 π] DiracDelta[ω]-5 ⅈ Sqrt[2 π] DiracDelta^′[ω]-3 Sqrt[2 π] DiracDelta^′′[ω]
```

Fourier transform of rational functions:

```wolfram
FourierTransform[1/(1+t^2),t,ω]
(* Output *)
ℯ^(-Abs[ω]) Sqrt[(π)/(2)]
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

```wolfram
FourierTransform[t/(1+t^3),t,ω]
(* Output *)
(1)/(3) Sqrt[(π)/(2)] (-2 (-1)^(5/6) ℯ^((1)/(2) (ⅈ+Sqrt[3]) ω) HeavisideTheta[-ω]+(6 ⅈ ℯ^((-1)^(5/6) ω) HeavisideTheta[ω])/((1+(-1)^(1/3))^2)-ⅈ ℯ^(-ⅈ ω) Sign[ω])
```

Plot the real and imaginary parts:

```wolfram
ReImPlot[%,{ω,-4,4}]
```

*([Graphics])*

Reciprocal of square root:

```wolfram
FourierTransform[1/Sqrt[Abs[t]],t,ω]
(* Output *)
(1)/(Sqrt[Abs[ω]])
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

Expressions involving trigonometric functions:

```wolfram
FourierTransform[Cos[t]^2,t,ω]
(* Output *)
(1)/(2) Sqrt[(π)/(2)] DiracDelta[-2+ω]+Sqrt[(π)/(2)] DiracDelta[ω]+(1)/(2) Sqrt[(π)/(2)] DiracDelta[2+ω]
```

```wolfram
FourierTransform[Exp[α t] Sin[ β t],t,ω]
(* Output *)
ⅈ Sqrt[(π)/(2)] DiracDelta[ⅈ α+β-ω]-ⅈ Sqrt[(π)/(2)] DiracDelta[-ⅈ α+β+ω]
```

```wolfram
FourierTransform[Exp[-t^2] Sin[t],t,ω]
(* Output *)
(ⅈ (-1+Cosh[ω]+Sinh[ω]) (Cosh[(1)/(4) (1+ω)^2]-Sinh[(1)/(4) (1+ω)^2]))/(2 Sqrt[2])
```

```wolfram
FourierTransform[(1)/(2α)(Sin[α t]+α t Cos[α t]),t,ω]
(* Output *)
(ⅈ Sqrt[(π)/(2)] DiracDelta[-α+ω])/(2 α)-(ⅈ Sqrt[(π)/(2)] DiracDelta[α+ω])/(2 α)-(1)/(2) ⅈ Sqrt[(π)/(2)] DiracDelta^′[-α+ω]-(1)/(2) ⅈ Sqrt[(π)/(2)] DiracDelta^′[α+ω]
```

```wolfram
FourierTransform[Sin[t]+Cos[t],t,ω]
(* Output *)
(1+ⅈ) Sqrt[(π)/(2)] DiracDelta[-1+ω]+(1-ⅈ) Sqrt[(π)/(2)] DiracDelta[1+ω]
```

Ratio of sine and linear function:

```wolfram
FourierTransform[(Sin[3 t])/(t),t,ω]
(* Output *)
(1)/(2) Sqrt[(π)/(2)] (Sign[3-ω]+Sign[3+ω])
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5},Exclusions->None]
```

*([Graphics])*

Composition of elementary functions:

```wolfram
FourierTransform[Sin[t^2],t,ω]
(* Output *)
(1)/(2) (Cos[(ω^2)/(4)]-Sin[(ω^2)/(4)])
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

Logarithmic function:

```wolfram
FourierTransform[Log[Abs[t]],t,ω]
(* Output *)
-(Sqrt[(π)/(2)])/(Abs[ω])-EulerGamma Sqrt[2 π] DiracDelta[ω]
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

#### Special Functions

[Sinc](https://reference.wolfram.com/language/ref/Sinc.html) function:

```wolfram
FourierTransform[Sinc[t],t,ω]
(* Output *)
(1)/(2) Sqrt[(π)/(2)] (Sign[1-ω]+Sign[1+ω])
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5},Exclusions->None]
```

*([Graphics])*

Expressions involving Bessel functions:

```wolfram
FourierTransform[BesselJ[2,t],t,ω]
(* Output *)
-(Sqrt[(2)/(π)] (-1+2 ω^2) (-HeavisideTheta[-1+ω]+HeavisideTheta[1+ω]))/(Sqrt[1-ω^2])
```

Plot the transform:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

[SinIntegral](https://reference.wolfram.com/language/ref/SinIntegral.html) function:

```wolfram
FourierTransform[SinIntegral[t],t,ω,Assumptions->-1<ω<1]
(* Output *)
(ⅈ Sqrt[(π)/(2)])/(ω)
```

Laguerre polynomial:

```wolfram
FourierTransform[LaguerreL[5,t],t,ω]
(* Output *)
Sqrt[2 π] DiracDelta[ω]+5 ⅈ Sqrt[2 π] DiracDelta^′[ω]-5 Sqrt[2 π] DiracDelta^′′[ω]-(5)/(3) ⅈ Sqrt[2 π] DiracDelta^(3)[ω]+(5)/(12) Sqrt[(π)/(2)] DiracDelta^(4)[ω]+(1)/(60) ⅈ Sqrt[(π)/(2)] DiracDelta^(5)[ω]
```

Airy function:

```wolfram
FourierTransform[AiryAi[t],t,ω]
(* Output *)
(ℯ^(-(ⅈ ω^3)/(3)))/(Sqrt[2 π])
```

Plot the magnitude of the Fourier transform for complex $\omega$:

```wolfram
Plot3D[Abs[%]/.ω->a+b ⅈ,{a,-5,5},{b,-5,5}]
(* Output *)
![image](img/image_001.png)
```

#### Piecewise Functions and Distributions

Fourier transform of a piecewise function:

```wolfram
f[t_]=Piecewise[{{t,0<= t<= 1},{1,t>1}}];
Plot[f[t],{t,0,3}]
```

*([Graphics])*

```wolfram
FourierTransform[f[t],t,ω]
(* Output *)
(-(1)/(ω^2)+(Cos[ω])/(ω^2)+π DiracDelta[ω]+(ⅈ Sinc[ω])/(ω))/(Sqrt[2 π])
```

Absolute value using [Sign](https://reference.wolfram.com/language/ref/Sign.html) function:

```wolfram
f[t_]=t Sign[t];
Plot[f[t],{t,-2,2}]
```

*([Graphics])*

```wolfram
FourierTransform[f[t],t,ω]
(* Output *)
-(Sqrt[(2)/(π)])/(ω^2)
```

Restriction of a sine function to a half-period:

```wolfram
Plot[Sin[α t]UnitBox[α t/π]/.α-> π,{t,-1,1}]
```

*([Graphics])*

```wolfram
Simplify[FourierTransform[Sin[α t]UnitBox[α t/π],t,ω],α>0]
(* Output *)
(ⅈ Sqrt[(π)/(2)] (Sinc[(π (α-ω))/(2 α)]-Sinc[(π (α+ω))/(2 α)]))/(2 α)
```

Triangular function:

```wolfram
f[t_]=Piecewise[{{t,0<= t <= 1},{2-t,1< t<= 2},{0,2<t }}];
Plot[f[t],{t,0,4}]
```

*([Graphics])*

```wolfram
FourierTransform[f[t],t,ω]
(* Output *)
-((-1+ℯ^(ⅈ ω))^2)/(Sqrt[2 π] ω^2)
```

[Ramp](https://reference.wolfram.com/language/ref/Ramp.html):

```wolfram
FourierTransform[Ramp[t-α],t,ω]
(* Output *)
-(ℯ^(ⅈ α ω))/(Sqrt[2 π] ω^2)-Sqrt[(π)/(2)] α DiracDelta[ω]-ⅈ Sqrt[(π)/(2)] DiracDelta^′[ω]
```

[UnitStep](https://reference.wolfram.com/language/ref/UnitStep.html):

```wolfram
FourierTransform[UnitStep[t-α],t,ω]
(* Output *)
(ⅈ ℯ^(ⅈ α ω))/(Sqrt[2 π] ω)+Sqrt[(π)/(2)] DiracDelta[ω]
```

Product of [UnitStep](https://reference.wolfram.com/language/ref/UnitStep.html) and cosine functions:

```wolfram
FourierTransform[UnitStep[t-π]Cos[t-π],t,ω]
(* Output *)
(ⅈ ℯ^(ⅈ π ω))/(2 Sqrt[2 π] (-1+ω))+(ⅈ ℯ^(ⅈ π ω))/(2 Sqrt[2 π] (1+ω))-(1)/(2) Sqrt[(π)/(2)] DiracDelta[-1+ω]-(1)/(2) Sqrt[(π)/(2)] DiracDelta[1+ω]
```

Plot the magnitude and phase:

```wolfram
GraphicsPlot[Abs[%],{ω,-5,5},PlotRange -> -0.55., AxesLabel -> ωNone, PlotLabel -> magnitude, ImageSize -> 200],Plot[Arg[%],{ω,-5,5},Exclusions -> None, PlotRange -> 1, /, 2, ), -1010, AxesLabel -> ωNone, PlotLabel -> phase, ImageSize -> 200]
```

*([Graphics])*

#### Periodic Functions

Fourier transform of [SquareWave](https://reference.wolfram.com/language/ref/SquareWave.html):

```wolfram
Plot[SquareWave[t],{t,0,4},Exclusions->None]
```

*([Graphics])*

```wolfram
FourierTransform[SquareWave[t],t,ω]
(* Output *)
∑_{K[1]=-∞}^{∞}Sqrt[2 π] DiracDelta[ω+π K[1]] ({, {{0, K[1]==0}, {-(ⅈ (-1)^K[1] (1+(-1)^K[1]) (-1+ⅈ^K[1])^2)/(2 π K[1]), True}}})
```

[TriangleWave](https://reference.wolfram.com/language/ref/TriangleWave.html):

```wolfram
Plot[TriangleWave[t],{t,0,4}]
```

*([Graphics])*

```wolfram
FourierTransform[TriangleWave[t],t,ω]
(* Output *)
∑_{K[1]=-∞}^{∞}Sqrt[2 π] DiracDelta[ω+π K[1]] ({, {{0, K[1]==0}, {(1)/(2 π^2)ⅈ ((2 π (4 (-1)^Abs[K[1]]+Cos[(1)/(4) π Abs[K[1]]]-6 Cos[(3)/(4) π Abs[K[1]]]) Sign[K[1]])/(Abs[K[1]])-(1)/(K[1]^2)(-1)^K[1] (8 π K[1]-6 ℯ^((1)/(4) ⅈ π K[1]) π K[1]-6 ℯ^((7)/(4) ⅈ π K[1]) π K[1]+ℯ^((5)/(4) ⅈ π K[1]) (-4 ⅈ+π K[1])+ℯ^((3)/(4) ⅈ π K[1]) (4 ⅈ+π K[1])+8 (-1)^K[1] Sign[K[1]] (Sin[(1)/(4) π Abs[K[1]]]-2 Sin[(3)/(4) π Abs[K[1]]]))), True}}})
```

[SawtoothWave](https://reference.wolfram.com/language/ref/SawtoothWave.html):

```wolfram
Plot[SawtoothWave[t],{t,0,6},Exclusions->None]
```

*([Graphics])*

```wolfram
FourierTransform[SawtoothWave[t],t,ω]
(* Output *)
∑_{K[1]=-∞}^{∞}Sqrt[2 π] DiracDelta[ω+π K[1]] ({, {{(1)/(2), K[1]==0}, {-((-1)^K[1] (1+(-1)^K[1]) (-1+(-1)^K[1]-ⅈ π K[1]))/(2 π^2 K[1]^2), True}}})
```

Full-wave-rectified function with period $\pi$:

```wolfram
Plot[RealAbs[Sin[t]],{t,-2π,2π}]
```

*([Graphics])*

```wolfram
FourierTransform[RealAbs[Sin[t]],t,ω]
(* Output *)
∑_{K[1]=-∞}^{∞}-(2 (-1)^K[1] (1+(-1)^K[1]) (1+ⅈ^K[1])^2 Sqrt[(2)/(π)] DiracDelta[2 ω+K[1]])/(-4+K[1]^2)
```

Rectified wave:

```wolfram
f[t_]=Piecewise[{{Sin[t],0<= Mod[(t)/(π),2]<= 1},{0,1< Mod[(t)/(π),2]<2}}];
Plot[f[t],{t,0,10}]
```

*([Graphics])*

```wolfram
FourierTransform[f[t],t,ω]
(* Output *)
∑_{K[1]=-∞}^{∞}-(2 ⅈ^(-K[1]) (1+(-1)^K[1]) (1+ⅈ^K[1]) Sqrt[(2)/(π)] DiracDelta[2 ω+K[1]])/(-4+K[1]^2)
```

#### Generalized Functions

Fourier transform involving [HeavisideTheta](https://reference.wolfram.com/language/ref/HeavisideTheta.html):

```wolfram
FourierTransform[HeavisideTheta[t],t,ω]
(* Output *)
(ⅈ)/(Sqrt[2 π] ω)+Sqrt[(π)/(2)] DiracDelta[ω]
```

```wolfram
FourierTransform[HeavisideTheta[t] HeavisideTheta[1-t],t,ω]
(* Output *)
(ⅈ-ⅈ Cos[ω]+Sin[ω])/(Sqrt[2 π] ω)
```

Plot the magnitude and phase:

```wolfram
GraphicsPlot[Abs[%],{ω,-10,10},PlotRange -> -0.51., AxesLabel -> ωNone, PlotLabel -> magnitude, ImageSize -> 200],Plot[Arg[%],{ω,-10,10},Exclusions -> None, PlotRange -> -55, AxesLabel -> ωNone, PlotLabel -> phase, ImageSize -> 200]
```

*([Graphics])*

[DiracDelta](https://reference.wolfram.com/language/ref/DiracDelta.html):

```wolfram
FourierTransform[DiracDelta[t],t,ω]
(* Output *)
(1)/(Sqrt[2 π])
```

Derivative of [DiracDelta](https://reference.wolfram.com/language/ref/DiracDelta.html):

```wolfram
FourierTransform[DiracDelta'[t],t,ω]
(* Output *)
-(ⅈ ω)/(Sqrt[2 π])
```

[HeavisideLambda](https://reference.wolfram.com/language/ref/HeavisideLambda.html):

```wolfram
Plot[HeavisideLambda[t-1],{t,0,4}]
```

*([Graphics])*

```wolfram
FourierTransform[HeavisideLambda[t-1],t,ω]
(* Output *)
-((-1+ℯ^(ⅈ ω))^2)/(Sqrt[2 π] ω^2)
```

[HeavisidePi](https://reference.wolfram.com/language/ref/HeavisidePi.html):

```wolfram
Plot[HeavisidePi[t-(3)/(2)],{t,0,3},Exclusions->None]
```

*([Graphics])*

```wolfram
FourierTransform[HeavisidePi[t-(3)/(2)],t,ω]
(* Output *)
-(ⅈ ℯ^(ⅈ ω) (-1+ℯ^(ⅈ ω)))/(Sqrt[2 π] ω)
```

#### Multivariate Functions

Bivariate Fourier transform of a constant:

```wolfram
FourierTransform[1,{x,y},{u,v}]
(* Output *)
2 π DiracDelta[u] DiracDelta[v]
```

Exponential function:

```wolfram
FourierTransform[Exp[α x+β y],{x,y},{u,v}]
(* Output *)
2 π DiracDelta[u-ⅈ α] DiracDelta[v-ⅈ β]
```

Trivariate cosine:

```wolfram
FourierTransform[Cos[x+y+z],{x,y,z},{u,v,w}]
(* Output *)
Sqrt[2] π^(3/2) DiracDelta[-1+u] DiracDelta[-1+v] DiracDelta[-1+w]+Sqrt[2] π^(3/2) DiracDelta[1+u] DiracDelta[1+v] DiracDelta[1+w]
```

Product of power and exponential:

```wolfram
FourierTransform[(x y z)^3 Exp[-(x^2+y^2+z^2)],{x,y,z},{u,v,w}]
(* Output *)
(ⅈ ℯ^(-(u^2)/(4)-(v^2)/(4)-(w^2)/(4)) u (-6+u^2) v (-6+v^2) w (-6+w^2))/(1024 Sqrt[2])
```

Fourier transform of a product of exponential and [SquareWave](https://reference.wolfram.com/language/ref/SquareWave.html) functions:

```wolfram
Plot3D[E^(-Abs[x])SquareWave[y/(2π)],{x,-2π,2π},{y,-2π,2π},PlotRange->All,Mesh->None]
```

*([Graphics3D])*

```wolfram
FourierTransform[E^(-Abs[x]) SquareWave[y/(2 Pi)],{x,y},{u,v}]
(* Output *)
(1)/(1+u^2)Sqrt[(2)/(π)] ∑_{K[1]=-∞}^{∞}2 Sqrt[2 π] DiracDelta[2 v+K[1]] ({, {{0, K[1]==0}, {-(ⅈ (-1)^K[1] (1+(-1)^K[1]) (-1+ⅈ^K[1])^2)/(2 π K[1]), True}}})
```

#### Formal Properties

Fourier transform of a first-order derivative:

```wolfram
FourierTransform[f'[t],t,ω]
(* Output *)
-ⅈ ω FourierTransform[f[t],t,ω]
```

Fourier transform of a second-order derivative:

```wolfram
FourierTransform[f''[t],t,ω]
(* Output *)
-ω^2 FourierTransform[f[t],t,ω]
```

Fourier transform threads itself over equations:

```wolfram
FourierTransform[f'[t]==Log[t],t,ω]
(* Output *)
-ⅈ ω FourierTransform[f[t],t,ω]==(Sqrt[(π)/(2)])/(ω)-(Sqrt[(π)/(2)])/(Abs[ω])+(ⅈ π^(3/2) DiracDelta[ω])/(Sqrt[2])-EulerGamma Sqrt[2 π] DiracDelta[ω]
```

#### Numerical Evaluation

Calculate the Fourier transform at a single point:

```wolfram
FourierTransform[(1)/(Sqrt[t]),t,-.9]
(* Output *)
1.0540925533894598-1.0540925533894598 ⅈ
```

Alternatively, calculate the Fourier transform symbolically:

```wolfram
FourierTransform[(1)/(Sqrt[t]),t,ω]
(* Output *)
-(((1)/(2)-(ⅈ)/(2)) (-1+Sign[ω]))/(Sqrt[Abs[ω]])
```

Then evaluate it for the specific value of $\omega$:

```wolfram
N[%/.ω-> -.9]
(* Output *)
1.0540925533894598-1.0540925533894598 ⅈ
```

### Options

#### AccuracyGoal

The option [AccuracyGoal](https://reference.wolfram.com/language/ref/AccuracyGoal.html) sets the number of digits of accuracy:

```wolfram
exact=FourierTransform[HeavisideLambda[t-1],t,-9/10]
(* Output *)
-(50)/(81) (-1+ℯ^(9 ⅈ/10))^2 ℯ^(-9 ⅈ/5) Sqrt[(2)/(π)]
```

```wolfram
FourierTransform[HeavisideLambda[t-1],t,-.9,AccuracyGoal->5]-exact
(* Output *)
5.546899278785489×10^-8+5.2648017057066454×10^-8 ⅈ
```

With default settings:

```wolfram
FourierTransform[HeavisideLambda[t-1],t,-.9]-exact
(* Output *)
3.1922080812041287×10^-9+1.335343013941781×10^-8 ⅈ
```

#### Assumptions

Specify the range of a variable using [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html):

```wolfram
FourierTransform[BesselJ[3,t],t,ω,Assumptions->-1<ω<1 ]
(* Output *)
-(ⅈ Sqrt[(2)/(π)] ω (-3+4 ω^2))/(Sqrt[1-ω^2])
```

```wolfram
FourierTransform[BesselJ[3,t],t,ω,Assumptions->ω>1 ]
(* Output *)
0
```

```wolfram
FourierTransform[BesselJ[3, t],t,ω,Assumptions->ω<-1 ]
(* Output *)
0
```

#### FourierParameters

Fourier transform for the unit box function with different parameters:

```wolfram
params={{0,1},{1,1},{-1,1},{0,2 π}};
funs=funs=Table[FourierTransform[UnitBox[t-(1)/(2)],t,ω,FourierParameters->p],{p,params}]
(* Output *)
{(ⅈ-ⅈ Cos[ω]+ω Sinc[ω])/(Sqrt[2 π] ω),(ⅈ-ⅈ Cos[ω]+ω Sinc[ω])/(ω),(ⅈ-ⅈ Cos[ω]+ω Sinc[ω])/(2 π ω),-(ⅈ (-1+ℯ^(2 ⅈ π ω)))/(2 π ω)}
```

Create a nicely formatted table of the results:

```wolfram
header={ "Parameters",HoldForm@FourierTransform[UnitBox[t-(1)/(2)],t,ω]};
Grid[Prepend[Transpose[{params, funs}],header],Background -> NoneLightDarkSwitched[RGBColor[0.87, 0.94, 1]]None11 -> LightDarkSwitched[Hue[0.6, 0.4, 1]]12 -> LightDarkSwitched[Hue[0.6, 0.4, 1]], BaseStyle -> FontSize -> Larger, Frame -> All, FrameStyle -> StandardBlue, ItemStyle -> ScriptLevel -> 011 -> Bold12 -> Bold, Alignment -> CenterCenter, Spacings -> 31.5]//
(* Output *)
"Parameters" | ℱ_t[t-(1)/(2)](ω)
{0,1} | (ω sinc(ω)-ⅈ cos(ω)+ⅈ)/(Sqrt(2 π) ω)
{1,1} | (ω sinc(ω)-ⅈ cos(ω)+ⅈ)/(ω)
{-1,1} | (ω sinc(ω)-ⅈ cos(ω)+ⅈ)/(2 π ω)
{0,2 π} | -(ⅈ (-1+ℯ^(2 ⅈ π ω)))/(2 π ω)
```

#### GenerateConditions

Use [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html)->[True](https://reference.wolfram.com/language/ref/True.html) to get parameter conditions for when a result is valid:

```wolfram
FourierTransform[E^(-α t^2),t,ω,GenerateConditions->True]
(* Output *)
(ℯ^(-(ω^2)/(4 α)))/(Sqrt[2] Sqrt[α])
```

#### PrecisionGoal

The option [PrecisionGoal](https://reference.wolfram.com/language/ref/PrecisionGoal.html) sets the relative tolerance in the integration:

```wolfram
exact=FourierTransform[8/(t^2+1),t,1/2]
(* Output *)
4 Sqrt[(2 π)/(ℯ)]
```

```wolfram
FourierTransform[8/(t^2+1),t,.5, PrecisionGoal->15]-exact
(* Output *)
-0.008955638933840326+0. ⅈ
```

With default settings:

```wolfram
FourierTransform[8/(t^2+1),t,.5]-exact
(* Output *)
-7.219056641361021×10^-8+0. ⅈ
```

#### WorkingPrecision

If [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) is specified, the computation is done at that working precision:

```wolfram
exact=FourierTransform[8/(t^2+1),t,1/2]
(* Output *)
4 Sqrt[(2 π)/(ℯ)]
```

```wolfram
FourierTransform[8/(t^2+1),t,.5, WorkingPrecision->18]-exact
(* Output *)
-5.95001651314676726259145×10^-15
```

With default settings:

```wolfram
FourierTransform[8/(t^2+1),t,.5]-exact
(* Output *)
-7.219056641361021×10^-8+0. ⅈ
```

### Applications

#### Signals and Systems

Find the convolution of signals:

```wolfram
h[t_]:=UnitStep[t+1/2]-UnitStep[t-1/2];
x[t_]:=Cos[π t];
```

The product of their Fourier transforms:

```wolfram
FourierTransform[h[t],t,ω,FourierParameters->{1,-1}]*FourierTransform[x[t],t,ω,FourierParameters->{1,-1}]
(* Output *)
(2 (π DiracDelta[π-ω]+π DiracDelta[π+ω]) Sin[(ω)/(2)])/(ω)
```

Find the inverse transform:

```wolfram
InverseFourierTransform[%,ω,t,FourierParameters->{1,-1}]
(* Output *)
(2 Cos[π t])/(π)
```

Compare with [Convolve](https://reference.wolfram.com/language/ref/Convolve.html):

```wolfram
Convolve[h[τ],x[τ],τ,t]
(* Output *)
(2 Cos[π t])/(π)
```

Spectrum of the product of two signals, with one given in the frequency domain by:

```wolfram
Plot[UnitTriangle[ω],{ω,-1.2,1.2}]
```

*([Graphics])*

The Fourier transform of the signal $y(t)=cos(2 t)$:

```wolfram
FourierTransform[Cos[2 t],t,ω,FourierParameters->{1,-1}]
(* Output *)
π DiracDelta[-2+ω]+π DiracDelta[2+ω]
```

The Fourier transform of the product of $y(t)$ with the original signal is the convolution of its transforms:

```wolfram
1/2π Convolve[UnitTriangle[t],%/.ω->t,t,ω]
(* Output *)
(1)/(2) π (π UnitTriangle[2-ω]+π UnitTriangle[2+ω])
```

Its spectrum:

```wolfram
Plot[%,{ω,-5,5}]
```

*([Graphics])*

Frequency response of an LTI system defined by an ODE:

```wolfram
eqn=y'[t]+y[t]==x'[t];
```

Apply the Fourier transform over the equation:

```wolfram
FourierTransform[eqn,t,ω]
(* Output *)
FourierTransform[y[t],t,ω]-ⅈ ω FourierTransform[y[t],t,ω]==-ⅈ ω FourierTransform[x[t],t,ω]
```

Solve for the Fourier transform of $y(t)$:

```wolfram
Solve[%,FourierTransform[y[t],t,ω]]
(* Output *)
{{FourierTransform[y[t],t,ω]->(ω FourierTransform[x[t],t,ω])/(ⅈ+ω)}}
```

The frequency response of the LTI system is the ratio of the Fourier transforms of the output function $y(t)$ over the input function $x(t)$:

```wolfram
%[[1]][[1,2]]/FourierTransform[x[t],t,ω]
(* Output *)
(ω)/(ⅈ+ω)
```

#### Ordinary Differential Equations

Solve a differential equation using Fourier transforms:

```wolfram
eqn=y'[t]+y[t]==Sin[t];
```

Apply the Fourier transform over the equation:

```wolfram
FourierTransform[eqn,t,ω]
(* Output *)
FourierTransform[y[t],t,ω]-ⅈ ω FourierTransform[y[t],t,ω]==ⅈ Sqrt[(π)/(2)] DiracDelta[-1+ω]-ⅈ Sqrt[(π)/(2)] DiracDelta[1+ω]
```

Solve for the Fourier transform:

```wolfram
SolveValues[%,FourierTransform[y[t],t,ω]]
(* Output *)
{(-Sqrt[2 π] DiracDelta[-1+ω]+Sqrt[2 π] DiracDelta[1+ω])/(2 (ⅈ+ω))}
```

Find the inverse transform to get the solution:

```wolfram
InverseFourierTransform[%[[1]],ω,t]
(* Output *)
(1)/(2) (-Cos[t]+Sin[t])
```

Compare with [DSolveValue](https://reference.wolfram.com/language/ref/DSolveValue.html):

```wolfram
DSolveValue[{eqn,y[0]==-1/2},y[t],t]
(* Output *)
(1)/(2) (-Cos[t]+Sin[t])
```

#### Partial Differential Equations

Consider the heat equation: $\mathit{u}_{\mathit{t}}=\alpha^{2}\mathit{u}_{\mathit{x x}}$ with initial condition $\mathit{u}(\mathit{x},0)$:

```wolfram
eqn=D[u[x,t],t]==α^2 D[u[x,t],x,x];
```

Fourier transform with respect to $\mathit{x}$:

```wolfram
FourierTransform[eqn,x,ω]
(* Output *)
FourierTransform[u^(0,1)[x,t],x,ω]==-α^2 ω^2 FourierTransform[u[x,t],x,ω]
```

With $u1(\omega,t)=\mathcal{F}_{x}[u(x,t)](\omega)$ and $u1(\omega,0)=\omega_{0}$, solve this ODE:

```wolfram
DSolveValue[{D[u1[ω,t],t]==-α^2 ω^2 u1[ω,t],u1[ω,0]==ω_0},u1[ω,t],t]
(* Output *)
ℯ^(-t α^2 ω^2) ω_0
```

Compute the inverse Fourier transform:

```wolfram
InverseFourierTransform[1/Sqrt[2 π]E^(-ω^2α^2t),ω,x]
(* Output *)
(ℯ^(-(x^2)/(4 t α^2)))/(2 Sqrt[π] Sqrt[t α^2])
```

And convolution to get the solution:

```wolfram
Convolve[(ℯ^(-(x^2)/(4 t α^2)))/(2 Sqrt[π] Sqrt[t α^2])/.x->y,u[y,0],y,x]
(* Output *)
(Convolve[ℯ^(-(y^2)/(4 t α^2)),u[y,0],y,x])/(2 Sqrt[π] Sqrt[t α^2])
```

Consider the special case with initial condition $\mathit{u}(\mathit{x},0)=UnitBox[x]$ and $\alpha=1$:

```wolfram
%/.{u[y,0]->UnitBox[y],α->1}
(* Output *)
(1)/(2) (Erf[(1-2 x)/(4 Sqrt[t])]+Erf[(1+2 x)/(4 Sqrt[t])])
```

Compare with [DSolveValue](https://reference.wolfram.com/language/ref/DSolveValue.html):

```wolfram
DSolveValue[{eqn/.{α->1},{u[x,0]==UnitBox[x]}},u[x,t],{x,t}]
(* Output *)
(1)/(2) (-Erf[(-1+2 x)/(4 Sqrt[t])]+Erf[(1+2 x)/(4 Sqrt[t])])
```

Plot the initial conditions and solutions for different values of $t$:

```wolfram
Plot[{UnitBox[x],Evaluate@Table[%,{t,{.005,.1,.3,.8}}]},{x,-3,3},PlotLegends->{UnitBox[x],.005,.1,.3,.8},Exclusions->None]
(* Output *)
![image](img/image_003.png)
```

Plot the solution over the $x$-$t$ plane.

```wolfram
Plot3D[Evaluate[(1)/(2) (Erf[(1-2 x)/(4 Sqrt[t])]+Erf[(1+2 x)/(4 Sqrt[t])])], {x, -2, 2}, {t, 0, 1}, PlotRange -> All, PlotPoints -> 250, AxesLabel -> Automatic]
(* Output *)
![image](img/image_005.png)
```

#### Evaluation of Integrals

Calculate the following definite integral:

```wolfram
Inactive[Integrate][(u Cos[t u])/(1+u^2),{u,0,∞}]
(* Output *)
(u Cos[t u])/(1+u^2)
```

Compute the Fourier transform with respect to $t$ and interchange the order of transform and integration:

```wolfram
Inactive[Integrate][FourierTransform[(u Cos[t u])/(1+u^2),t,ω],{u,0,∞}]
(* Output *)
((Sqrt[(π)/(2)] u DiracDelta[-u+ω])/(1+u^2)+(Sqrt[(π)/(2)] u DiracDelta[u+ω])/(1+u^2))
```

Integrate over $u$:

```wolfram
Activate[%]
(* Output *)
(Sqrt[(π)/(2)] ω (-HeavisideTheta[-ω]+HeavisideTheta[ω]))/(1+ω^2)
```

Use the inverse Fourier transform to get the result:

```wolfram
FullSimplify[InverseFourierTransform[%,ω,t],Assumptions->t>0]
(* Output *)
-Cosh[t] CoshIntegral[t]+Sinh[t] SinhIntegral[t]
```

Compare with [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Integrate[(u Cos[t u])/(1+u^2),{u,0,∞},Assumptions->t>0]
(* Output *)
-Cosh[t] CoshIntegral[t]+Sinh[t] SinhIntegral[t]
```

#### Other Applications

The power spectrum of a damped sinusoid:

```wolfram
FourierTransform[Sin[10^3 t]Exp[-t/10]UnitStep[t],t,ω]
(* Output *)
(50000 Sqrt[(2)/(π)])/(100000001-20 ⅈ ω-100 ω^2)
```

```wolfram
LogLogPlot[Evaluate@Abs[%],{ω,1×10^-1,1×10^5},PlotRange->All]
```

*([Graphics])*

The Fourier transform of a radially symmetric function in the plane can be expressed as a Hankel transform. Verify this relation for the function defined by:

```wolfram
f[x_,y_]:=(x^2+y^2) E^(-Sqrt[x^2+y^2])
```

Plot the function:

```wolfram
Plot3D[f[x,y],{x,-3,3},{y,-3,3},PlotRange->All,Mesh->False]
```

*([Graphics3D])*

Compute its Fourier transform:

```wolfram
FourierTransform[f[x,y],{x,y},{u,v}]
(* Output *)
(6-9 (u^2+v^2))/((1+u^2+v^2)^(7/2))
```

Obtain the same result using [HankelTransform](https://reference.wolfram.com/language/ref/HankelTransform.html):

```wolfram
HankelTransform[f[x,y]/.{x->r Cos[m],y ->r Sin[m]}//Simplify,r,s]/.{s->Sqrt[u^2+v^2]}
(* Output *)
(6-9 (u^2+v^2))/((1+u^2+v^2)^(7/2))
```

Plot the Fourier transform:

```wolfram
Plot3D[%,{u,-3,3},{v,-3,3},PlotRange->All,Mesh->False]
```

*([Graphics3D])*

Generate a gallery of Fourier transforms for a list of radially symmetric functions:

```wolfram
flist={(1)/(Sqrt[r^2+1]),(Cos[2 π r])/(r),(1)/(r),ℯ^-r,BesselJ[0,r]^2,ℯ^(-r^2)};
```

Compute the Hankel transforms for these functions:

```wolfram
tlist=FullSimplify[HankelTransform[flist,r,s],Assumptions-> 0<s<2]
(* Output *)
{(ℯ^-s)/(s),0,(1)/(s),(1)/((1+s^2)^(3/2)),(2)/(π s Sqrt[4-s^2]),(1)/(2) ℯ^(-(s^2)/(4))}
```

Generate the gallery of Fourier transforms as required:

```wolfram
(Grid[#1,Alignment->{{Right,{Left}},Center},Spacings->0]&)[Table[Text[ReplaceAll[Part[flist, i], r -> Sqrt[x, ^, 2, +, y, ^, 2]]]Plot3D[Evaluate[N[ReplaceAll[Part[flist, i], r -> Sqrt[x, ^, 2, +, y, ^, 2]]]], x-22, y-22, PlotRange -> All, Exclusions -> None, Axes -> None, Ticks -> None, Boxed -> False]Style[[ShortRightArrow], 24, Bold]Plot3D[Evaluate[ReplaceAll[N[Part[tlist, i]], s -> Sqrt[u, ^, 2, +, v, ^, 2]]], u-22, v-22, PlotRange -> All, Exclusions -> None, Axes -> None, Ticks -> None, Boxed -> False]Text[ReplaceAll[Part[tlist, i], s -> Sqrt[u, ^, 2, +, v, ^, 2]]],{i,6}]]
(* Output *)
![image](img/image_007.png)
```

Calculate the power spectrum of a stationary [OrnsteinUhlenbeckProcess](https://reference.wolfram.com/language/ref/OrnsteinUhlenbeckProcess.html):

```wolfram
FourierTransform[CovarianceFunction[OrnsteinUhlenbeckProcess[μ,σ,θ],h],h,ω,FourierParameters->{1,1}, Assumptions->θ>0]
(* Output *)
(σ^2)/(θ^2+ω^2)
```

A quick look at the Heisenberg uncertainty principle:

Consider a fixed-area box function as the position space wavefunction of a particle. Its Fourier transform gives the momentum space wavefunction of the particle:

```wolfram
Simplify[FourierTransform[1/Sqrt[2α]UnitBox[x/(2α)],x,k,Assumptions->α∈PositiveIntegers],α∈PositiveIntegers]
(* Output *)
(Sqrt[α] Sinc[k α])/(Sqrt[π])
```

When $\alpha$ is small, the height of the fixed-area box is big, and the position of the particle is almost guaranteed. The momentum space wavefunction is approximately $\sqrt{\frac{\alpha}{\pi}}$ for values between its two roots closest to zero, which makes it almost impossible to find its momentum. Similarly, vice versa, as seen here:

```wolfram
Manipulate[GraphicsPlot[1/Sqrt[2α]UnitBox[x/(2α)],{x,-6.1,6.1},Exclusions -> None, AxesLabel -> Automatic, PlotLabel -> Position Space, PlotRange -> -0.11, Epilog -> RedArrowheads[-0.050.05]Arrow[-α0α0], ImageSize -> Small],Plot[Sqrt[, α], Sinc[k, , α], /, Sqrt[Pi], ), k-2Pi2Pi, PlotRange -> -0.41.5, AxesLabel -> Automatic, PlotLabel -> Momentum Space, Epilog -> RedArrowheads[-0.050.05]Arrow[-Pi/α)0Pi/α0], ImageSize -> Small],α16, ControlPlacement -> Bottom, SaveDefinitions -> True]
```

### Properties & Relations

By default, the Fourier transform of $f(t)$ is:

```wolfram
HoldForm[FourierTransform[f[t],t,ω]=HoldForm[1/Sqrt[2π]]*Integrate[f[t]ℯ^(ⅈ ω t),{t,-∞,∞}]]//
(* Output *)
ℱ_t[f(t)](ω)=(1)/(Sqrt(2 π)) ∫_-∞^∞f(t) ℯ^(ⅈ ω t)ⅆt
```

For $f(t)=e^{-t^{2}}cos(t)$, the definite integral becomes:

```wolfram
1/Sqrt[2π]Integrate[Exp[-t^2]Cos[t]Exp[ⅈ ω t],{t,-∞,∞}]//FullSimplify
(* Output *)
(ℯ^(-(1)/(4) (1+ω)^2) (1+ℯ^ω))/(2 Sqrt[2])
```

Compare with [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html):

```wolfram
FourierTransform[Exp[-t^2]Cos[t],t,ω]//TrigToExp//Factor//Simplify
(* Output *)
(ℯ^(-(1)/(4) (1+ω)^2) (1+ℯ^ω))/(2 Sqrt[2])
```

Use [Asymptotic](https://reference.wolfram.com/language/ref/Asymptotic.html) to compute an asymptotic approximation:

```wolfram
Asymptotic[Inactive[FourierTransform][E^(-t^2-t^4),t,x],x->0]
(* Output *)
(ℯ^(1/8) BesselK[(1)/(4),(1)/(8)])/(2 Sqrt[2 π])
```

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) and [InverseFourierTransform](https://reference.wolfram.com/language/ref/InverseFourierTransform.html) are mutual inverses:

```wolfram
InverseFourierTransform[FourierTransform[f[t],t,ω],ω,t]
(* Output *)
f[t]
```

```wolfram
FourierTransform[InverseFourierTransform[G[ ω],ω,t],t,ω]
(* Output *)
G[ω]
```

```wolfram
FourierTransform[1/(t^2+1),t,ω]
(* Output *)
ℯ^(-Abs[ω]) Sqrt[(π)/(2)]
```

```wolfram
InverseFourierTransform[%,ω,t]
(* Output *)
(1)/(1+t^2)
```

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) and [FourierCosTransform](https://reference.wolfram.com/language/ref/FourierCosTransform.html) are equal for even functions:

```wolfram
FourierTransform[Exp[-t^2],t,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

```wolfram
FourierCosTransform[Exp[-t^2],t,ω ]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) and [FourierSinTransform](https://reference.wolfram.com/language/ref/FourierSinTransform.html) differ by `ⅈ` for odd functions:

```wolfram
FourierSinTransform[t Exp[-Abs[t]],t,ω]
(* Output *)
(2 Sqrt[(2)/(π)] ω)/((1+ω^2)^2)
```

```wolfram
FourierTransform[t Exp[-Abs[t]],t,ω]
(* Output *)
(2 ⅈ Sqrt[(2)/(π)] ω)/((1+ω^2)^2)
```

### Possible Issues

The result from an inverse Fourier transform may not have the same form as the original:

```wolfram
FourierTransform[UnitStep[1+t]UnitStep[1-t],t,ω]
(* Output *)
(Sqrt[(2)/(π)] Sin[ω])/(ω)
```

```wolfram
InverseFourierTransform[%,ω,t]
(* Output *)
(1)/(2) (Sign[1-t]+Sign[1+t])
```

### Neat Examples

The Fourier transforms of weighted Hermite polynomials have a very simple form:

```wolfram
FourierTransform[Exp[-t^2]HermiteH[10,t],t,ω]
(* Output *)
-(ℯ^(-(ω^2)/(4)) ω^10)/(Sqrt[2])
```

```wolfram
InverseFourierTransform[%,ω,t]//Factor
(* Output *)
32 ℯ^(-t^2) (-945+9450 t^2-12600 t^4+5040 t^6-720 t^8+32 t^10)
```

Create a table of basic Fourier transforms:

```wolfram
flist={t^n,E^(a t),E^(-a t),Exp[-t^2],Sin[a t],t Sin[a t],Sinc[t],DiracDelta[t-a],Log[Abs[t]],UnitStep[t],UnitBox[t],BesselJ[2,a t],BesselY[0,Abs[t]]};
```

```wolfram
Grid[Prepend[{#,Assuming[{a>0},Simplify[FourierTransform[#1,t,ω]]]}&/@flist,{f[t],FourierTransform[f[t],t,ω]}],Background -> NoneLightDarkSwitched[RGBColor[0.87, 0.94, 1]]None11 -> LightDarkSwitched[Hue[0.6, 0.4, 1]]12 -> LightDarkSwitched[Hue[0.6, 0.4, 1]], BaseStyle -> FontSize -> Larger, Frame -> All, FrameStyle -> StandardBlue, ItemStyle -> ScriptLevel -> 011 -> Bold12 -> Bold, Alignment -> CenterCenter, Spacings -> 31.5]//
(* Output *)
f(t) | ℱ_t[f(t)](ω)
t^(n) | (ℯ^((ⅈ π n)/(2)) sin(π n) n+1 (ω-1) ω^(-n-1))/(Sqrt(2 π))
ℯ^(a t) | Sqrt(2 π) ω-ⅈ a
ℯ^(-a t) | Sqrt(2 π) ⅈ a+ω
ℯ^(-t^(2)) | (ℯ^(-(ω^(2))/(4)))/(Sqrt(2))
sin(a t) | ⅈ Sqrt((π)/(2)) ω-a-ⅈ Sqrt((π)/(2)) a+ω
t sin(a t) | Sqrt((π)/(2)) δ^(′)(ω-a)-Sqrt((π)/(2)) δ^(′)(a+ω)
sinc(t) | (1)/(2) Sqrt((π)/(2)) (1-ω+ω+1)
t-a | (cos(a ω)+ⅈ sin(a ω))/(Sqrt(2 π))
log(t) | -(Sqrt((π)/(2)))/(ω)-[DoubledGamma] Sqrt(2 π) ω
t | Sqrt((π)/(2)) ω+(ⅈ)/(Sqrt(2 π) ω)
t | (sinc((ω)/(2)))/(Sqrt(2 π))
2 |  | (Sqrt((2)/(π)) (a^(2)-2 ω^(2)))/(a^(2) Sqrt(a^(2)-ω^(2))) | (a^(2))/(ω^(2))>1
0 | True
0 |  | -(Sqrt((2)/(π)))/(Sqrt(ω^(2)-1)) | ω>1||ω<-1
0 | True
```

## Tech Notes ▪Integral Transforms and Related Operations

## Related Guides ▪Integral Transforms ▪Fourier Analysis ▪Signal Transforms ▪Generalized Functions ▪Calculus ▪Summation Transforms

## History Introduced in 1999 (4.0) | Updated in 2025 (14.2)
