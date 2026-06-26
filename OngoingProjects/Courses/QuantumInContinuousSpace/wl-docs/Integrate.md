# Integrate (∫)

> [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,*x*] — gives the indefinite integral $\int f d x$.
> [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,{*x*,*x*_*min*,*x*_*max*}] — gives the definite integral $\int_{x_{\mathit{min}}}^{x_{\mathit{max}}} f d x$.
> [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,{*x*,*x*_*min*,*x*_*max*},{*y*,*y*_*min*,*y*_*max*},…] — gives the multiple integral $\int_{x_{\mathit{min}}}^{x_{\mathit{max}}}d x \int_{y_{\mathit{min}}}^{y_{\mathit{max}}}d y \ldots f$.
> [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,{*x*,*y*,…}∈*reg*] — integrates over the geometric region `*reg*`.

## Details and Options

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,*x*] can be entered as `∫*f*ⅆ*x*`.

`∫` can be entered as esc int esc or [∫](https://reference.wolfram.com/language/ref/character/Integral.html).

`ⅆ` is not an ordinary `d`; it is entered as esc dd esc or [ⅆ](https://reference.wolfram.com/language/ref/character/DifferentialD.html).

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,{*x*,*y*,…}∈*reg*] can be entered as `*∫*_{*x*,*y*,…}∈*reg**f*`.

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,{*x*,*x*_*min*,*x*_*max*}] can be entered with `*x*_*min*` as a subscript and `*x*_*max*` as a superscript to `∫`.

Multiple integrals use a variant of the standard iterator notation. The first variable given corresponds to the outermost integral and is done last.

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) can evaluate integrals of rational functions. It can also evaluate integrals that involve exponential, logarithmic, trigonometric, and inverse trigonometric functions, so long as the result comes out in terms of the same set of functions.

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) can give results in terms of many special functions.

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) carries out some simplifications on integrals it cannot explicitly do.

You can get a numerical result by applying [N](https://reference.wolfram.com/language/ref/N.html) to a definite integral.

You can assign values to patterns involving [Integrate](https://reference.wolfram.com/language/ref/Integrate.html) to give results for new classes of integrals.

The integration variable can be a construct such as `*x*[*i*]` or any expression whose head is not a mathematical function.

For indefinite integrals, [Integrate](https://reference.wolfram.com/language/ref/Integrate.html) tries to find results that are correct for almost all values of parameters.

For definite integrals, the following options can be given:

| [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) | [$Assumptions](https://reference.wolfram.com/language/ref/$Assumptions.html) | assumptions to make about parameters |
| --- | --- | --- |
| [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to generate answers that involve conditions on parameters  |
| [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) | [None](https://reference.wolfram.com/language/ref/None.html) | how to name generated parameters |
| [PrincipalValue](https://reference.wolfram.com/language/ref/PrincipalValue.html) | [False](https://reference.wolfram.com/language/ref/False.html) | whether to find Cauchy principal values |

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) can evaluate essentially all indefinite integrals and most definite integrals listed in standard books of tables.

In [StandardForm](https://reference.wolfram.com/language/ref/StandardForm.html), [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*f*,*x*] is output as `∫*f*ⅆ*x*`.

## Examples

### Basic Examples

An indefinite integral:

```wolfram
Integrate[x^2+Sin[x],x]
(* Output *)
(x^3)/(3)-Cos[x]
```

Compute a definite integral:

```wolfram
Integrate[1/(x^3+1),{x,0,1}]
(* Output *)
(1)/(18) (2 Sqrt[3] π+Log[64])
```

Visualize the area given by this integral:

```wolfram
Plot[(1)/(x^3+1),{x,0,1},Filling->Axis]
```

*([Graphics])*

Use esc int esc to enter `∫` and esc dd esc to enter `ⅆ`:

```wolfram
∫Sqrt[x+Sqrt[x]]ⅆx
(* Output *)
(1)/(12) (Sqrt[Sqrt[x]+x] (-3+2 Sqrt[x]+8 x)+3 ArcTanh[(Sqrt[Sqrt[x]+x])/(1+Sqrt[x])])
```

```wolfram
[%]
(* Output *)
(1)/(12) (Sqrt(x+Sqrt(x)) (8 x+2 Sqrt(x)-3)+3 tanh^(-1)((Sqrt(x+Sqrt(x)))/(Sqrt(x)+1)))
```

Use ctrl to enter the lower limit, then ctrl for the upper limit:

```wolfram
∫_0^∞Log[x] Exp[-x^2]ⅆx
(* Output *)
-(1)/(4) Sqrt[π] (EulerGamma+Log[4])
```

### Scope

#### Basic Uses

Compute an indefinite integral:

```wolfram
int=Integrate[1/(x^3+1),x]
(* Output *)
(ArcTan[(-1+2 x)/(Sqrt[3])])/(Sqrt[3])+(1)/(3) Log[1+x]-(1)/(6) Log[1-x+x^2]
```

Verify the answer by differentiation:

```wolfram
D[int,x]//Simplify
(* Output *)
(1)/(1+x^3)
```

Use esc intt esc to enter a template $\int_\mathrm{d}_$ and tab to move between fields:

```wolfram
∫(a x^2 +b x +c)ⅆx
(* Output *)
c x+(b x^2)/(2)+(a x^3)/(3)
```

Include the constant of integration in an indefinite integral:

```wolfram
Integrate[x^2,x,GeneratedParameters->C]
(* Output *)
(x^3)/(3)+1
```

Compute a definite integral over a finite interval:

```wolfram
Integrate[(2x+3)/(x^2+5x+6),{x,0,2}]
(* Output *)
Log[(125)/(54)]
```

Infinite interval:

```wolfram
Integrate[BesselJ[2,x]/(1+x^2),{x,0,∞}]
(* Output *)
(1)/(6) (2-3 π BesselI[2,1]+3 π StruveL[2,1])
```

Doubly infinite interval:

```wolfram
Integrate[Sin[x]/x,{x,-∞,∞}]
(* Output *)
π
```

Use esc dintt esc to enter a template $\int_{_}^{_}_\mathrm{d}_$ and tab to move between fields:

```wolfram
∫_0^πCos[x]ⅆx
(* Output *)
0
```

Integrate a function with a symbolic parameter:

```wolfram
Integrate[(a x)/(x^3+2),{x,1,∞}]
(* Output *)
(a (Sqrt[3] π-2 Sqrt[3] ArcTan[(-1+2^(2/3))/(Sqrt[3])]-Log[4+2 2^(1/3)-2 2^(2/3)]+2 Log[2+2^(2/3)]))/(6 2^(1/3))
```

An integral that only converges for some values of parameters:

```wolfram
Integrate[Exp[-c x^2],{x,-∞,∞}]
(* Output *)
(Sqrt[π])/(Sqrt[c])
```

Specify alternate assumptions to use:

```wolfram
Integrate[Exp[-c x^2],{x,-∞,∞},Assumptions->Re[c]==0]
(* Output *)
((-c)^(3/4) Sqrt[(π)/(2)] (-1+Sign[c]))/(c^(5/4))
```

Multivariate integrals:

```wolfram
Integrate[(BesselJ[3,y])/(x+1),{x,0,5},{y,0,5}]
(* Output *)
-((-1+(2)/(5) BesselJ[1,5]+BesselJ[2,5]) Log[6])
```

```wolfram
Integrate[Exp[-x+y],{x,0,∞},{y,0,1}]
(* Output *)
-1+ℯ
```

Multiple integral with `*x*` integration last:

```wolfram
Integrate[Sin[x y],{x,0,1},{y,0,x}]
(* Output *)
(1)/(2) (EulerGamma-CosIntegral[1])
```

In [StandardForm](https://reference.wolfram.com/language/ref/StandardForm.html), the differential `ⅆ*y*` precedes `ⅆ*x*`:

```wolfram
∫_0^1∫_0^xSin[x y]ⅆyⅆx
(* Output *)
(1)/(2) (EulerGamma-CosIntegral[1])
```

Visualize the function over the domain of integration:

```wolfram
Plot3D[Sin[x y],{x,y}∈Triangle[{{0,0},{0,1},{1,0}}]]
```

*([Graphics3D])*

Integrals over standard regions:

```wolfram
Integrate[1,{x,y}∈Circle[]]
(* Output *)
2 π
```

The character `∈` can be entered as esc el esc or [∈](https://reference.wolfram.com/language/ref/character/Element.html):

```wolfram
Integrate[1,x∈Ball[3]]
(* Output *)
(4 π)/(3)
```

Enter a region specification $vars \in reg$ in an underscript using ctrl:

```wolfram
∫_{{x,y,z}∈Tetrahedron[]}(x^4 y^5+z^10)
(* Output *)
(1385)/(759103488 Sqrt[2])
```

Use esc rintt esc to enter a template $\int_{_\in_}_$ and tab to move between fields:

```wolfram
∫_{{x,y}∈Rectangle[]}x y^2
(* Output *)
(1)/(6)
```

Formal integrals:

```wolfram
Integrate[f'[x],x]
(* Output *)
f[x]
```

```wolfram
Integrate[((1-f[x])f'[x])/E^f[x], x]
(* Output *)
ℯ^(-f[x]) f[x]
```

Integrals of vector- and array-valued functions:

```wolfram
Integrate[{x,(1)/(Sqrt[x]),Sin[x]},{x,0,5}]
(* Output *)
{(25)/(2),2 Sqrt[5],1-Cos[5]}
```

```wolfram
Integrate[({{x, Sqrt[x], x^(1/3)}, {Sqrt[x], x^(1/3), x^(1/4)}, {x^(1/3), x^(1/4), x^(1/5)}}),{x,0,5}]//MatrixForm
(* Output *)
({{(25)/(2), (10 Sqrt[5])/(3), (15 5^(1/3))/(4)}, {(10 Sqrt[5])/(3), (15 5^(1/3))/(4), 4 5^(1/4)}, {(15 5^(1/3))/(4), 4 5^(1/4), (25 5^(1/5))/(6)}})
```

Invoke [NIntegrate](https://reference.wolfram.com/language/ref/NIntegrate.html) automatically if symbolic integration fails:

```wolfram
Integrate[E^(-x^x),{x,0,1}]
(* Output *)
∫_0^1ℯ^(-x^x)ⅆx
```

```wolfram
N[%]
(* Output *)
0.45854198647406746
```

#### Indefinite Integrals

Some basic integrals:

```wolfram
∫x^nⅆx
(* Output *)
(x^(1+n))/(1+n)
```

```wolfram
∫(1)/(x)ⅆx
(* Output *)
Log[x]
```

```wolfram
∫Sin[x]ⅆx
(* Output *)
-Cos[x]
```

```wolfram
∫Exp[x]ⅆx
(* Output *)
ℯ^x
```

Generate an answer with a constant of integration:

```wolfram
Integrate[Exp[x], x, GeneratedParameters->C]
(* Output *)
ℯ^x+1
```

Integrals of trigonometric functions:

```wolfram
∫Sec[x]^2ⅆx
(* Output *)
Tan[x]
```

```wolfram
∫Cot[x]ⅆx
(* Output *)
Log[Sin[x]]
```

```wolfram
∫Sec[x]ⅆx
(* Output *)
ArcCoth[Sin[x]]
```

Verify the previous answer via differentiation:

```wolfram
D[%,x]//Simplify
(* Output *)
Sec[x]
```

Create a nicely formatted table of integrals:

```wolfram
flist={x^n,1/x,b^x,Log[x],Sin[x],Cos[x],(1)/(x^2+a^2),(1)/(x^2-a^2)};
```

```wolfram
Grid[Prepend[Transpose[{flist,Integrate[flist,x]}],{f[x],∫f[x]ⅆx}],Background -> NoneLightDarkSwitched[RGBColor[0.87, 0.94, 1]]None11 -> LightDarkSwitched[Hue[0.6, 0.4, 1]]12 -> LightDarkSwitched[Hue[0.6, 0.4, 1]], BaseStyle -> FontSize -> Larger, Frame -> All, FrameStyle -> StandardBlue, ItemStyle -> ScriptLevel -> 011 -> Bold12 -> Bold, Alignment -> CenterCenter, Spacings -> 31.5]//
(* Output *)
f(x) | ∫f(x)ⅆx
x^(n) | (x^(n+1))/(n+1)
(1)/(x) | log(x)
b^(x) | (b^(x))/(log(b))
log(x) | x log(x)-x
sin(x) | -cos(x)
cos(x) | sin(x)
(1)/(a^(2)+x^(2)) | (tan^(-1)((x)/(a)))/(a)
(1)/(x^(2)-a^(2)) | -(tanh^(-1)((x)/(a)))/(a)
```

Rational functions can always be integrated in closed form:

```wolfram
Integrate[1/(x^4-1),x]
(* Output *)
-(ArcTan[x])/(2)+(1)/(4) Log[1-x]-(1)/(4) Log[1+x]
```

```wolfram
Integrate[( x + 1)/(a^(2 ) x^2 +c^2),x]
(* Output *)
(ArcTan[(a x)/(c)])/(a c)+(Log[c^2+a^2 x^2])/(2 a^2)
```

Sometimes they involve sums of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
Integrate[1/(x^5+2x+1),x]
(* Output *)
RootSum[1+2 #1+#1^5&,(Log[x-#1])/(2+5 #1^4)&]
```

Integrals of general elementary functions:

```wolfram
∫ArcSin[x]ⅆx
(* Output *)
Sqrt[1-x^2]+x ArcSin[x]
```

```wolfram
∫Cosh[x]ⅆx
(* Output *)
Sinh[x]
```

```wolfram
∫(x^2)/(Sqrt[a^2+x^2])ⅆx
(* Output *)
(1)/(2) x Sqrt[a^2+x^2]-(1)/(2) a^2 ArcTanh[(x)/(Sqrt[a^2+x^2])]
```

```wolfram
∫Exp[a x]Cos[b x]ⅆx
(* Output *)
(ℯ^(a x) (a Cos[b x]+b Sin[b x]))/(a^2+b^2)
```

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) returns antiderivatives valid in the complex plane where applicable:

```wolfram
intℂ[x_]=∫Tan[x]ⅆx
(* Output *)
-Log[Cos[x]]
```

```wolfram
D[intℂ[x],x]
(* Output *)
Tan[x]
```

A common antiderivative found in integral tables for $tan(x)$ is $log(sec(x))$:

```wolfram
intℝ[x_]=Log[RealAbs[Sec[x]]];
```

This is a valid antiderivative for real values of $x$:

```wolfram
Simplify[D[intℝ[x],x]]
(* Output *)
Tan[x]
```

On the real line, the two integrals have the same real part:

```wolfram
FullSimplify[Re[intℂ[x]]==Re[intℝ[x]]]
(* Output *)
True
```

```wolfram
Plot[{Re[intℂ[x]],Re[intℝ[x]]},{x,0,4Pi},PlotTheme->{"Detailed","NeonColors","DashedLines"}]
(* Output *)
![image](img/image_001.png)
```

But the imaginary parts differ by $\pi$ on any interval where $cos(x)$ is negative:

```wolfram
Plot[{Im[intℂ[x]],Im[intℝ[x]],Cos[x]},{x,0,4Pi},PlotTheme->{"Detailed","NeonColors","DashedLines"}]
(* Output *)
![image](img/image_003.png)
```

Similar integrals can lead to functions of different kinds:

```wolfram
Integrate[Sqrt[x]Sqrt[1+x],x]
(* Output *)
(1)/(4) (Sqrt[x] Sqrt[1+x] (1+2 x)-ArcTanh[(1)/(Sqrt[(x)/(1+x)])])
```

```wolfram
Integrate[Sqrt[x]Sqrt[1+x]Sqrt[2+x],x]
(* Output *)
(2)/(5) (Sqrt[x] (1+x)^(3/2) Sqrt[2+x]-(2 Sqrt[(1+x)/(2+x)] (2+x))/(Sqrt[x])-(2 ⅈ Sqrt[(1+x)/(2+x)] Sqrt[(2+x)/(x)] EllipticE[ⅈ ArcSinh[(1)/(Sqrt[x])],2])/(Sqrt[1+(1)/(x)]))
```

Many integrals can be done only in terms of special functions such as [Erf](https://reference.wolfram.com/language/ref/Erf.html):

```wolfram
Integrate[Exp[-x^2],x]
(* Output *)
(1)/(2) Sqrt[π] Erf[x]
```

Generalizations of [Log](https://reference.wolfram.com/language/ref/Log.html) such as [PolyLog](https://reference.wolfram.com/language/ref/PolyLog.html) and [LogIntegral](https://reference.wolfram.com/language/ref/LogIntegral.html):

```wolfram
Integrate[Log[1+x]^2/x,x]
(* Output *)
Log[-x] Log[1+x]^2+2 Log[1+x] PolyLog[2,1+x]-2 PolyLog[3,1+x]
```

```wolfram
Integrate[Log[Log[x]],x]
(* Output *)
x Log[Log[x]]-LogIntegral[x]
```

Hypergeometric functions such as [Hypergeometric2F1](https://reference.wolfram.com/language/ref/Hypergeometric2F1.html):

```wolfram
Integrate[Tan[x]^n,x]
(* Output *)
(Hypergeometric2F1[1,(1+n)/(2),(3+n)/(2),-Tan[x]^2] Tan[x]^(1+n))/(1+n)
```

Create a nicely-formatted table of special function integrals:

```wolfram
flist={BesselJ[0,x],Erf[x],EllipticE[x],Sqrt[x]PolyLog[2,x],AiryAi[x],SinhIntegral[2/(x+1)]};
```

```wolfram
Grid[Prepend[Transpose[{flist,Integrate[flist,x]}],{f[x],∫f[x]ⅆx}],Background -> NoneLightDarkSwitched[RGBColor[0.87, 0.94, 1]]None11 -> LightDarkSwitched[Hue[0.6, 0.4, 1]]12 -> LightDarkSwitched[Hue[0.6, 0.4, 1]], BaseStyle -> FontSize -> Larger, Frame -> All, FrameStyle -> StandardBlue, ItemStyle -> ScriptLevel -> 011 -> Bold12 -> Bold, Alignment -> CenterCenter, Spacings -> 31.5]//
(* Output *)
f(x) | ∫f(x)ⅆx
0 | x _1F_2((1)/(2);1,(3)/(2);-(x^(2))/(4))
erf(x) | x erf(x)+(ℯ^(-x^(2)))/(Sqrt(π))
x | (2 x x)/(3)-(2 x)/(3)+(2 x x)/(3)+(2 x)/(3)
Sqrt(x) 2 | (2)/(3) x^(3/2) 2+(4)/(27) (Sqrt(x) (3 x log(1-x)-2 (x+3))+6 tanh^(-1)(Sqrt(x)))
x | -(x ((3)^(1/3) x (2)/(3)^(2) _1F_2((2)/(3);(4)/(3),(5)/(3);(x^(3))/(9))-3 (1)/(3) (5)/(3) _1F_2((1)/(3);(2)/(3),(4)/(3);(x^(3))/(9))))/(9 3^(2/3) (2)/(3) (4)/(3) (5)/(3))
(2)/(x+1) | (x+1) ((2)/(x+1)+sinh((2)/(x+1)))-2 (2)/(x+1)
```

The variable of integration need not be a single symbol:

```wolfram
Integrate[p[x]Log[p[x]],p[x]]
(* Output *)
-(1)/(4) p[x]^2+(1)/(2) Log[p[x]] p[x]^2
```

#### Definite Integrals

Integrate a polynomial:

```wolfram
Integrate[x^4+ x^2+1,{x,1,3}]
(* Output *)
(886)/(15)
```

Visualize the area under the curve:

```wolfram
Plot[x^4+ x^2+1,{x,1,3},Filling->Axis]
```

*([Graphics])*

Integrate a symbolic polynomial:

```wolfram
Integrate[a x^2 + b x+c,{x,-2,2}]
(* Output *)
(16 a)/(3)+4 c
```

Integrate over a symbolic range:

```wolfram
Integrate[x^n + n x,{x,0,a},Assumptions->n>1&&n∈Integers]
(* Output *)
(1)/(2) a (a n+(2 a^n)/(1+n))
```

Rational functions:

```wolfram
Integrate[1/(x^4+x^2+1),{x,0,Infinity}]
(* Output *)
(π)/(2 Sqrt[3])
```

```wolfram
Plot[(1)/(x^4+x^2+1),{x,0,5},Filling->Axis,PlotRange->All,Epilog->{[Graphics],Arrow[{{0,0},{5,0}}]}]
```

*([Graphics])*

```wolfram
Integrate[(4x^2-7x-12)/((x+2)(x-3)),{x,-1,2}]
(* Output *)
12-(21 Log[4])/(5)
```

```wolfram
Plot[(4 x^2-7 x-12)/((x+2) (x-3)),{x,-1,2},Filling->0,PlotRange->All]
```

*([Graphics])*

Algebraic functions:

```wolfram
Integrate[(1+x^3) CubeRoot[x],{x,-1,1}]
(* Output *)
(6)/(13)
```

```wolfram
Plot[(1+x^3) CubeRoot[x],{x,-1,1},Filling->Axis,FillingStyle->{<|color -> RGBColor[0.95, 0.627, 0.1425, 0.2]|>,<|color -> RGBColor[0.24, 0.6, 0.8, 0.2]|>}]
```

*([Graphics])*

```wolfram
Integrate[1/((2+x^2)Sqrt[4+3x^2]),{x,-∞,∞}]
(* Output *)
ArcCosh[Sqrt[(3)/(2)]]
```

```wolfram
Plot[(1)/((2+x^2) Sqrt[4+3 x^2]),{x,-7,7},Filling->0,AxesOrigin->{0,0},Epilog->{[Graphics],Arrowheads[{-Medium,Medium}],Arrow[{{-7,0},{7,0}}]}]
```

*([Graphics])*

Trigonometric functions:

```wolfram
Integrate[Cos[2x]^4,{x,0,π}]
(* Output *)
(3 π)/(8)
```

```wolfram
Plot[Cos[2x]^4,{x,0,π/2},Filling->Axis]
```

*([Graphics])*

```wolfram
Integrate[x^2Sin[2x],{x,0,2π}]
(* Output *)
-2 π^2
```

```wolfram
Plot[x^2Sin[2x],{x,0,2π},Filling->Axis,FillingStyle->{<|color -> RGBColor[0.95, 0.627, 0.1425, 0.2]|>,<|color -> RGBColor[0.24, 0.6, 0.8, 0.2]|>}]
```

*([Graphics])*

Exponential and logarithmic functions:

```wolfram
Integrate[x^5Log[x],{x,1,2}]
(* Output *)
-(7)/(4)+(16 Log[64])/(9)
```

```wolfram
Plot[x^5Log[x],{x,1,2},Filling->Axis]
```

*([Graphics])*

```wolfram
Integrate[Exp[-x+Exp[-x]],{x,0,∞}]
(* Output *)
-1+ℯ
```

```wolfram
Plot[Exp[-x+Exp[-x]],{x,0,4},PlotRange->Full,Filling->Axis,Epilog->{[Graphics],Arrowheads[Medium],Arrow[{{0,0},{4,0}}]}]
```

*([Graphics])*

Hyperbolic trigonometric functions:

```wolfram
Integrate[x Cosh[x],{x,0,a}]
(* Output *)
1-Cosh[a]+a Sinh[a]
```

```wolfram
Integrate[x^2 Sech[x]^2,{x,-∞,∞}]
(* Output *)
(π^2)/(6)
```

Integrate a function with a vertical asymptote:

```wolfram
Integrate[x/Sqrt[1-x],{x,0,1}]
(* Output *)
(4)/(3)
```

```wolfram
Plot[x/Sqrt[1-x],{x,0,1},Filling->Axis]
```

*([Graphics])*

This can be viewed as a limit of the result of integration on a smaller interval:

```wolfram
[Limit]_{a->1^(-)}∫_0^a(x)/(Sqrt[1-x])ⅆx
(* Output *)
(4)/(3)
```

Compute the integral of a function with two vertical asymptotes:

```wolfram
area=∫_-1^1(x^2)/(Sqrt[1-x^4])ⅆx
(* Output *)
(2 Sqrt[π] Gamma[(7)/(4)])/(3 Gamma[(5)/(4)])
```

```wolfram
%//N
(* Output *)
1.198140234735592
```

```wolfram
Plot[(x^2)/(Sqrt[1-x^4]),{x,-1,1},Filling->0]
```

*([Graphics])*

This can be viewed as a multivariate limit of the result of integration on a smaller interval:

```wolfram
area==Limit[Integrate[(x^2)/(Sqrt[1-x^4]),{x,a,b},Assumptions->-1<a<b<1],{a,b}->{-1,1},Direction->{"FromAbove","FromBelow"}]
(* Output *)
True
```

Integrals over infinite intervals can be viewed as limits of integrals over finite domains:

```wolfram
Integrate[x E^(-x),{x,1,∞}]
(* Output *)
(2)/(ℯ)
```

The preceding is the limit as $a \to \infty$ of the integral from $1$ to $a$:

```wolfram
% == [Limit]_{a->∞}Integrate[x E^(-x),{x,1,a}]
(* Output *)
True
```

An integral over reals:

```wolfram
∫_-∞^∞ℯ^(-x^2)ⅆx
(* Output *)
Sqrt[π]
```

It is the bivariate limit of a finite integral:

```wolfram
%==[Limit]_{{a,b}->{-∞,∞}}∫_a^bℯ^(-x^2)ⅆx
(* Output *)
True
```

When there are parameters, conditions that ensure convergence may be reported:

```wolfram
Integrate[x^n,{x,0,1}]
(* Output *)
(1)/(1+n)
```

```wolfram
Integrate[E^(a x),{x,0,∞}]
(* Output *)
-(1)/(a)
```

Integrals of elementary functions may produce special function answers:

```wolfram
Integrate[Sin[t]/t,{t,a,b}]
(* Output *)
-SinIntegral[a]+SinIntegral[b]
```

```wolfram
Integrate[Cos[Sin[x]^2],{x,0,2π}]
(* Output *)
2 π BesselJ[0,(1)/(2)] Cos[(1)/(2)]
```

```wolfram
%//N
(* Output *)
5.174735523103564
```

```wolfram
Plot[Cos[Sin[x]^2],{x,0,2π},Filling->0,PlotRange->{0,1}]
```

*([Graphics])*

Create a formatted table of definite integrals over the positive reals of special functions:

```wolfram
flist={BesselJ[2,x]/x,BesselK[0,x]^2,AiryAi[x]^2,Exp[-x] Sinc[x],Sin[x]Erfc[x]};
```

```wolfram
Grid[Prepend[Transpose[{flist,∫_0^∞flistⅆx}],{f[x],∫_0^∞f[x]ⅆx}],Background -> NoneLightDarkSwitched[RGBColor[0.87, 0.94, 1]]None11 -> LightDarkSwitched[Hue[0.6, 0.4, 1]]12 -> LightDarkSwitched[Hue[0.6, 0.4, 1]], BaseStyle -> FontSize -> Larger, Frame -> All, FrameStyle -> StandardBlue, ItemStyle -> ScriptLevel -> 011 -> Bold12 -> Bold, Alignment -> CenterCenter, Spacings -> 31.5]//
(* Output *)
f(x) | ∫_0^∞f(x)ⅆx
(2)/(x) | (1)/(2)
0^(2) | (π^(2))/(4)
x^(2) | (1)/(3^(2/3) (1)/(3)^(2))
ℯ^(-x) sinc(x) | (π)/(4)
erfc(x) sin(x) | 1-(1)/((ℯ)^(1/4))
```

Integral along a complex line:

```wolfram
Integrate[Sqrt[x],{x,I,3-I}]
(* Output *)
-(2)/(3) ((-1)^(3/4)-(3-ⅈ)^(3/2))
```

Integrate along a piecewise linear contour in the complex plane:

```wolfram
Integrate[(1)/(z+1/2),{z,1,ℯ^((ⅈ π)/(3)),ℯ^((2 ⅈ π)/(3)),-1,ℯ^(-(2 ⅈ π)/(3)),ℯ^(-(ⅈ π)/(3)),1}]//FullSimplify
(* Output *)
2 ⅈ π
```

Integrate the same function along a circular contour:

```wolfram
Integrate[(ⅈ Exp[ⅈ ω])/(Exp[ⅈ ω]+1/2),{ω,0,2π}]
(* Output *)
2 ⅈ π
```

Plot the function and paths of integration:

```wolfram
ComplexPlot[(1)/(z+1/2),{z,-1.2-1.2I,1.2+1.2I},Epilog->{Thick,Arrowheads[Medium],White,Arrow[{{1,-.05},{1,0}}],Circle[],Dashed,Black,Arrow[ReIm/@Exp[2π ⅈ/6 Range[0,6]]]}]
```

![image](img/image_005.png)

#### Integrals of Piecewise and Generalized Functions

Compute the indefinite integral of a [Piecewise](https://reference.wolfram.com/language/ref/Piecewise.html) function:

```wolfram
f[x_]  := {, {{-x^2, x <0}, {x^2, True}}}
```

```wolfram
∫f[x]ⅆx
(* Output *)
{, {{-(x^3)/(3), x<=0}, {(x^3)/(3), True}}}
```

In this case, the derivative of the integral equals the original function:

```wolfram
Simplify[D[%,x]==f[x]]
(* Output *)
True
```

Integrate a discontinuous [Piecewise](https://reference.wolfram.com/language/ref/Piecewise.html) function:

```wolfram
f[x_]  := {, {{-x+1, x <0}, {x, 0<x<1}, {x^2, True}}}
```

```wolfram
g[x_]=∫f[x]ⅆx
(* Output *)
{, {{x-(x^2)/(2), x<=0}, {(x^2)/(2), 0<x<=1}, {(1)/(6)+(x^3)/(3), True}}}
```

Except at the point of discontinuity, the derivative of `g` equals `f`:

```wolfram
Resolve[ForAll[x,x!=0&&x∈Reals,g^′[x]==f[x]]]
(* Output *)
True
```

Visualize the function and its antiderivative:

```wolfram
Plot[{f[x],g[x]},{x,-2,2},PlotTheme->{"Detailed","DashedLines"}]
```

*([Graphics])*

Integrate functions that are piecewise-defined:

```wolfram
unitInt[x_]=Integrate[x UnitStep[x], x]
(* Output *)
(1)/(2) x^2 UnitStep[x]
```

```wolfram
Plot[{x UnitStep[x],unitInt[x]}, {x,-2,2},PlotTheme->{"Detailed","DashedLines"}]
```

*([Graphics])*

```wolfram
clipInt[x_]=Integrate[Clip[x]^2, x]
(* Output *)
{, {{x, x<=-1}, {-(2)/(3)+(x^3)/(3), -1<x<=1}, {-(4)/(3)+x, True}}}
```

```wolfram
Plot[{Clip[x]^2,clipInt[x]}, {x,-2,2},PlotTheme->{"Detailed","DashedLines"}]
(* Output *)
![image](img/image_006.png)
```

Integrate a piecewise function with infinitely many cases:

```wolfram
maxInt[x_]=Integrate[Max[Sin[x],Cos[x]],x]
(* Output *)
{, {{Sin[x], Cos[x]-Sin[x]>=0}, {-Cos[x], True}}}
```

Everywhere it is defined, the derivative of `maxInt` equals the original function:

```wolfram
Simplify[D[maxInt[x],x]==Max[Sin[x],Cos[x]]]
(* Output *)
True
```

However, `maxInt` itself is discontinuous:

```wolfram
Plot[{Max[Sin[x],Cos[x]],maxInt[x]},{x,0,4Pi},PlotTheme->{"Detailed","DashedLines"}]
(* Output *)
![image](img/image_008.png)
```

Compute a definite integral of a [Piecewise](https://reference.wolfram.com/language/ref/Piecewise.html) function:

```wolfram
f[x_]:={, {{1, x<2}, {x^3, 2<x<3}, {0, True}}}
```

```wolfram
Integrate[f[x],{x,0,3}]
(* Output *)
(73)/(4)
```

Compute the integral with a variable endpoint:

```wolfram
int[x_]=Integrate[f[t],{t,0,x},Assumptions->x>0]
(* Output *)
{, {{(73)/(4), x>3}, {x, x<=2}, {(1)/(4) (-8+x^4), True}}}
```

Visualize the function and its integral:

```wolfram
Plot[{f[x],int[x]},{x,0,4},PlotTheme->{"Detailed","DashedLines"}]
```

*([Graphics])*

Compute definite integrals of piecewise functions such as [Floor](https://reference.wolfram.com/language/ref/Floor.html):

```wolfram
Integrate[Floor[x^2],{x,0,3}]
(* Output *)
21-3 Sqrt[2]-Sqrt[3]-Sqrt[5]-Sqrt[6]-Sqrt[7]
```

```wolfram
Plot[Floor[x^2],{x,0,3},Filling->Axis]
```

*([Graphics])*

[PrimePi](https://reference.wolfram.com/language/ref/PrimePi.html):

```wolfram
Integrate[PrimePi[n],{n,0,100}]
(* Output *)
1440
```

```wolfram
Plot[PrimePi[x],{x,0,100},Filling->Axis]
```

*([Graphics])*

A composition of piecewise functions:

```wolfram
Integrate[Ceiling[x^2+TriangleWave[x]],{x,-2,2} ]
(* Output *)
42-3 Sqrt[2]-3 Sqrt[3]-Sqrt[5]-Sqrt[6]-Sqrt[7]-Sqrt[10]-Sqrt[11]-Sqrt[13]-Sqrt[14]-Sqrt[15]
```

Compute the definite integral with a variable upper limit:

```wolfram
f[x_]:=x^2 UnitStep[x+1]
```

```wolfram
int[a_]=Integrate[f[x],{x,-3,a}]
(* Output *)
(1)/(3) (1+a^3) UnitStep[1+a]
```

```wolfram
Plot[{f[a],int[a]},{a,-2,2},PlotTheme->{"Detailed","DashedLines"}]
(* Output *)
![image](img/image_010.png)
```

A function with an infinite number of cases:

```wolfram
f[x_]:=Max[Sin[x],Cos[x]]
```

```wolfram
f[x]//PiecewiseExpand
(* Output *)
{, {{Cos[x], Cos[x]-Sin[x]>=0}, {Sin[x], True}}}
```

Integrate over a finite number of cases using [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html):

```wolfram
maxInt[a_]=Integrate[f[x],{x,0,a},Assumptions->0<a<2Pi]
(* Output *)
{, {{(3)/(Sqrt[2]), a==2 π-2 ArcTan[1+Sqrt[2]]}, {1+Sqrt[2], a==π}, {Sqrt[2]-Cos[a], π<a<2 π-2 ArcTan[1+Sqrt[2]]||-2 ArcTan[1-Sqrt[2]]<a<π}, {Sin[a], a<=-2 ArcTan[1-Sqrt[2]]}, {2 Sqrt[2]+Sin[a], True}}}
```

The integral is a continuous function of the upper limit over the domain of integration:

```wolfram
Plot[{f[a],maxInt[a]},{a,0,2Pi},PlotTheme->{"Detailed","DashedLines"}]
(* Output *)
![image](img/image_012.png)
```

Integrate generalized functions:

```wolfram
Integrate[DiracDelta[x],{x,-∞,∞}]
(* Output *)
1
```

```wolfram
Integrate[Sin[x]^2 HeavisideTheta[x+π]HeavisideTheta[π-x],{x,-∞,∞}]
(* Output *)
π
```

```wolfram
Integrate[DiracDelta[Cos[x^2-1]],{x,-2,2}]
(* Output *)
Sqrt[(2)/(2+π)]
```

Indefinite integrals of generalized functions return generalized functions:

```wolfram
Integrate[DiracDelta[x],x]
(* Output *)
HeavisideTheta[x]
```

A nested integral:

```wolfram
Integrate[DiracDelta[x],x,x]
(* Output *)
x HeavisideTheta[x]
```

Integrate generalized functions over subsets of the reals:

```wolfram
Integrate[f[x] DiracDelta[x-a],{x,-1,1},Assumptions->a∈Reals]
(* Output *)
f[a] HeavisideTheta[1-a] HeavisideTheta[1+a]
```

```wolfram
Integrate[Cos[x-a] HeavisideLambda[x],{x,0,Pi},Assumptions->a∈Reals]
(* Output *)
-Cos[1-a]+Cos[a]+Sin[a]
```

Integrate an interpolating function:

```wolfram
f=Interpolation[{10,-5,10,-10,20}]
(* Output *)
InterpolatingFunction[...]
```

```wolfram
g=Head[Integrate[f[x],x]]
(* Output *)
InterpolatingFunction[...]
```

Test that `g` is a correct antiderivative at `*x*==3.5:

```wolfram
g'[3.5]==f[3.5]
(* Output *)
True
```

Visualize the function and its antiderivative:

```wolfram
Plot[{f[x],g[x]},{x,1,5},PlotTheme->{"Detailed","DashedLines"}]
(* Output *)
![image](img/image_014.png)
```

#### Nested Integrals

Compute a second antiderivative of a function:

```wolfram
Integrate[a x^2 + b x +c, x, x]//Expand
(* Output *)
(c x^2)/(2)+(b x^3)/(6)+(a x^4)/(12)
```

```wolfram
D[%,{x,2}]
(* Output *)
c+b x+a x^2
```

Compute the third antiderivative:

```wolfram
∫∫∫(1)/(x^2+4)ⅆxⅆxⅆx
(* Output *)
(1)/(4) ((-4+x^2) ArcTan[(x)/(2)]-2 x (-1+Log[4+x^2]))
```

```wolfram
FullSimplify[D[%,{x,3}]]
(* Output *)
(1)/(4+x^2)
```

Integrate a function with respect to two different variables:

```wolfram
Integrate[x^2Sin[y], y, x]
(* Output *)
-(1)/(3) x^3 Cos[y]
```

The mixed partial derivative gives the original function:

```wolfram
D[%,x,y]
(* Output *)
x^2 Sin[y]
```

Generate a constant of integration for a single integral:

```wolfram
Integrate[x,x,GeneratedParameters->C]
(* Output *)
(x^2)/(2)+1
```

Generate constants for a nested integral with respect to the same variable:

```wolfram
Integrate[x,x,x,GeneratedParameters->C]
(* Output *)
(x^3)/(6)+1+x 2
```

This is the most general second antiderivative of the integrand:

```wolfram
D[%,x,x]
(* Output *)
x
```

Generate two functions of integration for a nested integral with respect to two variables:

```wolfram
Integrate[x,x,y,GeneratedParameters->C]
(* Output *)
(x^2 y)/(2)+1[y]+2[x]
```

This is the most general mixed antiderivative of the integrand:

```wolfram
D[%,x,y]
(* Output *)
x
```

Integrate over the rectangle from $\{0,0 \}$ to $\{a,b \}$:

```wolfram
Integrate[Exp[x]Cos[y], {x,0,a}, {y,0,b}]
(* Output *)
(-1+ℯ^a) Sin[b]
```

Equivalently:

```wolfram
∫_0^a∫_0^bExp[x]Cos[y]ⅆyⅆx
(* Output *)
(-1+ℯ^a) Sin[b]
```

Integrate in the opposite order:

```wolfram
∫_0^b∫_0^aExp[x]Cos[y]ⅆxⅆy
(* Output *)
(-1+ℯ^a) Sin[b]
```

Combine indefinite and definite integration:

```wolfram
Integrate[(x+y), {y,0,1}, x]
(* Output *)
(x)/(2)+(x^2)/(2)
```

Compute a rational double integral over a rectangular region:

```wolfram
Integrate[x^2y/(x+y),{x,0,3},{y,1,2}]
(* Output *)
(39)/(8)+36 Log[2]-(65 Log[5])/(4)
```

This gives the volume of the shaded region:

```wolfram
Plot3D[(x^2 y)/(x+y),{x,0,3},{y,1,2},Filling->Axis,FillingStyle->Opacity[0.6],Mesh->None]
```

*([Graphics3D])*

Compute a trigonometric double integral over a rectangular region:

```wolfram
Integrate[y Sin[x y],{x,1,4},{y,0,π}]
(* Output *)
0
```

There is as much positive volume (dark gray) as negative (light blue):

```wolfram
Plot3D[{y Sin[x y],0},{x,1,4},{y,0,π},PlotRange -> All, Mesh -> None, PlotStyle -> GrayLevel[0.3]FaceForm, PlotTheme -> SolidGrid, BoxRatios -> 1, Lighting -> Neutral, ViewPoint -> 2.25-1.91.5, Filling -> 1 -> AxisDirective[Specularity[1], Hue[0.58, 1, 0.5, 0.5]]GrayLevel[0.2]]
```

*([Graphics3D])*

Compute a polynomial double integral over the area between two curves:

```wolfram
∫_-2^2∫_2 x^2^4+x^2(8-(3 y)/(3)-(3 x y)/(3)+(x^2 y)/(3)+(x^3 y)/(3))ⅆyⅆx
(* Output *)
(20224)/(315)
```

Visualize the domain of integration and the volume corresponding to the integral:

```wolfram
Plot[{2x^2,4+x^2},{x,-2,2},Filling->{1->{2}},PlotLegends->"Expressions"]
(* Output *)
![image](img/image_016.png)
```

```wolfram
Plot3D[8-(3 y)/(3)-(3 x y)/(3)+(x^2 y)/(3)+(x^3 y)/(3),{x,y}∈ImplicitRegion[-2<x<2&&2 x^2<y<4+x^2,{x,y}],FillingStyle->Opacity[.5],Filling->Axis,PlotRange->All,ViewPoint->{.3,3.1,1.3},Mesh->None]
```

*([Graphics3D])*

Compute a triple integral over a rectangular prism:

```wolfram
Integrate[Cos[z]Exp[x]y,{x,0,1},{y,-1,2},{z,0,3}]
(* Output *)
(3)/(2) (-1+ℯ) Sin[3]
```

Visualize the region of integration:

```wolfram
RegionPlot3D[0<x<1&&-1<y<2&&0<z<3,{x,-0.5,1.5},{y,-1.5,2.5},{z,-0.5,3.5}]
```

*([Graphics3D])*

Compute an algebraic triple integral over a cone-shaped 3D region:

```wolfram
Integrate[Sqrt[x^2+z^2]y,{x,-2,2},{y,x^2,4},{z,-Sqrt[y-x^2],Sqrt[y-x^2]}]//FullSimplify
(* Output *)
(512 π)/(21)
```

Visualize the region of integration:

```wolfram
RegionPlot3D[y>x^2+z^2&&0<y<4,{x,-3,3},{y,-0.5,4.5},{z,-3,3}]
```

*([Graphics3D])*

Integrate a multivariate function over a five-dimensional cube:

```wolfram
Integrate[a b^2 -c/d^e,{a,0,1},{b,1,2},{c,2,3},{d,3,4},{e,4,5}]
(* Output *)
(7)/(6)+(5)/(2) (-ExpIntegralEi[-4 Log[3]]+ExpIntegralEi[-3 Log[3]]+ExpIntegralEi[-4 Log[4]]-ExpIntegralEi[-3 Log[4]])
```

```wolfram
N[%]
(* Output *)
1.1563591256166947
```

Integrate $f(r,\theta)=\frac{sin(\theta)}{r}$ over the unit ball in 4 dimensions:

```wolfram
f[r_,θ_]:=(Sin[θ])/(r)
```

Look up the coordinate ranges for hyperspherical coordinates in [CoordinateChartData](https://reference.wolfram.com/language/ref/CoordinateChartData.html):

```wolfram
CoordinateChartData["Hyperspherical","CoordinateRangeAssumptions",{r,θ,φ,ψ}]
(* Output *)
r>0&&0<θ<π&&0<φ<π&&-π<ψ<=π
```

Also look up the volume factor:

```wolfram
CoordinateChartData["Hyperspherical","VolumeFactor",{r,θ,φ,ψ}]
(* Output *)
r^3 Sin[θ]^2 Sin[φ]
```

Do the integral:

```wolfram
∫_-π^π∫_0^π∫_0^π∫_0^Rf[r,θ](r^3 Sin[θ]^2 Sin[φ])ⅆrⅆθⅆφⅆψ
(* Output *)
(16 π R^3)/(9)
```

#### Region Integrals

Integrate a constant over a unit disk:

```wolfram
Integrate[1,{x,y}∈Disk[]]
(* Output *)
π
```

Enter the integral in typeset form:

```wolfram
∫_{{x,y}∈Disk[]}1
(* Output *)
π
```

Equivalently, integrate over a rectangular region and restrict to a disk using [Boole](https://reference.wolfram.com/language/ref/Boole.html):

```wolfram
Integrate[Boole[x^2+y^2<=1],{x,-1,1},{y,-1,1}]
(* Output *)
π
```

An integral over a unit disk:

```wolfram
Subscript[∫, {x,y}∈Disk[]](Cos[3x]Sin[5y]+x^2)
(* Output *)
(π)/(4)
```

The same integral expressed using [Boole](https://reference.wolfram.com/language/ref/Boole.html):

```wolfram
∫_-1^1∫_-1^1(Cos[3x]Sin[5y]+x^2)Boole[x^2+y^2<=1]ⅆyⅆx
(* Output *)
(π)/(4)
```

The same integral reduced to an iterated integral with bounds depending on the previous variables:

```wolfram
∫_-1^1∫_-Sqrt[1-x^2]^Sqrt[1-x^2](Cos[3x]Sin[5y]+x^2)ⅆyⅆx
(* Output *)
(π)/(4)
```

Plot the integrand over the integration region:

```wolfram
Plot3D[Cos[3x]^2Sin[5y]+1,{x,y}∈Disk[],Filling->0,Mesh->None,FillingStyle->Automatic]
```

*([Graphics3D])*

Express a normal definite integral using region notation:

```wolfram
Integrate[Sin[b x],{x}∈Interval[{0,Pi}]]
(* Output *)
(1-Cos[b π])/(b)
```

Compare with the list notation:

```wolfram
Integrate[Sin[b x],{x,0,Pi}]
(* Output *)
(1-Cos[b π])/(b)
```

With [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) -> [True](https://reference.wolfram.com/language/ref/True.html), assumptions are generated so that the region is non-degenerate:

```wolfram
Integrate[x^2y,{x,y}∈Rectangle[{0,0},{a,b}],GenerateConditions->True]
(* Output *)
(a^3 b^2)/(6)
```

Compare with an iterated integral:

```wolfram
Integrate[x^2y,{x,0,a},{y,0,b},GenerateConditions->True]
(* Output *)
(a^3 b^2)/(6)
```

Integrate over the unit circle:

```wolfram
Subscript[∫, {x,y}∈Circle[]](x^2 y+y^2)/(2)
(* Output *)
(π)/(2)
```

Express the same integral as a one-dimensional integral using polar coordinates:

```wolfram
∫_-Pi^Pi(Cos[θ]^2Sin[θ]+Sin[θ]^2)/(2)ⅆθ
(* Output *)
(π)/(2)
```

Integrate over a sphere of radius $r$:

```wolfram
∫_{{x,y,z}∈Sphere[{x0,y0,z0},r]}(1)/((x-x0)^2+(y-y0)^2+(z-z0)^2)
(* Output *)
4 π
```

Integrate over a finite set of points:

```wolfram
Integrate[1,{x,y,z}∈Point[{{1,2,3},{4,5,6},{7,8,9}}]]
(* Output *)
3
```

Regions can be given as logical combinations of inequalities:

```wolfram
int1=Integrate[Boole[1<x^2-y^2<4&&x y<1&&x>0&&y>0],{x,-∞,∞}, {y,-∞,∞} ]
(* Output *)
(1)/(4) (8 ArcCoth[2+Sqrt[5]]+ArcSinh[2]+4 Log[2]-2 Log[4]-2 Log[(1)/(2) (1+Sqrt[5])])
```

Define the region as an [ImplicitRegion](https://reference.wolfram.com/language/ref/ImplicitRegion.html) and integrate directly over the region:

```wolfram
reg=ImplicitRegion[1<x^2-y^2<4&&x y<1&&x>0&&y>0,{x,y}];
```

```wolfram
int2=Integrate[1,{x,y}∈reg]
(* Output *)
(1)/(2) (-1+2 ArcCsch[2]+ArcSinh[2]+Log[2]-Log[4]+Log[2 (-1+Sqrt[5])]-Log[2 (<|icon -> Root, small -> "1.27", approx -> 1.272019649514069, interp -> Root[-1, -, #, ^, 2, +, #, ^, 4, &, 2, 0], head -> Root, big -> -1-#1^2+#1^4&, degree -> 4, shortDegree -> 4, number -> 2|>+Sqrt[-1+<|icon -> Root, small -> "1.27", approx -> 1.272019649514069, interp -> Root[-1, -, #, ^, 2, +, #, ^, 4, &, 2, 0], head -> Root, big -> -1-#1^2+#1^4&, degree -> 4, shortDegree -> 4, number -> 2|>^2])]+<|icon -> Root, small -> "1.27", approx -> 1.272019649514069, interp -> Root[-1, -, #, ^, 2, +, #, ^, 4, &, 2, 0], head -> Root, big -> -1-#1^2+#1^4&, degree -> 4, shortDegree -> 4, number -> 2|> Sqrt[-1+<|icon -> Root, small -> "1.27", approx -> 1.272019649514069, interp -> Root[-1, -, #, ^, 2, +, #, ^, 4, &, 2, 0], head -> Root, big -> -1-#1^2+#1^4&, degree -> 4, shortDegree -> 4, number -> 2|>^2])
```

The integrals are equivalent:

```wolfram
FullSimplify[int1-int2]
(* Output *)
0
```

Visualize the domain of integration:

```wolfram
RegionPlot[reg]
(* Output *)
![image](img/image_018.png)
```

Integral over a three-dimensional region defined by inequalities:

```wolfram
reg=ImplicitRegion[x+2y+3z<2&&-1<x<y<z<1,{x,y,z}];
```

```wolfram
∫_{x,y,z}∈reg(x^2+2y z)
(* Output *)
(53833)/(151875)
```

Visualize 3D regions using [RegionPlot3D](https://reference.wolfram.com/language/ref/RegionPlot3D.html):

```wolfram
RegionPlot3D[reg,PlotPoints->65,PlotTheme->"FrameGrid"]
(* Output *)
![image](img/image_020.png)
```

Integrate over a solid cone:

```wolfram
Integrate[ (x^2+y^2) Boole[0<=z<=1&&x^2+y^2<=z^2], {x,-∞,∞} ,{y,-∞,∞}, {z,-∞,∞} ]
(* Output *)
(π)/(10)
```

Visualize the domain of integration:

```wolfram
RegionPlot3D[0<=z<=1&&x^2+y^2<=z^2,{x,-1,1},{y,-1,1},{z,0,1}]
(* Output *)
![image](img/image_022.png)
```

Integrate a function with parameters, getting a piecewise result:

```wolfram
Integrate[Boole[ax<y],{x,0,1},{y,0,1}]
(* Output *)
{, {{1, a<=0}, {(2-a)/(2), 0<a<=1}, {(1)/(2 a), True}}}
```

A region with infinitely many components:

```wolfram
reg=ImplicitRegion[Sin[x]>1/2,{{x,0,∞}}];
```

```wolfram
∫_{{x}∈reg}Exp[-x]
(* Output *)
(ℯ^(7 π/6))/(1+ℯ^(2 π/3)+ℯ^(4 π/3))
```

#### Symbolic Features of Integrals

Integrals involving unknown functions are done when possible:

```wolfram
Integrate[f''[x]+2 a f'[x],x]
(* Output *)
2 a f[x]+f^′[x]
```

Differentiate with respect to an endpoint, yielding the fundamental theorem of calculus:

```wolfram
Integrate[f[t], {t,0,x}]
(* Output *)
∫_0^xf[t]ⅆt
```

```wolfram
D[%, x]
(* Output *)
f[x]
```

A generalization:

```wolfram
∂_x∫_a[x]^b[x]f[t]ⅆt
(* Output *)
-f[a[x]] a^′[x]+f[b[x]] b^′[x]
```

Symbolic integrals can be differentiated with respect to parameters:

```wolfram
Integrate[f[x, a], {x, 0, 1}]
(* Output *)
∫_0^1f[x,a]ⅆx
```

```wolfram
D[%,a]
(* Output *)
∫_0^1f^(0,1)[x,a]ⅆx
```

Differentiate with respect to a parameter that appears in both integrand and endpoints:

```wolfram
∂_a∫_p[a]^q[a]f[x,a]ⅆx
(* Output *)
∫_p[a]^q[a]f^(0,1)[x,a]ⅆx-f[p[a],a] p^′[a]+f[q[a],a] q^′[a]
```

Use the [Inactive](https://reference.wolfram.com/language/ref/Inactive.html) form of [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Inactive[Integrate][1/(1+x),x]
(* Output *)
(1)/(1+x)
```

Differentiate:

```wolfram
D[%,x]
(* Output *)
(1)/(1+x)
```

Illustrate indefinite integral identities:

```wolfram
Inactivate[∫x^nⅆx,Integrate]==∫x^nⅆx
(* Output *)
x^n==(x^(1+n))/(1+n)
```

```wolfram
Inactivate[∫(1)/(Sqrt[1-x^2])ⅆx,Integrate]==∫(1)/(Sqrt[1-x^2])ⅆx
(* Output *)
(1)/(Sqrt[1-x^2])==ArcSin[x]
```

Verify the identities starting from the inactive forms:

```wolfram
Activate[{x^n==(x^(1+n))/(1+n),(1)/(Sqrt[1-x^2])==ArcSin[x]}]
(* Output *)
{True,True}
```

Illustrate the basic commutation trick for differentiating under the integral sign:

```wolfram
Inactivate[Integrate[D[f[x,a],a],x]==((D[Integrate[f[x,a],x],a]))]
(* Output *)
Inactive==Inactive
```

```wolfram
Activate[Inactive==Inactive]
(* Output *)
True
```

Compute the [LaplaceTransform](https://reference.wolfram.com/language/ref/LaplaceTransform.html) of an integral:

```wolfram
LaplaceTransform[Integrate[f[u],{u,0,t}],t,s]
(* Output *)
(LaplaceTransform[f[t],t,s])/(s)
```

### Options

#### Assumptions

By default, conditions are generated on parameters that guarantee convergence:

```wolfram
Integrate[x^n,{x,0,1}]
(* Output *)
(1)/(1+n)
```

With [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html), a result valid under the given assumptions is given:

```wolfram
Integrate[x^n,{x,0,1},Assumptions->n>0]
(* Output *)
(1)/(1+n)
```

Manually specify [Assumptions](https://reference.wolfram.com/language/ref/Assumptions.html) to test values outside the automatically generated conditions:

```wolfram
Integrate[Exp[-a x^2],{x,-∞,∞}]
(* Output *)
(Sqrt[π])/(Sqrt[a])
```

This integral is also convergent for purely imaginary $a$:

```wolfram
Integrate[Exp[-a x^2],{x,-∞,∞},Assumptions->Re[a]==0]
(* Output *)
((-a)^(3/4) Sqrt[(π)/(2)] (-1+Sign[a]))/(a^(5/4))
```

Specify assumptions to evaluate a piecewise indefinite integral:

```wolfram
Integrate[Floor[x],x]
(* Output *)
∫Floor[x]ⅆx
```

```wolfram
Integrate[Floor[x],x,Assumptions->0<x<5]
(* Output *)
{, {{0, x<=1}, {-1+x, 1<x<=2}, {-3+2 x, 2<x<=3}, {-6+3 x, 3<x<=4}, {-10+4 x, True}}}
```

#### GenerateConditions

By default, univariate definite integrals generate conditions on parameters that ensure convergence:

```wolfram
Integrate[Exp[-a x^2],{x,0,∞}]
(* Output *)
(Sqrt[π])/(2 Sqrt[a])
```

Generate a result without conditions:

```wolfram
Integrate[Exp[-a x^2],{x,0,∞},GenerateConditions->False]
(* Output *)
(Sqrt[π])/(2 Sqrt[a])
```

Use [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html)->[False](https://reference.wolfram.com/language/ref/False.html) to speed up integration:

```wolfram
Integrate[x^n,{x,-1,b},GenerateConditions->False]//Timing
(* Output *)
{0.139881,((-1)^n+b^(1+n))/(1+n)}
```

```wolfram
Integrate[ x^n,{x,-1,b},GenerateConditions->True]//Timing
(* Output *)
{0.488809,((-1)^n+b^(1+n))/(1+n)}
```

#### GeneratedParameters

By default a particular antiderivative is returned:

```wolfram
Integrate[x,x]
(* Output *)
(x^2)/(2)
```

Specify a value of [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) to obtain the general antiderivative:

```wolfram
Integrate[x,x,GeneratedParameters->C]
(* Output *)
(x^2)/(2)+1
```

One parameter is generated for each indefinite integral:

```wolfram
Integrate[Exp[2x],x,x,GeneratedParameters->C]
(* Output *)
(ℯ^(2 x))/(4)+1+x 2
```

```wolfram
Integrate[Exp[2x],x,x,x,GeneratedParameters->C]
(* Output *)
(ℯ^(2 x))/(8)+1+x 2+x^2 3
```

If the input expression already contains a generated parameter, the next available index will be used:

```wolfram
Integrate[(ℯ^(2 x))/(8)+ 3,x,GeneratedParameters->C]
(* Output *)
(ℯ^(2 x))/(16)+x 3+4
```

For nested integrals with multiple variables, the antiderivative contains arbitrary functions:

```wolfram
Integrate[Exp[2x],x,y,GeneratedParameters->C]
(* Output *)
(1)/(2) ℯ^(2 x) y+1[y]+2[x]
```

```wolfram
Integrate[Exp[2x],x,y,z,GeneratedParameters->C]
(* Output *)
(1)/(2) ℯ^(2 x) y z+1[y,z]+2[x,z]+3[x,y]
```

This is the most general antiderivative of the integrand:

```wolfram
D[%,x,y,z]
(* Output *)
ℯ^(2 x)
```

The value of [GeneratedParameters](https://reference.wolfram.com/language/ref/GeneratedParameters.html) is applied to the index of each generated parameter:

```wolfram
Integrate[Cos[x],x,x,GeneratedParameters->f]
(* Output *)
-Cos[x]+f[1]+x f[2]
```

The value can be a pure function:

```wolfram
Integrate[Cos[x],x,x,GeneratedParameters->(const_#&)]
(* Output *)
-Cos[x]+const_1+x const_2
```

A value of [None](https://reference.wolfram.com/language/ref/None.html) disables generated parameters:

```wolfram
Integrate[Cos[x],x,x,GeneratedParameters->None]
(* Output *)
-Cos[x]
```

#### PrincipalValue

The ordinary Riemann definite integral is divergent:

```wolfram
Integrate[1/x, {x, -1, 2}]
(* Output *)
Integrate
(* Output *)
∫_-1^2(1)/(x)ⅆx
```

The Cauchy principal value integral is finite:

```wolfram
Integrate[1/x, {x, -1, 2}, PrincipalValue->True]
(* Output *)
Log[2]
```

The value is the limit of removing a symmetric region about the singularity:

```wolfram
[Limit]_{a->0^(+)}(∫_-1^-a(1)/(x)ⅆx+∫_a^2(1)/(x)ⅆx)
(* Output *)
Log[2]
```

The ordinary Riemann definite integral is divergent:

```wolfram
Integrate[Tan[x],{x,0,π}]
(* Output *)
Integrate
(* Output *)
∫_0^πTan[x]ⅆx
```

Regularize the divergence at $x=\pi/2$:

```wolfram
Integrate[Tan[x],{x,0,π},PrincipalValue->True]
(* Output *)
0
```

```wolfram
Plot[Tan[x],{x,0,Pi},Filling->Axis, FillingStyle->{<|color -> RGBColor[0.95, 0.627, 0.1425, 0.2]|>,<|color -> RGBColor[0.24, 0.6, 0.8, 0.2]|>}]
```

*([Graphics])*

### Applications

#### The Geometry of Integrals

The integral $\int_{a}^{b}c \mathrm{d}x$ of a constant function $c$ is the signed area of the rectangle of height $c$ and width $b-a$:

```wolfram
∫_1^42ⅆx
(* Output *)
6
```

```wolfram
%==2×(4-1)
(* Output *)
True
```

```wolfram
∫_1^4(-1)ⅆx
(* Output *)
-3
```

```wolfram
%==(-1)×(4-1)
(* Output *)
True
```

Visualize the two rectangles:

```wolfram
Plot[{2,-1},{x,1,4},Filling->Axis]
```

*([Graphics])*

The integral $\int_{a}^{b}f(x)\mathrm{d}x$ of a piecewise-constant function is the sum of the signed areas of the rectangles defined by its plot:

```wolfram
f[x_]:={, {{0, x<0}, {1, 0<=x<1}, {-2, 1<=x<3}, {5, 3<=x}}}
```

```wolfram
∫_0^4f[x]ⅆx
(* Output *)
2
```

```wolfram
%==1×(1-0)+(-2)(3-1)+5×(4-3)
(* Output *)
True
```

Visualize the rectangles:

```wolfram
Plot[f[x],{x,0,4},Filling->Axis,FillingStyle->{<|color -> RGBColor[0.95, 0.627, 0.1425, 0.2]|>,<|color -> RGBColor[0.24, 0.6, 0.8, 0.2]|>}]
```

*([Graphics])*

The integral $\int_{a}^{b}f(x)\mathrm{d}x$ of a general function is the signed area between its plot and the horizontal axis:

```wolfram
f[x_]:= x^2+1
```

```wolfram
Integrate[f[x],{x,0,2}]
(* Output *)
(14)/(3)
```

```wolfram
Plot[f[x],{x,0,2},Filling->Axis,AxesOrigin->{0,0}]
```

*([Graphics])*

This can be related to the piecewise-constant case by considering rectangles defined by its plot:

```wolfram
rectangles[f_,a_,b_,n_]:={Opacity[0.75],StandardOrange,N@Table[Rectangle[{a+k(b-a)/n,0},{a+(k+1)(b-a)/n,f[a+k(b-a)/n]}],{k,0,n-1}]}
```

For `*n*==5 on the interval `[0,2]`, the rectangles are the following:

```wolfram
Plot[f[x],{x,0,2},Epilog->rectangles[f,0,2,5],AxesOrigin->{0,0}]
```

*([Graphics])*

The area of these rectangles defines a Riemann sum that approximates the area under the curve:

```wolfram
RiemannSum[f_,a_,b_,n_]:=
Sum[f[(a+(b-a)k/n)],{k,0,n-1}] (b-a)/n
```

```wolfram
RiemannSum[f,0,2,n]
(* Output *)
(2 (2-6 n+7 n^2))/(3 n^2)
```

Using [DiscreteLimit](https://reference.wolfram.com/language/ref/DiscreteLimit.html) to obtain the exact answer as $n \to \infty$ gives the same answer as [Integrate](https://reference.wolfram.com/language/ref/Integrate.html) did:

```wolfram
DiscreteLimit[%,n->∞]
(* Output *)
(14)/(3)
```

Visualize the process for this function as well as three others:

```wolfram
label[f_,a_,b_,n_]:=StringForm["\`A_curve == ``\\ \\ \\ A_rects == ``",NumberForm[N@∫_a^bf[x]ⅆx,{4,3}],NumberForm[N@RiemannSum[f,a,b,n],{4,3}]]
```

```wolfram
With[{a=0,b=2},Manipulate[Plot[fn[x],{x,a,b},Filling->Axis,AxesOrigin->{0,0},Epilog->rectangles[fn,a,b,2^n],PlotLabel->label[fn,a,b,2^n]],{{n,1,Dynamic@StringForm["n == ``",2^n]}, 1,12,1},{{fn,f},{f,CubeRoot,Sin,Abs[2#-1]&->HoldForm@2x-1}},SaveDefinitions->True]]
```

The Fundamental Theorem of Calculus relates a function to its integral from a fixed lower limit to a variable upper limit:

```wolfram
f[x_] := Sqrt[x^2+1]
```

Consider the definite integral of the this from from $0$ to $x$:

```wolfram
g[x_]=Integrate[f[t], {t, 0, x},Assumptions->x>0]
(* Output *)
(1)/(2) (x Sqrt[1+x^2]+ArcSinh[x])
```

The Fundamental Theorem of Calculus states that $g^{'}(x)=f(x)$:

```wolfram
Simplify[g'[x]==f[x]]
(* Output *)
True
```

This can be seen from the limit definition of derivative:

```wolfram
Limit[(g[x+h]-g[x])/(h),h->0]
(* Output *)
Sqrt[1+x^2]
```

Note that $g(x+h)-g(x)=\int_{x}^{x+h}f(t)\mathrm{d}t$ is an area consisting of a rectangle of height $f(x)$ and width $h$ plus a small correction that vanishes as $h \to 0$, as illustrated by the following table for $x=1$:

```wolfram
TableForm[Table[{h,(g[1+h]-g[1])/(h),(g[1+h]-g[1])/(h)-f[1]},{h,{.1,.05,.001,0.0005,0.00001,0.00000005}}],TableHeadings->{{},{h,HoldForm[(g[1+h]-g[1])/(h)],HoldForm[(g[1+h]-g[1]-f[1]h)/(h)]}}]
(* Output *)
{{, h, (g[1+h]-g[1])/(h), (g[1+h]-g[1]-f[1] h)/(h)}, {"", 0.1, 1.4501367131676601, 0.035923150794564984}, {"", 0.05, 1.4320358248151699, 0.017822262442074743}, {"", 0.001, 1.4145671746670363, 0.00035361229394115234}, {"", 0.0005, 1.4143903537968683, 0.00017679142377313717}, {"", 0.00001, 1.414217097917003, 3.5355439078621487×10^-6}, {"", 5.×10^-8, 1.4142135773553832, 1.498228807683688×10^-8}}
```

Hence, the limit can be seen geometrically to equal $f(x)$, as illustrated in the following visualization:

```wolfram
title[u_,k_]:=Column[HoldForm[Style[1, /, h, ), Style[g[x, +, h], -, g[x], FontColor -> DefaultPlotColors, [, 2]], ScriptLevel -> 0]]==g[u, +, k]-g[u])/kHoldForm[f[x]]==Sqrt[1, +, u, ^, 2], Spacings -> 0.4, Alignment -> ==]
```

```wolfram
labels[u_,k_]:=Text[HoldForm[g[x], ==, Integrate[f[t], t0x]], u/20.5, BaseStyle -> GrayLevel[1, 0.5]]Text[\`(x, f(x)), uSqrt[1, +, u, ^, 2], 1.2-1]Text[\`(x + h, f(x)), k+uSqrt[1, +, u, ^, 2], -1.21]PointSize[Large]Line[u0uSqrt[1, +, u, ^, 2]k+uSqrt[1, +, u, ^, 2]k+u0]Text[Style[f, 16], 2.52.9]Point[uSqrt[1, +, u, ^, 2]k+uSqrt[1, +, u, ^, 2]]
```

```wolfram
Manipulate[Plot[{Piecewise[{{f[s],0<=s<=u},{Undefined,True}}],Piecewise[{{f[s],u<=s<=u+k},{Undefined,True}}],f[s]},{s,0,3},PlotLabel->title[u,k],Epilog->labels[u,k],Filling -> 1 -> 02 -> 0, FillingStyle -> Opacity[1], AxesOrigin -> 00, PlotStyle -> NoneNoneLightDarkSwitched[GrayLevel[0]], ImageSize -> Medium, BaseStyle -> 15, PerformanceGoal -> Quality],{{u,1.5,Style[x,Italic,14]},1,2,Appearance->{"Labeled"},LabelStyle->14},{{k,.25,Style[h,Italic,14]},0.00001,.5,Appearance->"Labeled",LabelStyle->14},SaveDefinitions->True]
```

Integrate a discrete set of data with [Interpolation](https://reference.wolfram.com/language/ref/Interpolation.html):

```wolfram
data=Table[{x,Sin[(1)/(x+1/2)]},{x,0,10,0.25}];
```

```wolfram
Integrate[Interpolation[data][x],{x,0,10}]
(* Output *)
2.7408984919556962
```

Compare with the exact area under the curved that was interpolated:

```wolfram
Integrate[Sin[(1)/(x+1/2)],{x,0,10}]//N
(* Output *)
2.7432473941510094
```

Visualize the two areas:

```wolfram
{ListLinePlot[data,Filling->Axis,ImageSize->Small,InterpolationOrder->0],Plot[Interpolation[data][x],{x,0,10},Filling->Axis,ImageSize->Small]}
(* Output *)
{[Graphics],[Graphics]}
```

#### Area Between Curves

Compute the area under the curve of $f(x)=\sqrt{9-x^{2}}+1$ from $-3$ to $0$:

```wolfram
Plot[1+Sqrt[9-x^2],{x,-3,0},Filling->Axis,AxesOrigin->{0,0}]
```

*([Graphics])*

```wolfram
Integrate[1+Sqrt[9-x^2],{x,-3,0}]
(* Output *)
3+(9 π)/(4)
```

Find the area under the curve of $f(x)=sin^{2}(x) cos^{4}(x)$ from $-\pi$ to $\pi$:

```wolfram
Plot[Sin[x]^2 Cos[x]^4,{x,-π,π},Filling->Axis]
```

*([Graphics])*

```wolfram
Integrate[Sin[x]^2 Cos[x]^4,{x,-π,π}]
(* Output *)
(π)/(8)
```

Determine the total area enclosed between of $f(x)=sin(x)$ and the $x$-axis:

```wolfram
Plot[Sin[x],{x,-π,π},Filling->Axis]
```

*([Graphics])*

The total area is given by the integral of the absolute value:

```wolfram
Integrate[RealAbs[Sin[x]],{x,-π,π}]
(* Output *)
4
```

Equivalently compute this as the sum of two integrals of the difference between the top and bottom:

```wolfram
∫_-π^0(0-Sin[x])ⅆx+∫_0^π(Sin[x]-0)ⅆx
(* Output *)
4
```

Compute the area between $f(x)=5 x-x^{2}$ and $g(x)=x$ from $0$ to $4$:

```wolfram
Plot[{5x-x^2,x},{x,0,4},PlotLegends->"Expressions",Filling->{1->{2}}]
(* Output *)
![image](img/image_024.png)
```

```wolfram
Integrate[(5x-x^2)-(x),{x,0,4}]
(* Output *)
(32)/(3)
```

Find the area enclosed by $f(x)=\frac{2}{x^{4}+1}$ and $g(x)=x^{2}$:

```wolfram
pts=Solve[y==2/(1+x^4)&&y==x^2,{x,y},Reals]
(* Output *)
{{x->-1,y->1},{x->1,y->1}}
```

Since $f(0)>g(0)$, $f$ will be above $g$ in the interval of interest and the area will equal:

```wolfram
Integrate[2/(x^4+1)-x^2,{x,-1,1}]
(* Output *)
-(2)/(3)+(π+2 ArcCoth[Sqrt[2]])/(Sqrt[2])
```

```wolfram
%//N
(* Output *)
2.8012252826929775
```

Visualize the region of interest and the two functions:

```wolfram
Show[RegionPlot[x^2<y<(2)/(x^4+1),{x,-1,1},{y,0,2}],Plot[{2/(1+x^4),x^2},{x,-2,2},PlotLegends->"Expressions"],PlotRange->All]
(* Output *)
![image](img/image_026.png)
```

Compute the area enclosed by $y=cos(x)$ and $y=1-\frac{2 x}{\pi}$:

```wolfram
pts=Solve[y==Cos[x]&&y==1-2x/π&&0<=x<=π,{x,y}]
(* Output *)
{{x->0,y->1},{x->(π)/(2),y->0},{x->π,y->-1}}
```

Find the area as the integral of the absolute value of the difference over the entire interval:

```wolfram
Integrate[RealAbs[Cos[x]-(1-2x/π)],{x,0,π}]//Expand
(* Output *)
2-(π)/(2)
```

Visualize the two functions and the area between them:

```wolfram
Show[RegionPlot[y<=Cos[x]&&y>=1-2x/π || y>=Cos[x]&&y<=1-2x/π,{x,0,π},{y,-2,2},PlotPoints->60,PlotStyle->Opacity[.75]],Plot[{Cos[x],1-2x/π},{x,-1,4},PlotTheme->"Detailed"],Graphics[{PointSize[Large],Point[{x,y}/.pts]}],PlotRange->All,Axes->True]
(* Output *)
![image](img/image_028.png)
```

Use the plot to split the integral into two equivalent integrals with no absolute value:

```wolfram
∫_0^(π)/(2)(Cos[x]-(1-(2 x)/(π)))ⅆx+∫_(π)/(2)^π((1-(2 x)/(π))-Cos[x])ⅆx
(* Output *)
2-(π)/(2)
```

To compute the area enclosed by $y=\frac{1}{x^{2}}$, $y=x$, and $y=\frac{x}{8}$, first find the points of intersection:

```wolfram
pt1=Solve[y==1/x^2&&y==x,{x,y},Reals]
(* Output *)
{{x->1,y->1}}
```

```wolfram
pt2=Solve[y==1/x^2&&y==(1/8)x,{x,y},Reals]
(* Output *)
{{x->2,y->(1)/(4)}}
```

```wolfram
pt3=Solve[y==(1/8)x&&y==x,{x,y},Reals]
(* Output *)
{{x->0,y->0}}
```

Visualize the three curves over an area containing the points:

```wolfram
ContourPlot[{y==1/x^2,y==x,y==(1/8)x},{x,-1,3},{y,-1,2},PlotLegends->"Expressions",Axes->True,AspectRatio->Automatic,Epilog->{PointSize[Large],Point[{x,y}/.Join[pt1,pt2,pt3]]}]
```

*([Graphics])*

From the plot, it is clear $y$ is above the line $y=\frac{x}{8}$ and below the other two curves:

```wolfram
reg = ImplicitRegion[y<=1/x^2&&y<=x&&y>=(1/8)x,{x,y}]
(* Output *)
ImplicitRegion[y<=(1)/(x^2)&&y<=x&&y>=(x)/(8),{x,y}]
```

```wolfram
RegionPlot[reg,AspectRatio->Automatic]
(* Output *)
![image](img/image_030.png)
```

Area can be found using two integrals, one for each "top function":

```wolfram
Integrate[x-(1/8)x,{x,0,1}]+Integrate[1/x^2-(1/8)x,{x,1,2}]
(* Output *)
(3)/(4)
```

This can be reduced to a single integral using [Min](https://reference.wolfram.com/language/ref/Min.html):

```wolfram
Integrate[Min[x,1/x^2]-(1/8)x,{x,0,2}]
(* Output *)
(3)/(4)
```

Compare with the answer returned by [Area](https://reference.wolfram.com/language/ref/Area.html):

```wolfram
Area[reg]
(* Output *)
(3)/(4)
```

#### Regions of Revolution

Compute the volume enclosed when $y=sin^{2}(x)$ for $\pi/2<x<\pi$ is rotated about the $x$-axis:

```wolfram
Integrate[π(Sin[x]^2)^2,{x,0,π}]
(* Output *)
(3 π^2)/(8)
```

Visualize the solid:

```wolfram
RevolutionPlot3D[Sin[x]^2,{x,0,π},RevolutionAxis->"X"]
```

*([Graphics3D])*

Use cylindrical shells to find the volume enclosed when $y=x (x-1)^{2}$, $0<x<1$, is rotated about the $y$-axis:

```wolfram
Integrate[2π x(x(x-1)^2),{x,0,1}]
(* Output *)
(π)/(15)
```

Visualize the solid, adding the cap at $y=0$:

```wolfram
Show[RevolutionPlot3D[{x(x-1)^2},{x,0,1}],RevolutionPlot3D[0,{x,0,1}]]
```

*([Graphics3D])*

Find the volume of the region formed by rotating the area between $f(x)=x$ and $g(x)=x^{2}$ about the $y$-axis:

```wolfram
f[x_]:=x
g[x_]:=x^2
```

Find where the curves intersect:

```wolfram
Solve[f[x]==g[x],x]
(* Output *)
{{x->0},{x->1}}
```

Between these two values of $x$, $f$ is above $g$:

```wolfram
Plot[{g[x],f[x]},{x,0,1},Filling->{2->{1}}]
```

*([Graphics])*

Visualize the volume:

```wolfram
RevolutionPlot3D[{{x,f[x]},{x,g[x]}},{x,0,1},PlotStyle->Opacity[.6]]
(* Output *)
![image](img/image_032.png)
```

Integrate cylindrical shells of height $f(x)-g(x)$ and circumference $2 \pi x$ to find the volume:

```wolfram
∫_0^12π x(f[x]-g[x])ⅆx
(* Output *)
(π)/(6)
```

Determine the volume the region above $f(x)=csc(x)$ and below $y=2$ for $0 \leq x \leq \pi$, rotated about the $x$-axis:

```wolfram
f[x_]:=Csc[x]
```

Find where the curves intersect, adding the constraint on the range of $x$:

```wolfram
Solve[f[x]==2&&0<=x<=π,x]
(* Output *)
{{x->(π)/(6)},{x->(5 π)/(6)}}
```

The relevant range of $x$ values is between these two points:

```wolfram
Plot[{2,f[x]},{x,0,π},Filling->{2->{1}},PlotRange->{0,2}]
```

*([Graphics])*

Visualize the volume:

```wolfram
Show[RevolutionPlot3D[f[x],{x,(π)/(6),(5 π)/(6)},RevolutionAxis->{1, 0, 0}],Graphics3D[{Opacity[.5],Cylinder[{{0,0,0},{π,0,0}},2]}],ViewPoint->{2,-2.5,1},PlotRange->{{π/6,5π/6},{-2,2},{-2,2}}]
```

*([Graphics3D])*

Integrate washers of area $\pi(2^{2}-f(x)^{2})$ to find the volume:

```wolfram
∫_(π)/(6)^(5 π)/(6)π(2^2-f[x]^2)ⅆx
(* Output *)
-2 Sqrt[3] π+(8 π^2)/(3)
```

```wolfram
%//N
(* Output *)
15.436148884166313
```

Compute the surface area when $y=x^{3}$ for $0<x<1$ is rotated about the $y$-axis:

```wolfram
RevolutionPlot3D[x^3,{x,0,1}]
```

*([Graphics3D])*

Apply the formula of the infinitesimal width of each strip:

```wolfram
w=Sqrt[1+D[x^3,x]^2]
(* Output *)
Sqrt[1+9 x^4]
```

Multiply the width by the circumference $2 \pi x$ of each circle and integrate:

```wolfram
∫_0^12 π x wⅆx
(* Output *)
(1)/(6) π (3 Sqrt[10]+ArcTanh[(3)/(Sqrt[10])])
```

Find the area when $y=\sqrt{x^{2}+1}$ for -$\sqrt{3}<x<\sqrt{3}$ is rotated about the $x$-axis:

```wolfram
RevolutionPlot3D[Sqrt[x^2+1],{x,-Sqrt[3],Sqrt[3]},RevolutionAxis->"X"]
```

*([Graphics3D])*

The infinitesimal width of each strip is given by the following:

```wolfram
w=Sqrt[1+D[Sqrt[x^2+1],x]^2]
(* Output *)
Sqrt[1+(x^2)/(1+x^2)]
```

Multiplying the width by the circumference $2 \pi \sqrt{x^{2}+1}$ and integrating yields the answer:

```wolfram
∫_-Sqrt[3]^Sqrt[3]2 π Sqrt[x^2+1] wⅆx
(* Output *)
π (2 Sqrt[21]+Sqrt[2] ArcSinh[Sqrt[6]])
```

Determine the surface area when $y=\frac{1}{x}$ for $1 \leq x \leq 3$ is rotated about line $x=4$:

```wolfram
f[x_]:=1/x
```

The infinitesimal width of each strip is given by the following:

```wolfram
w=Sqrt[1+D[f[x],x]^2]
(* Output *)
Sqrt[1+(1)/(x^4)]
```

Since $x<4$ for the curve in question, each strip has radius $4-x$ and width $w \mathrm{d}x$:

```wolfram
∫_1^32π(4-x)wⅆx
(* Output *)
8 Sqrt[π] Gamma[(3)/(4)]^2+π (Sqrt[2]-Sqrt[82]-ArcTanh[Sqrt[2]]+ArcTanh[Sqrt[82]]-(656)/(3) Sqrt[82] Hypergeometric2F1[1,(5)/(4),(3)/(4),-81])
```

Find the numerical approximation of this value:

```wolfram
%//N//Chop
(* Output *)
27.524486788827765
```

Visualize the surface using modified cylindrical coordinates based on the line $x=4$, $z=0$:

```wolfram
ParametricPlot3D[{Cos[θ](x-4)+4,1/x,Sin[θ](x-4)},{x,1,3},{θ,0,2π},BoxRatios -> 1GoldenRatio^-11, ViewVertical -> 010, ViewPoint -> -22.5-1.5]
```

*([Graphics3D])*

#### Arc Length, Surface Area, and Volume

Compute the arc length of the plot $f(x)=\frac{x^{4}}{8}+\frac{1}{4 x^{2}}$ from $1$ to $2$:

```wolfram
f[x_]:=(x^4)/(8)+(1)/(4 x^2)
```

```wolfram
Plot[f[x],{x,1,2},PlotRange->{0,2}]
```

*([Graphics])*

Apply the formula for infinitesimal arc length:

```wolfram
ds=Sqrt[1+D[f[x],x]^2]//Simplify
(* Output *)
(1)/(2) Sqrt[2+(1)/(x^6)+x^6]
```

Integrate to find the arc length:

```wolfram
Integrate[ds,{x,1,2}]
(* Output *)
(33)/(16)
```

Compare with the answer returned by [ArcLength](https://reference.wolfram.com/language/ref/ArcLength.html):

```wolfram
ArcLength[f[x],{x,1,2}]
(* Output *)
(33)/(16)
```

Compute the arc length of the plot $f(x)=\sqrt{x-x^{2}}+sin^{-1}(\sqrt{x})$ from $0$ to $1$:

```wolfram
f[x_]:=Sqrt[x-x^2]+ArcSin[Sqrt[x]]
```

```wolfram
Plot[f[x],{x,0,1}]
```

*([Graphics])*

Apply the formula for infinitesimal arc length:

```wolfram
ds=Sqrt[1+D[f[x],x]^2]//Simplify
(* Output *)
Sqrt[(1)/(x)]
```

Integrate to find the arc length:

```wolfram
Integrate[ds,{x,0,1}]
(* Output *)
2
```

Compare with the answer returned by [ArcLength](https://reference.wolfram.com/language/ref/ArcLength.html):

```wolfram
ArcLength[f[x],{x,0,1}]
(* Output *)
2
```

Length of a parametrically defined circle:

```wolfram
{x,y}={2Cos[u],2Sin[u]};
```

The relevant parameter range is $0$ to $2 \pi$:

```wolfram
ParametricPlot[{x,y},{u,0,2π}]
```

*([Graphics])*

The infinitesimal arc length is constant:

```wolfram
ds=Sqrt[D[x,u]^2+D[y,u]^2]//Simplify
(* Output *)
2
```

Integrate to find the total arc length:

```wolfram
Integrate[ds,{u,0,2π}]
(* Output *)
4 π
```

Compare with the answer returned by [ArcLength](https://reference.wolfram.com/language/ref/ArcLength.html):

```wolfram
ArcLength[{x,y},{u,0,2π}]
(* Output *)
4 π
```

Length of a 3D parametrically defined ellipse:

```wolfram
c[u_]={1,0,0}2Cos[u]+Normalize[{0,1,1}]Sin[u]
(* Output *)
{2 Cos[u],(Sin[u])/(Sqrt[2]),(Sin[u])/(Sqrt[2])}
```

Visualize the ellipse:

```wolfram
ParametricPlot3D[c[u],{u,0,2π}]
```

*([Graphics3D])*

The infinitesimal arc length is non-constant:

```wolfram
ds=Sqrt[c'[u].c'[u]]
(* Output *)
Sqrt[Cos[u]^2+4 Sin[u]^2]
```

Integrate to find the total arc length:

```wolfram
Integrate[ds,{u,0,2π}]
(* Output *)
4 EllipticE[-3]
```

```wolfram
%//N
(* Output *)
9.688448220547675
```

Compare with the answer returned by [ArcLength](https://reference.wolfram.com/language/ref/ArcLength.html):

```wolfram
ArcLength[c[u],{u,0,2π}]
(* Output *)
4 EllipticE[-3]
```

Find the surface area of the plot $f(x,y)=x^{2}+y^{2}$ over the rectangle $[0,1]\times[0,1]$:

```wolfram
f[x_,y_]:=x^2+y^2
```

```wolfram
Plot3D[f[x,y],{x,0,1},{y,0,1},BoxRatios->Automatic]
```

*([Graphics3D])*

Apply the formula for infinitesimal surface area of a plot:

```wolfram
dA=Sqrt[1+D[f[x,y],x]^2+D[f[x,y],y]^2]//Simplify
(* Output *)
Sqrt[1+4 x^2+4 y^2]
```

Integrate to find the arc length:

```wolfram
Integrate[dA,{x,0,1},{y,0,1}]
(* Output *)
(1)/(24) (24-6 ArcCot[2]+ArcTan[(4)/(3)]+14 Log[5])
```

```wolfram
%//N
(* Output *)
1.8615641807530907
```

Compare with the answer returned by [Area](https://reference.wolfram.com/language/ref/Area.html):

```wolfram
Area[f[x,y],{x,0,1},{y,0,1}]
(* Output *)
(1)/(24) (24-6 ArcCot[2]+ArcTan[(4)/(3)]+14 Log[5])
```

Find the area of the surface $s(u,v)=\{cos(u) (cos(v)+2),sin(u) (cos(v)+2),sin(v)\}$ where $0 \leq u,v \leq 2 \pi$:

```wolfram
s[u_,v_]:={(2+Cos[v])Cos[u],(2+Cos[v])Sin[u],Sin[v]}
```

This surface is a torus:

```wolfram
ParametricPlot3D[s[u,v],{u,0,2π},{v,0,2π}]
```

*([Graphics3D])*

Apply the formula for infinitesimal surface area of a parametric surface:

```wolfram
dA=Simplify[Norm[∂_us[u,v]⨯∂_vs[u,v]],{u,v}∈Reals]
(* Output *)
2+Cos[v]
```

Integrate to find the total surface area:

```wolfram
Integrate[dA,{u,0,2π},{v,0,2π}]
(* Output *)
8 π^2
```

Compare with the answer returned by [Area](https://reference.wolfram.com/language/ref/Area.html):

```wolfram
Area[s[u,v],{u,0,2π},{v,0,2π}]
(* Output *)
8 π^2
```

Find the volume of the following parametric region, where $0 \leq r \leq 1$, $0 \leq u,v \leq 2 \pi$:

```wolfram
ρ[r_,u_,v_]:={(2+r Cos[v])Cos[u],(2+r Cos[v])Sin[u],Sin[v]}
```

This region is a solid torus:

```wolfram
ParametricPlot3D[ρ[1,u,v],{u,0,2Pi},{v,0,2Pi}]
```

*([Graphics3D])*

Compute the Jacobian determinant:

```wolfram
j=Simplify[RealAbs[Det[Grad[ρ[r,u,v],{r,u,v}]]],0<=r<=1&&0<=v<=2π]
(* Output *)
Cos[v]^2 (2+r Cos[v])
```

Integrate to find the volume:

```wolfram
∫_0^1∫_0^2 π∫_0^2πjⅆuⅆvⅆr
(* Output *)
4 π^2
```

Compare with the answer returned by [Volume](https://reference.wolfram.com/language/ref/Volume.html):

```wolfram
Volume[ρ[r,u,v],{r,0,1},{u,0,2Pi},{v,0,2Pi}]
(* Output *)
4 π^2
```

Find the volume of the following parametric region, where $0 \leq r \leq 1$, $0 \leq \theta \leq 2 \pi$, and $0 \leq \psi \leq \pi$:

```wolfram
ρ[r_,θ_,ψ_]:={3r Cos[ψ],2r Sin[θ]Sin[ψ],r Cos[θ]Sin[ψ]}
```

The region is an ellipsoid:

```wolfram
ParametricPlot3D[ρ[1,θ,ψ],{θ,0,2Pi},{ψ,0,Pi}]
(* Output *)
![image](img/image_034.png)
```

Compute the Jacobian determinant:

```wolfram
j=Simplify[RealAbs[Det[Grad[ρ[r,θ,ψ],{r,θ,ψ}]]],0<=r<=1&&0<=ψ<=π]
(* Output *)
6 r^2 Sin[ψ]
```

Integrate to find the volume:

```wolfram
∫_0^1∫_0^2 π∫_0^πjⅆψⅆθⅆr
(* Output *)
8 π
```

Compare with the answer returned by [Volume](https://reference.wolfram.com/language/ref/Volume.html):

```wolfram
Volume[ρ[r,θ,ψ],{r,0,1},{θ,0,2Pi},{ψ,0,Pi}]
(* Output *)
8 π
```

#### Line Integrals

Compute the line integral $\oint f(s)\mathrm{d}s$ of $f(x,y)=x+y^{2}$ over the origin-centered ellipse with semi-major axes $2$ and $1$:

```wolfram
f[{x_,y_}]:=x +y^2
```

```wolfram
reg=Region[Circle[{0,0},{2,1}],Axes->True]
```

*([Graphics])*

Parameterize the ellipse:

```wolfram
c[t_]:={2Cos[t],Sin[t]}
```

```wolfram
ParametricPlot[c[t],{t,0,2Pi}]
```

*([Graphics])*

Perform the integral using the fact that $\mathrm{d}s=\left\|c^{'}(t)\right\|\mathrm{d}t$:

```wolfram
∫_0^2 πf[c[t]] Norm[c^′[t]]ⅆt
(* Output *)
(4)/(9) (7 EllipticE[-3]-4 EllipticK[-3])
```

```wolfram
%//N
(* Output *)
5.618556929315176
```

Compare the direct integral over the ellipse:

```wolfram
∫_{[DoubleStruckX]∈reg}f[[DoubleStruckX]]
(* Output *)
(4)/(9) (7 EllipticE[-3]-4 EllipticK[-3])
```

Calculate the closed line integral $\oint \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{r}$ of $v(x,y)=(x y,3 y^{2})$ over the following parametric curve:

```wolfram
c[t_]:={Cos[t],Sin[2t]}
```

The curve forms an infinity figure, traversed from red to purple as shown in the following plot:

```wolfram
Show[ParametricPlot[{Cos[x],Sin[2x]},{x,0,2Pi},PlotTheme->"Detailed",ColorFunction->(Hue[0.8#3]&)]]
```

*([Graphics])*

Define the vector field $\overset{⇀}{v}$:

```wolfram
v[{x_,y_}]:={x y,3y^2}
```

Perform the calculation using the definition $\oint \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{r}=\oint \overset{⇀}{v}\cdot \overset{⇀}{c^{'}}(t)\mathrm{d}t$:

```wolfram
Integrate[v[c[t]].c'[t],{t,0,2π}]
(* Output *)
-(π)/(2)
```

To calculate `∫x^(4)dx+x yⅆy` over the triangle with vertices $(0,0)$, $(1,0)$, and $(0,1)$, define the associated vector field:

```wolfram
v[{x_,y_}]:={x^4,x y}
```

Parametrize the triangle using a piecewise-linear parametrization:

```wolfram
c[t_]={, {{{0,0}+t({1,0}-{0,0}), 0<=t<=1}, {{1,0}+(t-1)({0,1}-{1,0}), 1<=t<=2}, {{0,1}+(t-2)({0,0}-{0,1}), 2<=t<=3}, {Indeterminate, True}}}
(* Output *)
{, {{{t,0}, 0<=t<=1}, {{2-t,-1+t}, 1<=t<=2}, {{0,3-t}, 2<=t<=3}, {Indeterminate, True}}}
```

The parametrization is oriented counter-clockwise:

```wolfram
ParametricPlot[c[t],{t,0,3},PlotTheme->"FrameGrid",ColorFunction->(Hue[0.8#3]&)]
```

*([Graphics])*

Compute the line integral from the definition $\oint \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{r}=\oint \overset{⇀}{v}\cdot \overset{⇀}{c^{'}}(t)\mathrm{d}t$:

```wolfram
Integrate[v[c[t]].c'[t],{t,0,3}]
(* Output *)
(1)/(6)
```

Calculate the work done by the force $F(x,y,z)=-\frac{G m M (x \overset{^}{i} y \overset{^}{j}+z \overset{^}{k})}{(x^{2}+y^{2}+z^{2})^{3/2}}$ as a particle takes the following path from $(a,0,0)$, $a>0$, to $(0,b,0)$, $b>0$:

```wolfram
c[t_]:={a Cos[t],b Sin[t], t(t-Pi/2)}
```

Define the force field as function from points to vectors:

```wolfram
F[{x_,y_,z_}]:=-(G M m{x,y,z})/((x^2+y^2+z^2)^(3/2))
```

The work done is the line integral $\int \overset{⇀}{F}\cdot \mathrm{d}\overset{⇀}{r}$:

```wolfram
Integrate[F[c[t]].c'[t],{t,0,Pi/2}, Assumptions->{a>0,b>0}]//Expand
(* Output *)
-(G m M)/(a)+(G m M)/(b)
```

Find a potential function for the following vector field:

```wolfram
v[{x_,y_,z_}]:={x^2,y^2,z^2}
```

This is possible because the vector field is conservative:

```wolfram
Curl[v[{x,y,z}],{x,y,z}]
(* Output *)
{0,0,0}
```

Define a family of straight-line curves that go from the origin at time $0$ to $(x,y,z)$ at time $1$:

```wolfram
c_{x_,y_,z_}[t_]:=t{x,y,z}
```

Let $\phi(x,y,z)$ be the line integral of $v$ from the origin to the point $(x,y,z)$:

```wolfram
φ[{x_,y_,z_}]=∫_0^1v[c_{x,y,z}[t]].D[c_{x,y,z}[t],t]ⅆt
(* Output *)
(1)/(3) (x^3+y^3+z^3)
```

Verify that $\phi$ is a potential function for $v$ using Grad

```wolfram
Grad[φ[{x,y,z}],{x,y,z}]==v[{x,y,z}]
(* Output *)
True
```

Use Green's Theorem to find the area of the area enclosed by the following curve:

```wolfram
c[u_]:={Cos[u]+1/7Cos[7u+π/3],Sin[u]+1/7Sin[7u]};
```

```wolfram
ParametricPlot[c[u],{u,0,2π}]
```

*([Graphics])*

The following vector-field has a two-dimensional [Curl](https://reference.wolfram.com/language/ref/Curl.html) of $1$:

```wolfram
v[{x_,y_}]:={0,x}
```

```wolfram
Curl[v[{x,y}],{x,y}]
(* Output *)
1
```

Apply Green's theorem in the form $\int \int 1 \mathrm{d}A=\oint \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{r}=\oint \overset{⇀}{v}\cdot \overset{⇀}{c^{'}}(u)\mathrm{d}u$ to compute the area:

```wolfram
∫_0^2 πv[c[u]].c^′[u]ⅆu
(* Output *)
(15 π)/(14)
```

#### Surface and Volume Integrals

Use Green's Theorem to compute $\int((3 y-e^{sin(x)})\mathrm{d}x+(7 x+\sqrt{y^{4}+1})\mathrm{d}y$ over the circle centered at the origin with radius 3:

```wolfram
v[x_,y_]:={3y^3 - Exp[Sin[x]],7x^2+Sqrt[y^4+1]}
```

Visualize the vector field and circle for the line integral:

```wolfram
VectorPlot[v[x,y],{x,-3,3},{y,-3,3},Epilog->Circle[{0,0},3],PlotRangePadding->None]
```

*([Graphics])*

The circulation of the vector field can be computed using [Curl](https://reference.wolfram.com/language/ref/Curl.html):

```wolfram
circ=Curl[v[x,y],{x,y}]
(* Output *)
14 x-9 y^2
```

Integrate over the interior of the circle:

```wolfram
∫_-3^3∫_-Sqrt[9-y^2]^Sqrt[9-y^2]circⅆxⅆy
(* Output *)
-(729 π)/(4)
```

Perform the integral using region notation:

```wolfram
∫_{{x,y}∈Disk[{0,0},3]}circ
(* Output *)
-(729 π)/(4)
```

Compute the integral over the unit sphere of $f(x,y,z)=(x+y+z)^{4}$:

```wolfram
f[{x_,y_,z_}]:=(x+y+z)^4
```

Parameterize the sphere:

```wolfram
s[θ_,φ_]:={Sin[θ]Cos[φ],Sin[θ]Sin[φ],Cos[θ]}
```

Determine infinitesimal surface area:

```wolfram
dS=Simplify[Norm[∂_θs[θ,φ]⨯∂_φs[θ,φ]],0<= θ<=π && -π<φ<π]
(* Output *)
Sin[θ]
```

Perform the integral $\text{DoubleContourIntegral}f \mathrm{d}S$:

```wolfram
∫_0^2π∫_0^πf[s[θ,φ]]dSⅆθⅆφ
(* Output *)
(36 π)/(5)
```

Compare with a region integral:

```wolfram
Subscript[∫, 𝕏∈Sphere[]]f[𝕏]
(* Output *)
(36 π)/(5)
```

Verify Stoke's theorem for $v(x,y,z)=\{-y,x,x^{2}-y^{2}\}$ for the upper unit hemisphere:

```wolfram
v[{x_,y_,z_}]:={-y,x,x^2-y^2}
```

Parameterize the surface using standard spherical coordinates:

```wolfram
s[u_,v_]:= {Sin[u]Cos[v],Sin[u]Sin[v],Cos[u]}
```

Visualize the surface and the vector field:

```wolfram
Show[ParametricPlot3D[s[u,v],{u,0,Pi/2},{v,0,2Pi}],VectorPlot3D[v[{x,y,z}],{x,-1,1},{y,-1,1},{z,0,1}],PlotRange->{{-1,1},{-1,1},{0,1}}]
```

*([Graphics3D])*

The boundary of the surface is the unit circle in the $x y$-plane:

```wolfram
c[t_]:={Cos[t],Sin[t],0}
```

Compute the curl of the vector field:

```wolfram
curl[{x_,y_,z_}]=Curl[v[{x,y,z}],{x,y,z}]
(* Output *)
{-2 y,-2 x,2}
```

Compute the oriented surface area element on the hemisphere:

```wolfram
dS=∂_us[u,v]⨯∂_vs[u,v]//Simplify
(* Output *)
{Cos[v] Sin[u]^2,Sin[u]^2 Sin[v],Cos[u] Sin[u]}
```

Stoke's theorem, $\int \int(\nabla \times \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{S})=\oint \overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{r}$, states that line integral of $\overset{⇀}{v}$ on boundary equals the flux integral of its curl through the surface:

```wolfram
∫_0^(π)/(2)∫_0^2 πcurl[s[u,v]].dSⅆvⅆu==∫_0^2 πv[c[t]].c^′[t]ⅆt
(* Output *)
True
```

Use the divergence theorem to compute the flux of $v(x,y,z)=(x y,e^{x z^{2}}+y^{2},sin(x y))$ through the surface bounded above by $z=1-x^{2}$, below by $z=0$, and on the side by $y+z=2$ and $y=0$:

```wolfram
ContourPlot3D[{y+z==2,z==1-x^2,z==0,y==0},{x,-1.5,1.5},{y,0,2},{z,0,1},ContourStyle->Opacity[0.5],Mesh->None]
```

*([Graphics3D])*

The divergence theorem, $\text{DoubleContourIntegral}\overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{S}=\int \int \int \nabla \cdot \overset{⇀}{v}\mathrm{d}V$, relates the flux to the volume integral of the divergence:

```wolfram
∫_-1^1∫_0^1-x^2∫_0^2-z{x y,y^2+E^(x z^2),Sin[x y]}ⅆyⅆzⅆx
(* Output *)
(184)/(35)
```

Use Gauss's Theorem to find the volume enclosed by the following parametric surface:

```wolfram
s[u_,v_]:={(2+Cos[v])Cos[u],(2+Cos[v])Sin[u],Sin[v](5/4-Cos[u])}
```

```wolfram
ParametricPlot3D[s[u,v],{u,0,2π},{v,0,2π}]
```

*([Graphics3D])*

The oriented area element on the surface is given by the following:

```wolfram
dA=∂_us[u,v]⨯∂_vs[u,v]//FullSimplify
(* Output *)
{-(1)/(4) Cos[u] (-5+4 Cos[u]) Cos[v] (2+Cos[v])+Sin[u]^2 Sin[v]^2,(1)/(4) (5 Cos[v] (2+Cos[v])-4 Cos[u] (1+2 Cos[v])) Sin[u],(2+Cos[v]) Sin[v]}
```

The following vector-field has a divergence equal $1$:

```wolfram
v[{x_,y_,z_}]:=(1)/(3){x,y,z}
```

```wolfram
Div[v[{x,y,z}],{x,y,z}]
(* Output *)
1
```

Apply Gauss's Theorem in the form$\int \int \int 1 \mathrm{d}V=\text{DoubleContourIntegral}\overset{⇀}{v}\cdot \mathrm{d}\overset{⇀}{A}$ to compute the volume:

```wolfram
Integrate[v[s[u,v]].dA,{u,0,2π},{v,0,2π}]
(* Output *)
5 π^2
```

Given a mass density $\rho(x,y,z)=(2-y^{2}) (x^{2}+z^{2})$, find the mass of region given by the following:

```wolfram
φ[u_,v_,w_]:={(5+2w Cos[v])Cos[u],2Sin[v],(5+2w Cos[v])Sin[u]};
```

The ranges of the parameters are $u,v \in[0,2 \pi]$ and $w \in[0,1]$, producing a filled torus:

```wolfram
ParametricPlot3D[φ[u,v,1],{u,0,2π},{v,0,2π}]
```

*([Graphics3D])*

Enter the mass density function:

```wolfram
ρ[{x_,y_,z_}] := (2- y^2)(x^2+z^2)
```

Compute the Jacobian determinate:

```wolfram
jac=Simplify[RealAbs@Det[Grad[φ[u,v,w],{u,v,w}]],0<=w<=1&&v∈Reals]
(* Output *)
4 Cos[v]^2 (5+2 w Cos[v])
```

Integrate to find the total mass:

```wolfram
Integrate[ρ[φ[u,v,w]]jac,{u,0,2π},{v,0,2π},{w,0,1}]
(* Output *)
1160 π^2
```

Derive a formula for the integral of $\sum_{i=1}^{n}x_{i}^{2}$ over an $n$-dimensional unit ball:

```wolfram
Table[Integrate[Sum[Indexed[x,i]^2,{i,n}],x∈Ball[n]],{n,10}]
(* Output *)
{(2)/(3),(π)/(2),(4 π)/(5),(π^2)/(3),(8 π^2)/(21),(π^3)/(8),(16 π^3)/(135),(π^4)/(30),(32 π^4)/(1155),(π^5)/(144)}
```

```wolfram
FindSequenceFunction[%,n]
(* Output *)
(n π^(n/2))/(2 Gamma[2+(n)/(2)])
```

Verify the formula:

```wolfram
Table[(n π^(n/2))/(2 Gamma[2+(n)/(2)]),{n,10}]
(* Output *)
{(2)/(3),(π)/(2),(4 π)/(5),(π^2)/(3),(8 π^2)/(21),(π^3)/(8),(16 π^3)/(135),(π^4)/(30),(32 π^4)/(1155),(π^5)/(144)}
```

#### Average Values and Centroids

Compute the average value of $f(x)=sin(x) cos^{4}(x)$ between $0$ and $\pi$:

```wolfram
(∫_0^πCos[x]^4 Sin[x]ⅆx)/(π-0)
(* Output *)
(2)/(5 π)
```

Visualize the function and its average value:

```wolfram
Plot[{Cos[x]^4Sin[x],2/(5π)},{x,0,π}]
```

*([Graphics])*

Find the mean of $f(x,y)=x+y$ over the parallelogram based at the origin with sides $\{0,1 \}$ and $\{1,1 \}$:

```wolfram
reg=Region[Parallelogram[{0,0},{{1,0},{1,2}}],ImageSize->Small,Axes->True]
```

*([Graphics])*

As $\frac{y}{2}\leq x \leq \frac{y}{2}+1$, the mean is given by the following ratio of integrals:

```wolfram
(Integrate[x+y,{y,0,2},{x,y/2,y/2+1}])/(Integrate[1,{y,0,2},{x,y/2,y/2+1}])
(* Output *)
2
```

Express the integrals using region notation:

```wolfram
(Integrate[x+y,{x,y}∈reg])/(Area[reg])
(* Output *)
2
```

Visualize the function and its mean value:

```wolfram
Plot3D[{x+y,2},{x,y}∈reg,BoxRatios->{1, 1, 1}]
```

*([Graphics3D])*

To compute the centroid of the region under the curve of $f(x)=2 x$ from $0$ to $1$, first find the area:

```wolfram
a=Integrate[2x,{x,0,1}]
(* Output *)
1
```

The centroid equals the average value of the coordinates:

```wolfram
xcentroid=(1/a)Integrate[x,{x,0,1},{y,0,2x}]
(* Output *)
(2)/(3)
```

```wolfram
ycentroid=(1/a)Integrate[y,{x,0,1},{y,0,2x}]
(* Output *)
(2)/(3)
```

Compare with the answer given by [RegionCentroid](https://reference.wolfram.com/language/ref/RegionCentroid.html):

```wolfram
p=RegionCentroid[ImplicitRegion[0<=y<=2x&&0<=x<=1,{x,y}]]
(* Output *)
{(2)/(3),(2)/(3)}
```

```wolfram
Show[Plot[2x,{x,0,1},Filling->Axis],Graphics[{PointSize[Large],Point[p]}]]
```

*([Graphics])*

Determine the centroid of the region between the curves $f(x)=x^{2}$ and $g(x)=\sqrt{x}$ from $0$ to $1$:

```wolfram
a=Integrate[Sqrt[x]-x^2,{x,0,1}]
(* Output *)
(1)/(3)
```

```wolfram
xcentroid=(1/a)Integrate[x(Sqrt[x]-x^2),{x,0,1}]
(* Output *)
(9)/(20)
```

```wolfram
ycentroid=(1/a)Integrate[(1/2)((Sqrt[x])^2-(x^2)^2),{x,0,1}]
(* Output *)
(9)/(20)
```

Compare with the answer returned by [RegionCentroid](https://reference.wolfram.com/language/ref/RegionCentroid.html):

```wolfram
p=RegionCentroid[ImplicitRegion[x^2<=y<=Sqrt[x]&&0<=x<=1,{x,y}]]
(* Output *)
{(9)/(20),(9)/(20)}
```

Visualize the region and its centroid:

```wolfram
Show[Plot[{x^2,Sqrt[x]},{x,0,1},Filling->{1->{2}},PlotLegends->"Expressions"],Graphics[{PointSize[Large],Point[p]}]]
(* Output *)
![image](img/image_036.png)
```

Derive general formulas for the centroid of the region under the curve $y=f(x)\geq 0$ from $a$ to $b$ using the fact that the integral gives the area under the curve:

```wolfram
A=∫_a^bf[x]ⅆx;
```

The $x$ centroid is the mean value of $x$ over the region from $a$ to $b$ and from $0$ to $f(x)$:

```wolfram
(1)/(A)∫_a^b∫_0^f[x]xⅆyⅆx
(* Output *)
(∫_a^bx f[x]ⅆx)/(∫_a^bf[x]ⅆx)
```

The $y$ centroid is similarly the mean value of $y$:

```wolfram
(1)/(A)∫_a^b∫_0^f[x]yⅆyⅆx
(* Output *)
(∫_a^b(f[x]^2)/(2)ⅆx)/(∫_a^bf[x]ⅆx)
```

Find the center of mass of the origin-centered hemisphere of radius $3$ with $y \geq 0$:

```wolfram
reg=Region[ImplicitRegion[x^2+y^2+z^2<=3&&y>=0,{x,y,z}]]
```

*([Graphics3D])*

First compute the volume of the region:

```wolfram
v=Integrate[1,𝕏∈reg]
(* Output *)
2 Sqrt[3] π
```

The center of mass is the average value of the position vector:

```wolfram
cm=(1)/(v)Integrate[𝕏,𝕏∈reg]
(* Output *)
{0,(3 Sqrt[3])/(8),0}
```

Visualize the center of mass:

```wolfram
Show[RegionPlot3D[DiscretizeRegion@reg,Mesh->None,Axes->True,PlotStyle->Opacity[1/2]],Graphics3D[{StandardCyan,Sphere[cm,0.1]}]]
(* Output *)
![image](img/image_038.png)
```

#### Probability, Expectation, and Standard Deviation

Compute the probability that $(x-2)^{2}<1$ when $x$ follows a standard normal distribution:

```wolfram
pdf=PDF[NormalDistribution[0,1],x]
(* Output *)
(ℯ^(-(x^2)/(2)))/(Sqrt[2 π])
```

```wolfram
Integrate[Boole[(x-2)^2<1]pdf,{x,-∞,∞}]
(* Output *)
(1)/(2) (Erfc[(1)/(Sqrt[2])]-Erfc[(3)/(Sqrt[2])])
```

Compare with the value returned by [Probability](https://reference.wolfram.com/language/ref/Probability.html):

```wolfram
Probability[(x-2)^2<1,x∼NormalDistribution[0,1]]
(* Output *)
(1)/(2) (Erfc[-(3)/(Sqrt[2])]-Erfc[-(1)/(Sqrt[2])])
```

Compute the probability that $x>4$ for an exponential distribution with mean $\frac{2}{5}$. The mean is the reciprocal of the parameter:

```wolfram
Mean[ExponentialDistribution[μ]]
(* Output *)
(1)/(μ)
```

Thus, use the probability density function ([PDF](https://reference.wolfram.com/language/ref/PDF.html)) for [ExponentialDistribution](https://reference.wolfram.com/language/ref/ExponentialDistribution.html)[5/2]:

```wolfram
PDF[ExponentialDistribution[(5)/(2)],x]
(* Output *)
{, {{(5)/(2) ℯ^(-5 x/2), x>=0}, {0, True}}}
```

The probability that $x>4$:

```wolfram
Integrate[(5/2)E^(-5t/2),{t,4,∞}]
(* Output *)
(1)/(ℯ^10)
```

Since the PDF vanishes for $x<0$, the probability that $x<2$ is the integral from $0$ to $2$:

```wolfram
Integrate[(5/2)E^(-5t/2),{t,0,2}]
(* Output *)
1-(1)/(ℯ^5)
```

The corresponding probabilistic statements:

```wolfram
Probability[x>4,x∼ExponentialDistribution[5/2]]
(* Output *)
(1)/(ℯ^10)
```

```wolfram
Probability[x<2,x∼ExponentialDistribution[5/2]]//Expand
(* Output *)
1-(1)/(ℯ^5)
```

Compute the probability that a value is within two standard deviations of the mean in a normal distribution:

```wolfram
(1)/(σ Sqrt[2 π])∫_μ-2 σ^μ+2 σExp[-((x-μ)^2)/(2 σ^2)]ⅆx
(* Output *)
Erf[Sqrt[2]]
```

Compare with the answer returned by [Probability](https://reference.wolfram.com/language/ref/Probability.html):

```wolfram
Probability[μ-2σ<=x<= μ+2σ,x∼NormalDistribution[μ,σ]]
(* Output *)
Erf[Sqrt[2]]
```

The value is approximate $95%$:

```wolfram
N[Erf[Sqrt[2]]]
(* Output *)
0.9544997361036416
```

This can be interpreted as saying that about $95%$ of the entire area under the curve lies between $x=-2$ and $x=2$ in the following plot:

```wolfram
p[x_]=PDF[NormalDistribution[0,1],x];
tr[x_]={, {{p[x], -2<=x<=2}, {Indeterminate, True}}};
Plot[{p[x],tr[x]},{x,-3,3},Filling->Axis, AxesLabel->{x/σ,p}]
```

*([Graphics])*

Compute the expectation of $\sqrt{\left|x\right|}$ when $x$ follows a standard Cauchy distribution:

```wolfram
pdf=PDF[CauchyDistribution[0,1],x]
(* Output *)
(1)/(π (1+x^2))
```

```wolfram
Integrate[Sqrt[Abs[x]] pdf,{x,-∞,∞}]
(* Output *)
Sqrt[2]
```

Compare with the answer returned by [Expectation](https://reference.wolfram.com/language/ref/Expectation.html):

```wolfram
Expectation[Sqrt[Abs[x]],x∼CauchyDistribution[0,1]]
(* Output *)
Sqrt[2]
```

Mean and variance of the normal distribution:

```wolfram
Integrate[x PDF[NormalDistribution[μ,σ],x], {x,-∞,∞},Assumptions ->σ>0]
(* Output *)
μ
```

```wolfram
Integrate[(x-μ)^2PDF[NormalDistribution[μ,σ],x], {x,-∞,∞},Assumptions->σ>0]
(* Output *)
σ^2
```

Compare with the built in functions [Mean](https://reference.wolfram.com/language/ref/Mean.html) and [Variance](https://reference.wolfram.com/language/ref/Variance.html):

```wolfram
{Mean[NormalDistribution[μ,σ]],Variance[NormalDistribution[μ,σ]]}
(* Output *)
{μ,σ^2}
```

Show that the standard deviation of an exponential distribution with mean μ is also μ:

```wolfram
Refine[Sqrt[(1)/(μ)∫_0^∞ℯ^(-x /μ)(x-μ)^2ⅆx],μ>0]
(* Output *)
μ
```

Compare with the answers returned by [Mean](https://reference.wolfram.com/language/ref/Mean.html) and [StandardDeviation](https://reference.wolfram.com/language/ref/StandardDeviation.html):

```wolfram
dist=ExponentialDistribution[1/μ];
```

```wolfram
{Mean[dist],StandardDeviation[dist]}
(* Output *)
{μ,μ}
```

Compute the cumulative distribution function ([CDF](https://reference.wolfram.com/language/ref/CDF.html)) from the probability density function ([PDF](https://reference.wolfram.com/language/ref/PDF.html)):

```wolfram
pdf[x_]=PDF[NormalDistribution[0,1],x]
(* Output *)
(ℯ^(-(x^2)/(2)))/(Sqrt[2 π])
```

```wolfram
cdf[x_]=Integrate[pdf[x],{x,-∞,x}]
(* Output *)
(1)/(2) (1+Erf[(x)/(Sqrt[2])])
```

The CDF gives the area under the PDF curve from $-\infty$ to $x$:

```wolfram
Manipulate[Plot[{pdf[x],cdf[x]},{x,-3,u},Filling->{1->Axis},PlotTheme->"Detailed",PlotRange->{{-3,3},{0,1}},PlotLabel->x==u],{{u,0,x},-3.01,3},SaveDefinitions->True]
```

#### Integral Transforms

Compute a Fourier transform:

```wolfram
ft[f_,x_,ω_]:=(1)/(Sqrt[2π])Integrate[f Exp[ⅈ ω x],{x,-∞,∞}]
```

```wolfram
ft[Exp[-x^2],x,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

Compare with [FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html):

```wolfram
FourierTransform[Exp[-x^2],x,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

Find a Laplace transform:

```wolfram
lt[f_,x_,s_]:=Integrate[f Exp[-s x],{x,0,∞}]
```

```wolfram
lt[x^2,x,s]
(* Output *)
(2)/(s^3)
```

Compare with [LaplaceTransform](https://reference.wolfram.com/language/ref/LaplaceTransform.html):

```wolfram
LaplaceTransform[x^2,x,s]
(* Output *)
(2)/(s^3)
```

Define the Hartley transform:

```wolfram
ht[f_,x_,ω_]:=(1)/(Sqrt[2π])Integrate[f(Cos[ω x]+Sin[ω x]),{x,-∞,∞}]
```

```wolfram
ht[Exp[-x^2],x,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

Since the function is even, the Hartley transform is equivalent to [FourierCosTransform](https://reference.wolfram.com/language/ref/FourierCosTransform.html):

```wolfram
FourierCosTransform[Exp[-x^2],x,ω]
(* Output *)
(ℯ^(-(ω^2)/(4)))/(Sqrt[2])
```

Find the Fourier coefficients of a function on `[0,1]`:

```wolfram
a[0]=2Integrate[x^2,{x,0,1}]
(* Output *)
(2)/(3)
```

```wolfram
a[i_]=2Integrate[x^2 Cos[2π i x],{x,0,1},Assumptions->i∈Integers]
(* Output *)
(2 i π Cos[2 i π]+(-1+2 i^2 π^2) Sin[2 i π])/(2 i^3 π^3)
```

```wolfram
b[i_]=2Integrate[x^2 Sin[2π i x],{x,0,1}]
(* Output *)
(-1+(1-2 i^2 π^2) Cos[2 i π]+2 i π Sin[2 i π])/(2 i^3 π^3)
```

Define the partial sums of the transform:

```wolfram
s[n_]:=a[0]/2+∑_i=1^na[i]Cos[2π i x]+∑_i=1^nb[i]Sin[2π i x]
```

Visualize the partial sums, which exhibit the Gibbs phenomenon due to the aperiodicity of $x^{2}$:

```wolfram
Plot[{x^2,s[#]},{x,0,1},PlotTheme->{"Minimal"},PlotLabel->n==#,ImageSize->225]&/@Range[2,8,2]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

Compute a Mellin transform:

```wolfram
mt[f_,x_,s_]:=Integrate[f x^(s-1),{x,0,∞}]
```

```wolfram
mt[Sin[x],x,s]
(* Output *)
Gamma[s] Sin[(π s)/(2)]
```

Compare with [MellinTransform](https://reference.wolfram.com/language/ref/MellinTransform.html):

```wolfram
MellinTransform[Sin[x],x,s,GenerateConditions->True]
(* Output *)
Gamma[s] Sin[(π s)/(2)]
```

Find a Hankel transform:

```wolfram
ht[f_,r_,s_,ν_]:=Integrate[r f BesselJ[ν,r s],{r,0,∞},Assumptions->s>0 &&ν>-1/2]
```

```wolfram
ht[1/r,r,s,0]
(* Output *)
(1)/(s)
```

Compare with [HankelTransform](https://reference.wolfram.com/language/ref/HankelTransform.html):

```wolfram
HankelTransform[1/r,r,s,ν]
(* Output *)
(1)/(s)
```

Compute a quadratic fractional Fourier transform in closed form:

```wolfram
fract[α_, w_]=Sqrt[(1-I Cot[α])/(2π)]E^(I Cot[α] w^2/2) Integrate[E^(-I Csc[α] w t+I Cot[α] t^2/2) UnitBox[t], {t,-∞, ∞}]
(* Output *)
-((((1)/(2)-(ⅈ)/(2)) ℯ^((1)/(2) ⅈ w^2 Cot[α]-ⅈ w^2 Csc[2 α]) Sqrt[1-ⅈ Cot[α]] (Erfi[(((1)/(2)+(ⅈ)/(2)) (w-(Cos[α])/(2)) Sqrt[Csc[α]])/(Sqrt[Cos[α]])]-Erfi[(((1)/(4)+(ⅈ)/(4)) (2 w+Cos[α]) Sqrt[Csc[α]])/(Sqrt[Cos[α]])]))/(Sqrt[2] Sqrt[Cos[α]] Sqrt[Csc[α]]))
```

Visualize the real and imaginary parts of the transform for different values of α:

```wolfram
Table[Plot[{Re[fract[α,w]], Im[fract[α,w]]}, {w,-2,2}, Filling->Axis,PlotRange->All,ImageSize->225, PlotLabel->α==α], {α,{0.02π,0.09π,3.2π,20.9π}}]
(* Output *)
{[Graphics],[Graphics],[Graphics],[Graphics]}
```

```wolfram
s
```

#### Real and Complex Analysis

Define the standard $L^{p}(\mathbb{R})$ norm of a univariate function:

```wolfram
lp[p_,f_,x_]:=Integrate[Abs[f]^p,{x,-∞,∞},Assumptions->p>=1]^(1/p)
```

Also define a formatting for this function:

```wolfram
Format[HoldPattern[lp[p_,f_,x_]]] := HoldForm[Norm[f,L^p]]
```

Compute the norms as a function of $p$ for three different functions:

```wolfram
f:=Exp[-x^2];
g:=Exp[-5 x^2];
h:=HeavisideTheta[x] Exp[-5 x^2];
```

```wolfram
{lp[p,f,x],lp[p,g,x],lp[p,h,x]}
(* Output *)
{p^(-(1)/(2)/p) π^((1)/(2)/p),p^(-(1)/(2)/p) ((π)/(5))^((1)/(2)/p),2^(-1/p) p^(-(1)/(2)/p) ((π)/(5))^((1)/(2)/p)}
```

The norm is always eventually an increasing function of $p$, but it may be initially decreasing:

```wolfram
Plot[{lp[p,f,x],lp[p,g,x],lp[p,h,x]},{p,1,12},PlotTheme->"Detailed",Evaluated->True,PlotRange->All]
```

*([Graphics])*

The Fourier transform is an $L^{2}$ isomorphism (the norm of the function and its transform are equal):

```wolfram
Table[lp[2,fn,x]==lp[2,FourierTransform[fn,x,k],k],{fn,{f,g}}]//FullSimplify
(* Output *)
{True,True}
```

It is not, however, an $L^{p}$ isomorphism for any other value, for example for $p=1$:

```wolfram
Table[lp[1,fn,x]==lp[1,FourierTransform[fn,x,k],k],{fn,{f,g,h}}]
(* Output *)
{False,False,False}
```

Define the weighted inner product for $L_{w}^{2}(a,b)$, with weight $w$ for functions defined on $a<x<b$:

```wolfram
ip[w_,f_,g_,{x_,a_,b_}]:=Integrate[w f Conjugate[g],{x,a,b}]
```

Orthogonality of Legendre polynomials $n$ on $-1<x<1$ with weight function $1$:

```wolfram
Table[ip[1,LegendreP[n,x],LegendreP[m,x],{x,-1,1}],{n,0,4},{m,0,4}]//MatrixForm
(* Output *)
({{2, 0, 0, 0, 0}, {0, (2)/(3), 0, 0, 0}, {0, 0, (2)/(5), 0, 0}, {0, 0, 0, (2)/(7), 0}, {0, 0, 0, 0, (2)/(9)}})
```

Orthogonality of Chebyshev polynomials $T_{n}(x)$ on $-1<x<1$ with weight function $\frac{1}{\sqrt{1-x^{2}}}$:

```wolfram
Table[ip[(1)/(Sqrt[1-x^2]),ChebyshevT[n,x],ChebyshevT[m,x],{x,-1,1}],{n,0,4},{m,0,4}]//MatrixForm
(* Output *)
({{π, 0, 0, 0, 0}, {0, (π)/(2), 0, 0, 0}, {0, 0, (π)/(2), 0, 0}, {0, 0, 0, (π)/(2), 0}, {0, 0, 0, 0, (π)/(2)}})
```

Orthogonality of Hermite polynomials $H_{n}(x)$ on $-\infty<x<\infty$ with weight function $e^{-x^{2}}$:

```wolfram
Table[ip[Exp[-x^2],HermiteH[n,x],HermiteH[m,x],{x,-∞,∞}],{n,0,4},{m,0,4}]//MatrixForm
(* Output *)
({{Sqrt[π], 0, 0, 0, 0}, {0, 2 Sqrt[π], 0, 0, 0}, {0, 0, 8 Sqrt[π], 0, 0}, {0, 0, 0, 48 Sqrt[π], 0}, {0, 0, 0, 0, 384 Sqrt[π]}})
```

Define an inner product on functions using [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
dot[f_,g_] := (2)/(π)∫_-1^1Conjugate[f]gSqrt[1-x^2]ⅆx
```

Construct an orthonormal basis using using [Orthogonalize](https://reference.wolfram.com/language/ref/Orthogonalize.html):

```wolfram
Orthogonalize[{1,x,x^2,x^3,x^4},dot]//Expand
(* Output *)
{1,2 x,-1+4 x^2,-4 x+8 x^3,1-12 x^2+16 x^4}
```

This inner product produces the Gegenbauer polynomials:

```wolfram
Table[GegenbauerC[n,1,x],{n,0,4}]
(* Output *)
{1,2 x,-1+4 x^2,-4 x+8 x^3,1-12 x^2+16 x^4}
```

Compute the residue of $f$ at $z=c$ as an integral over a contour enclosing $c$:

```wolfram
res[f_,{z_,c_}]:=[Limit]_{r->0^(+)}(1)/(2 π ⅈ)∫_0^2 π(f Dt[z,t]/.z->c+Exp[ⅈ t])ⅆt
```

```wolfram
res[(2)/(z-1),{z,1}]
(* Output *)
2
```

```wolfram
res[(Sin[z])/((z-2)^3),{z,2}]
(* Output *)
-(Sin[2])/(2)
```

Compare with the answers returned by [Residue](https://reference.wolfram.com/language/ref/Residue.html):

```wolfram
Residue[(2)/(z-1),{z,1}]
(* Output *)
2
```

```wolfram
Residue[(Sin[z])/((z-2)^3),{z,2}]
(* Output *)
-(Sin[2])/(2)
```

#### Integral Representation of Special Functions

Represent [HermiteH](https://reference.wolfram.com/language/ref/HermiteH.html) in terms of [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Assuming[{Re[z]>0,n∈Integers},FullSimplify[HermiteH[n,z]==(2^n ∫_-∞^∞((z+ⅈ t)^n)/(ℯ^(t^2))ⅆt)/(Sqrt[π])]]
(* Output *)
True
```

Visualize the first five Hermite polynomials:

```wolfram
Plot[Table[HermiteH[n,z],{n,0,4}]//Evaluate,{z,0,4},PlotTheme->{"Detailed","DashedLines"}]
```

*([Graphics])*

Express [Gamma](https://reference.wolfram.com/language/ref/Gamma.html) in terms of a logarithmic integral:

```wolfram
Simplify[Gamma[z]==Integrate[Log[1/t]^(z-1),{t,0,1}],Re[z]>0]
(* Output *)
True
```

Visualize the function:

```wolfram
Plot[Gamma[z],{z,0,5}]
```

*([Graphics])*

Represent [Zeta](https://reference.wolfram.com/language/ref/Zeta.html) in terms of [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
FullSimplify[Zeta[s]==(1)/(Gamma[s])∫_0^∞(t^(s-1))/(ℯ^t-1)ⅆt,Re[s]>1]
(* Output *)
True
```

### Properties & Relations

Integration is a linear operator:

```wolfram
f[x_]=E^x;
```

```wolfram
g[x_]=x^2+4x+17;
```

```wolfram
∫(f[x]+b g[x])ⅆx==∫f[x]ⅆx+b ∫g[x]ⅆx//Simplify
(* Output *)
True
```

Indefinite integration is the inverse of differentiation:

```wolfram
f[x_]:=x^n+x^-n
```

```wolfram
Integrate[f[x],x]
(* Output *)
(x^(1-n))/(1-n)+(x^(1+n))/(1+n)
```

```wolfram
f[x]==D[%,x]
(* Output *)
True
```

Definite integration can be defined in terms of [DiscreteLimit](https://reference.wolfram.com/language/ref/DiscreteLimit.html) and [Sum](https://reference.wolfram.com/language/ref/Sum.html):

```wolfram
f[x_]:=1/x
```

```wolfram
∫_1^3f[x]ⅆx==[Limit]_{n->_{Integers}∞}∑_{k=0}^{n-1}f[1+k((3-1))/(n)](3-1)/(n)
(* Output *)
True
```

Evaluate integrals numerically using [N](https://reference.wolfram.com/language/ref/N.html):

```wolfram
Integrate[Sin[Sin[x]],{x,0,1}]
(* Output *)
∫_0^1Sin[Sin[x]]ⅆx
```

```wolfram
N[%]
(* Output *)
0.43060610312069103
```

This effectively calls [NIntegrate](https://reference.wolfram.com/language/ref/NIntegrate.html):

```wolfram
NIntegrate[Sin[Sin[x]],{x,0,1}]
(* Output *)
0.43060610312069103
```

[Derivative](https://reference.wolfram.com/language/ref/Derivative.html) with a negative integer order does integrals:

```wolfram
Derivative[-2][Function[x,x^n]][x]
(* Output *)
(x^(2+n))/((1+n) (2+n))
```

```wolfram
Integrate[x^n,x,x]//Factor
(* Output *)
(x^(2+n))/((1+n) (2+n))
```

[ArcLength](https://reference.wolfram.com/language/ref/ArcLength.html) is the integral of 1 over a one-dimensional region:

```wolfram
∫_{{x,y}∈Circle[]}1
(* Output *)
2 π
```

```wolfram
ArcLength[Circle[]]
(* Output *)
2 π
```

[Area](https://reference.wolfram.com/language/ref/Area.html) is the integral of 1 over a two-dimensional region:

```wolfram
∫_{{x,y,z}∈Sphere[]}1
(* Output *)
4 π
```

```wolfram
Area[Sphere[]]
(* Output *)
4 π
```

[Volume](https://reference.wolfram.com/language/ref/Volume.html) is the integral of 1 over a three-dimensional region:

```wolfram
∫_{{x,y,z}∈Ball[]}1
(* Output *)
(4 π)/(3)
```

```wolfram
Volume[Ball[]]
(* Output *)
(4 π)/(3)
```

[RegionMeasure](https://reference.wolfram.com/language/ref/RegionMeasure.html) for a region $\mathcal{R}$ is given by the integral $\int_{x \in \mathcal{R}}1$:

```wolfram
ℛ=Circle[];
```

```wolfram
{RegionMeasure[ℛ],Integrate[1,{x,y}∈ℛ]}
(* Output *)
{2 π,2 π}
```

```wolfram
ℛ=Ball[];
```

```wolfram
{RegionMeasure[ℛ],Integrate[1,{x,y,z}∈ℛ]}
(* Output *)
{(4 π)/(3),(4 π)/(3)}
```

[RegionCentroid](https://reference.wolfram.com/language/ref/RegionCentroid.html) is equivalent to [Integrate](https://reference.wolfram.com/language/ref/Integrate.html)[*p*,*p*∈ℛ]/*m* with `*m*=[RegionMeasure](https://reference.wolfram.com/language/ref/RegionMeasure.html)[ℛ]`:

```wolfram
ℛ=Sphere[{1,2,3}];
m=RegionMeasure[ℛ];
```

```wolfram
{RegionCentroid[ℛ],Integrate[p,p∈ℛ]/m}
(* Output *)
{{1,2,3},{1,2,3}}
```

Solve a simple differential equation:

```wolfram
Integrate[Sin[x], x]
(* Output *)
-Cos[x]
```

[DSolveValue](https://reference.wolfram.com/language/ref/DSolveValue.html) returns a solution with the constant of integration:

```wolfram
DSolveValue[y'[x]==Sin[x],y[x], x]
(* Output *)
1-Cos[x]
```

[DSolve](https://reference.wolfram.com/language/ref/DSolve.html) returns a substitution rule for the solution:

```wolfram
DSolve[y'[x]==Sin[x],y, x]
(* Output *)
{{y->Function[{x},1-Cos[x]]}}
```

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html) computes the integral in closed form:

```wolfram
closed=Integrate[Sin[x t]^2,{t,0,1}]
(* Output *)
(1)/(2)-(Sin[2 x])/(4 x)
```

[AsymptoticIntegrate](https://reference.wolfram.com/language/ref/AsymptoticIntegrate.html) gives series approximating the exact result:

```wolfram
AsymptoticIntegrate[Sin[x t]^2,{t,0,1},{x,0,8}]
(* Output *)
(x^2)/(3)-(x^4)/(15)+(2 x^6)/(315)-(x^8)/(2835)
```

```wolfram
Series[closed,{x,0,8}]//Normal
(* Output *)
(x^2)/(3)-(x^4)/(15)+(2 x^6)/(315)-(x^8)/(2835)
```

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html) is defined in terms of an integral:

```wolfram
FourierTransform[Exp[-t^2] Sin[t],t,ω]==∫_-∞^∞Exp[-t^2] Sin[t]( Exp[ⅈ ω t])/(Sqrt[2 π])ⅆt//Simplify
(* Output *)
True
```

[LaplaceTransform](https://reference.wolfram.com/language/ref/LaplaceTransform.html) is defined in terms of an integral:

```wolfram
LaplaceTransform[Exp[-t^2] Sin[t],t,s]==∫_0^∞Exp[-t^2] Sin[t] Exp[-s t]ⅆt//FullSimplify
(* Output *)
True
```

### Possible Issues

#### Indefinite Integrals

Many simple integrals cannot be evaluated in terms of standard mathematical functions:

```wolfram
Integrate[Sin[x]/Log[x], x]
(* Output *)
∫(Sin[x])/(Log[x])ⅆx
```

The indefinite integral of a continuous function can be discontinuous:

```wolfram
f[x_]:=1/(2+Sin[x])
```

```wolfram
Plot[{f[x],Integrate[f[x],x]},{x,0,2π},PlotTheme->"Detailed",Evaluated->True]
```

*([Graphics])*

Using a definite integral with a variable upper limit can smooth the discontinuity. Find the discontinuity to the left of the plot:

```wolfram
Solve[Integrate[f[x],x]==0 && -π<x<0,x]
(* Output *)
{{x->-2 ArcTan[(1)/(2)]}}
```

In this case, integrating from previous discontinuity reproduces the values of the indefinite integral in the range $0 \leq x<\pi$:

```wolfram
Assuming[{x∈Reals},Plot[{f[x],Integrate[f[t],{t,-2 ArcTan[(1)/(2)],x}]},{x,0,2π},PlotTheme->"Detailed",Evaluated->True,Filling->{1->0}]]
(* Output *)
![image](img/image_040.png)
```

The derivative of an integral may not come out in the same form as the original function:

```wolfram
Integrate[1/Sqrt[a^2-x^2],x]
(* Output *)
-ArcTan[(x Sqrt[a^2-x^2])/(-a^2+x^2)]
```

```wolfram
D[%,x]
(* Output *)
-(-(2 x^2 Sqrt[a^2-x^2])/((-a^2+x^2)^2)-(x^2)/(Sqrt[a^2-x^2] (-a^2+x^2))+(Sqrt[a^2-x^2])/(-a^2+x^2))/(1+(x^2 (a^2-x^2))/((-a^2+x^2)^2))
```

[Simplify](https://reference.wolfram.com/language/ref/Simplify.html) and related constructs can often show equivalence:

```wolfram
Simplify[%]
(* Output *)
(1)/(Sqrt[a^2-x^2])
```

Different forms of the same integrand can give integrals that differ by constants of integration:

```wolfram
Integrate[1+(x+1)^3,x]
(* Output *)
x+(1)/(4) (1+x)^4
```

```wolfram
Integrate[Expand[1+(x+1)^3],x]
(* Output *)
2 x+(3 x^2)/(2)+x^3+(x^4)/(4)
```

```wolfram
Simplify[%%-%]
(* Output *)
(1)/(4)
```

Parameters like $n$ are assumed to be generic inside indefinite integrals:

```wolfram
Integrate[x^n,x]
(* Output *)
(x^(1+n))/(1+n)
```

```wolfram
Integrate[x^-1,x]
(* Output *)
Log[x]
```

Use definite integration with a variable upper limit to generate conditions:

```wolfram
Integrate[t^n,{t,0,x}]
(* Output *)
(x^(1+n))/(1+n)
```

When part of a sum cannot be integrated explicitly, the whole sum will stay unintegrated:

```wolfram
Integrate[f[x]+f'[x],x]
(* Output *)
∫(f[x]+f^′[x])ⅆx
```

Compare with:

```wolfram
Integrate[f[x],x]+Integrate[f'[x],x]
(* Output *)
f[x]+∫f[x]ⅆx
```

#### Definite Integrals

Substituting limits into an indefinite integral may not give the correct result for a definite integral:

```wolfram
Integrate[1/(2+Cos[x]),{x,0,2π}]
(* Output *)
(2 π)/(Sqrt[3])
```

```wolfram
expr=Integrate[1/(2+Cos[x]),x]
(* Output *)
(2 ArcTan[(Tan[(x)/(2)])/(Sqrt[3])])/(Sqrt[3])
```

```wolfram
(expr/.x->2π)-(expr/.x->0)
(* Output *)
0
```

The presence of a discontinuity in the expression for the indefinite integral leads to the anomaly:

```wolfram
Plot[expr,{x,0,2π}]
```

*([Graphics])*

Specifying integer assumptions may not give a simpler result:

```wolfram
Integrate[k Cos[k x],{x,-π,π},Assumptions->k∈Integers]
(* Output *)
2 Sin[k π]
```

Use [Simplify](https://reference.wolfram.com/language/ref/Simplify.html) and related functions to obtain the expected result:

```wolfram
Simplify[%,k∈Integers]
(* Output *)
0
```

One can specify a real range as an assumption:

```wolfram
Integrate[Exp[-x] LaguerreL[n,x]^2,{x,0,Infinity},Assumptions->Element[n,PositiveReals]]
(* Output *)
Integrate
(* Output *)
Integrate[ℯ^-x LaguerreL[n,x]^2,{x,0,∞},Assumptions->n∈Reals&&n>0]
```

Numeric integration agrees with this for `n=π`:

```wolfram
NIntegrate[Exp[-x] LaguerreL[Pi,x]^2,{x,0,Infinity}]
(* Output *)
NIntegrate
(* Output *)
1.7363435304473117×10^297
```

It happens that this integral does converge on a subset of zero measure, the positive integers:

```wolfram
Integrate[Exp[-x] LaguerreL[3,x]^2,{x,0,Infinity}]
(* Output *)
1
```

A discrete subcase specification such as integrality will still give the generic result:

```wolfram
Integrate[Exp[-x] LaguerreL[n,x]^2,{x,0,Infinity},Assumptions->Element[n,PositiveIntegers]]
(* Output *)
Integrate
(* Output *)
Integrate[ℯ^-x LaguerreL[n,x]^2,{x,0,∞},Assumptions->n∈Integers&&n>0]
```

A definite integral may have a closed form only over an infinite interval:

```wolfram
Integrate[BesselJ[2,x]/(1+x^2),{x,0,1}]
(* Output *)
∫_0^1(BesselJ[2,x])/(1+x^2)ⅆx
```

```wolfram
Integrate[BesselJ[2,x]/(1+x^2),{x,0,∞}]
(* Output *)
(1)/(6) (2-3 π BesselI[2,1]+3 π StruveL[2,1])
```

Integrals over regions do not test whether an integrand is absolutely integrable:

```wolfram
Integrate[(x^2-y^2)/((x^2+y^2)^2),{x,y}∈Rectangle[{0,0},{1,1}]]
(* Output *)
(π)/(4)
```

```wolfram
[Limit]_{a->0^(+)}Integrate[Abs[(x^2-y^2)/((x^2+y^2)^2)],{x,y}∈Rectangle[{a,a},{1,1}]]
(* Output *)
∞
```

Answers may then depend on how the region was decomposed for integration:

```wolfram
Integrate[(x^2-y^2)/(x^2+y^2)^2,{x,0,1},{y,0,1}]
(* Output *)
(π)/(4)
```

```wolfram
Integrate[(x^2-y^2)/(x^2+y^2)^2,{y,0,1},{x,0,1}]
(* Output *)
-(π)/(4)
```

Integrals over zero-dimensional regions use the counting measure:

```wolfram
ℛ=RegionIntersection[Ball[{0,0,0},Sqrt[14]],Ball[{2,4,6},Sqrt[14]]]
(* Output *)
Point[{1,2,3}]
```

```wolfram
Integrate[x+y+z,{x,y,z}∈ℛ]
(* Output *)
6
```

To use the measure of the ambient space, integrate over all space with the added condition $\{x,y,z \}\in \mathcal{R}$:

```wolfram
Integrate[(x+y+z) Boole[{x,y,z}∈ℛ],{x,y,z}∈FullRegion[3]]
(* Output *)
0
```

Setting [GenerateConditions](https://reference.wolfram.com/language/ref/GenerateConditions.html) to [False](https://reference.wolfram.com/language/ref/False.html) may produce unexpected answers:

```wolfram
Integrate[1/x,{x,0,a},GenerateConditions->False]
(* Output *)
Log[a]
```

In this case, the condition that the integral is divergent was lost:

```wolfram
Integrate[1/x,{x,0,a}]
(* Output *)
Integrate
(* Output *)
∫_0^a(1)/(x)ⅆx
```

### Interactive Examples

Consider Gabriel's horn, the interior of rotating $y=\frac{1}{x}$ around the $x$ axis for $x>1$:

```wolfram
f[x_]:=1/x
```

Compute the volume for arbitrary endpoint $a$:

```wolfram
v[a_]=Integrate[π f[x]^2,{x,1,a},Assumptions->a>1]
(* Output *)
((-1+a) π)/(a)
```

Compute the surface area for arbitrary endpoint $a$:

```wolfram
s[a_]=Integrate[2π f[x]Sqrt[1+D[f[x],x]^2],{x,1,a},Assumptions->a>1]
(* Output *)
π (Sqrt[2]-(Sqrt[1+a^4])/(a^2)-ArcSinh[1]+ArcSinh[a^2])
```

The limit as $a \to \infty$ of the volume is finite, but the surface area is infinite:

```wolfram
[Limit]_{a->∞}{v[a],s[a]}
(* Output *)
{π,∞}
```

Visualize the horn along with its volume and surface area as functions of $a$:

```wolfram
Manipulate[Column[{Style[a==b//,FontSize->14],RevolutionPlot3D[1/s,{s,1,b},RevolutionAxis->{1,0,0},PlotRange->{{1,15},{-1,1},{-1,1}},Mesh->None,ImageSize->300,Ticks->None,Boxed->False],Show[Plot[{s[x],v[x]},{x,1,16},ImageSize->Small,PlotLegends->{"Surface Area"==NumberForm[s[b],{5,3}],"Volume"==NumberForm[v[b],{4,3}]}],Graphics[{StandardOrange,PointSize[Large],Point[{b,v[b]}],StandardBlue,Point[{b,s[b]}]}]]},Alignment->Center],{{b,5.,a},1.01,15},SaveDefinitions->True]
```

### Neat Examples

The first six Borwein-type integrals are all exactly $\frac{\pi}{2}$:

```wolfram
Table[∫_0^∞(∏_{k=0}^{max}(Sin[(x)/(2 k+1)])/((x)/(2 k+1)))ⅆx,{max,6}]
(* Output *)
{(π)/(2),(π)/(2),(π)/(2),(π)/(2),(π)/(2),(π)/(2)}
```

From the seventh onward, they differ from $\frac{\pi}{2}$ by small amounts, for example the eighth:

```wolfram
∫_0^∞(∏_{k=0}^{8}(Sin[(x)/(2 k+1)])/((x)/(2 k+1)))ⅆx
(* Output *)
(17708695183056190642497315530628422295569865119 π)/(35417390788301195294898352987527510935040000000)
```

```wolfram
N[%-Pi/2]
(* Output *)
-1.8724491655091225×10^-8
```

A logarithmic integral from Srinivasa Ramanujan's notebooks:

```wolfram
Integrate[Log[(1+Sqrt[1+4x])/2]/x, {x,0,1}]
(* Output *)
(π^2)/(15)
```

## Tech Notes ▪Symbolic Mathematics: Basic Operations ▪Simplifying Algebraic Expressions ▪Integration ▪Indefinite Integrals ▪Integrals over Regions ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Calculus ▪Functions of Complex Variables ▪Geometric Computation ▪Mathematical Data ▪Solvers over Regions ▪Symbolic Vectors, Matrices and Arrays ▪Precollege Education ▪Integral Transforms ▪Scientific Models ▪Analytic Number Theory ▪Multiplicative Number Theory ▪Polygons ▪Polyhedra

## Related Links [MathWorld](https://mathworld.wolfram.com/Integral.html) [NKS|Online](http://www.wolframscience.com/nks/search/?q=Integrate) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 2003 (5.0) ▪ 2014 (10.0) ▪ 2019 (12.0) ▪ 2026 (15.0)
