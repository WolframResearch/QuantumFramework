# PowerSpectralDensity | [SpanFromLeft]

> [PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html)[*data*,ω] — estimates the power spectral density for `*data*`.
> [PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html)[*data*,ω,*sspec*] — estimates the power spectral density for `*data*` with smoothing specification `*sspec*`.
> [PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html)[*tproc*,ω]  — represents the power spectral density of a time series process `*tproc*`.

## Details and Options

[PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html) is also known as the energy spectral density.

[PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html)[*tproc*,ω] is defined for weakly stationary time series processes as $\sum_{h=-\infty}^{\infty}\mathit{r}(h)e^{-i h \omega}$, where $r(h)$ denotes [CovarianceFunction](https://reference.wolfram.com/language/ref/CovarianceFunction.html)[*proc*,*h*].

The following smoothing specifications `*sspec*` can be given:

*c* | use `*c*` as a cutoff
*w* | use a window function `*w*`
{*c*,*w**}* | use both a cutoff and a window function

For a window function `*w*` and positive integer `*c*`, [PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html)[*data*,ω,{*c*,*w*}] is computed as $\sum_{h=-c}^{c}r(h) w(h/2 c)e^{-i h \omega}$, where $r(h)$ is defined as [CovarianceFunction](https://reference.wolfram.com/language/ref/CovarianceFunction.html)[*data*,*h*].

By default, the cutoff `*c*` is chosen to be $n-1$, where $n$ is the length of `*data*`, and the window function is [DirichletWindow](https://reference.wolfram.com/language/ref/DirichletWindow.html).

A window function $w(x)$ is an even function such that $w(0)=1$, $\left|w(x)\right|\leq 1$, $w(x)=0$ for $\left|x\right|>1/2$, including standard windows such as [HammingWindow](https://reference.wolfram.com/language/ref/HammingWindow.html), [ParzenWindow](https://reference.wolfram.com/language/ref/ParzenWindow.html), etc.

A window function can be given as a list of values `{*w*_0,…}`, where $-1 \leq w_{i}\leq 1$, and it will be applied symmetrically in the vector case.

[PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html) takes the [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) option. Common settings for [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) include:

| `{1,1}` | $\sum_{h=-\infty}^{\infty}r(h) e^{-i h \omega}$ | default setting |
| --- | --- | --- |
| `{-1,1}` | $\frac{1}{2 \pi}\sum_{h=-\infty}^{\infty}r(h) e^{-i h \omega}$ | often used for time series |
| `{*a*,*b*}` | $|\frac{b}{2 \pi}|^{\frac{1-a}{2}}\sum_{h=-\infty}^{\infty}r(h) e^{-i b h \omega}$ | general setting |

## Examples

### Basic Examples

Estimate the power spectral density for some data:

```wolfram
PowerSpectralDensity[Range[10],ω]
(* Output *)
(33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω]-(31)/(10) Cos[6 ω]-(139)/(40) Cos[7 ω]-(63)/(20) Cos[8 ω]-(81)/(40) Cos[9 ω])
```

Calculate the power spectral density for a univariate time series:

```wolfram
PowerSpectralDensity[ARProcess[{a},σ^2],ω]
(* Output *)
(σ^2)/(1+a^2-2 a Cos[ω])
```

The sample power spectral density for a random sample from autoregressive time series:

```wolfram
data=RandomFunction[ARProcess[{.2},.1],{1,1000}];
```

Calculate power spectral density with cutoff:

```wolfram
cutoffs={5,10,17};
spec=Table[PowerSpectralDensity[data,ω,i],{i,cutoffs}];
```

```wolfram
Plot[Evaluate@spec,{ω,-π,π},PlotLegends->cutoffs]
(* Output *)
![image](img/image_001.png)
```

### Scope

#### Empirical Estimates

Estimate the power spectral density for a univariate time series:

```wolfram
data=TimeSeries[...];
```

```wolfram
PowerSpectralDensity[data,ω]
(* Output *)
4.2333486366419555+2 (3.032420687074626 Cos[ω]+1.3528813179050188 Cos[2 ω]-0.2952346003460752 Cos[3 ω]-1.5494784514774493 Cos[4 ω]-1.9409770982927996 Cos[5 ω]-1.4885625928479227 Cos[6 ω]-0.758294483073318 Cos[7 ω]-0.40856823037355844 Cos[8 ω]-0.21245971591442414 Cos[9 ω]+0.1515988490249254 Cos[10 ω])
```

```wolfram
Plot[%,{ω,-π,π},Filling->Axis,PlotRange->All]
```

*([Graphics])*

Power spectral density for a vector time series:

```wolfram
data=TimeSeries[...];
```

```wolfram
psd=PowerSpectralDensity[data,ω]//Simplify
(* Output *)
{{0.6389830242336754-0.1787206515756371 Cos[ω]-0.31653666950425013 Cos[2 ω]-0.14822089061086602 Cos[3 ω]+0.16354821403096959 Cos[4 ω]-0.14855195410236993 Cos[5 ω]+0.227057090491808 Cos[6 ω]-0.11310694308731999 Cos[7 ω]-0.15937794387745835 Cos[8 ω]+0.03168874492672294 Cos[9 ω]+0.0032379790747259145 Cos[10 ω],-0.3813327762627239+0.9753778716944479 Cos[ω]-0.1609699805477443 Cos[2 ω]-0.19050463464711348 Cos[3 ω]+0.22670477139923198 Cos[4 ω]+0.214069412409922 Cos[5 ω]-0.8704202697613537 Cos[6 ω]+0.4552743597499026 Cos[7 ω]-0.027959974919700037 Cos[8 ω]-0.1954684488060646 Cos[9 ω]-0.044770330308804306 Cos[10 ω]-(0.+0.45527501255932096 ⅈ) Sin[ω]+(0.+0.14285376816859008 ⅈ) Cos[ω] Sin[ω]-(0.+0.15050418238038116 ⅈ) Sin[3 ω]-(0.+0.21680581680032238 ⅈ) Sin[4 ω]+(0.+0.12734037754980665 ⅈ) Sin[5 ω]+(0.+0.3022025227943162 ⅈ) Sin[6 ω]-(0.+0.15499325255205612 ⅈ) Sin[7 ω]-(0.+0.30654640104266445 ⅈ) Sin[8 ω]+(0.+0.24582449526526684 ⅈ) Sin[9 ω]+(0.+0.05061435530365383 ⅈ) Sin[10 ω]},{-0.3813327762627239+0.9753778716944479 Cos[ω]-0.1609699805477443 Cos[2 ω]-0.19050463464711348 Cos[3 ω]+0.22670477139923198 Cos[4 ω]+0.214069412409922 Cos[5 ω]-0.8704202697613537 Cos[6 ω]+0.4552743597499026 Cos[7 ω]-0.027959974919700037 Cos[8 ω]-0.1954684488060646 Cos[9 ω]-0.044770330308804306 Cos[10 ω]+(0.+0.45527501255932096 ⅈ) Sin[ω]-(0.+0.14285376816859008 ⅈ) Cos[ω] Sin[ω]+(0.+0.15050418238038116 ⅈ) Sin[3 ω]+(0.+0.21680581680032238 ⅈ) Sin[4 ω]-(0.+0.12734037754980665 ⅈ) Sin[5 ω]-(0.+0.3022025227943162 ⅈ) Sin[6 ω]+(0.+0.15499325255205612 ⅈ) Sin[7 ω]+(0.+0.30654640104266445 ⅈ) Sin[8 ω]-(0.+0.24582449526526684 ⅈ) Sin[9 ω]-(0.+0.05061435530365383 ⅈ) Sin[10 ω],2.24237317095136-1.9336705897477882 Cos[ω]+0.2405091350011714 Cos[2 ω]+1.3574691595651456 Cos[3 ω]-1.0820767525327741 Cos[4 ω]-0.2897249345476386 Cos[5 ω]+0.7726153074270636 Cos[6 ω]-1.3840774143779089 Cos[7 ω]+0.8437943172680488 Cos[8 ω]-0.5950575728142518 Cos[9 ω]-0.17215382619242697 Cos[10 ω]}}
```

Power spectral density for each component:

```wolfram
Table[Plot[psd[[i,i]],{ω,-π,π}],{i,1,2}]
(* Output *)
![image](img/image_003.png)
```

Cross power spectral density between components:

```wolfram
Table[ParametricPlot[{Re[#],Im[#]}&@psd[[i,Mod[i,2]+1]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])]],{i,1,2}]
(* Output *)
![image](img/image_005.png)
```

Estimate the power spectral density for an ensemble of paths:

```wolfram
data=RandomFunction[ARProcess[{.8},1],{0,100},5];
```

```wolfram
psd=PowerSpectralDensity[data,ω];
```

```wolfram
Plot[psd,{ω,-π,π},PlotRange->All,Filling->Axis]
(* Output *)
![image](img/image_007.png)
```

Compare empirical and theoretical power spectral densities functions:

```wolfram
proc=MAProcess[{.4,.3,.5,.6,.3},1.];
```

```wolfram
data=RandomFunction[proc,{0,10^2}];
```

```wolfram
Plot[#,{ω,-π,π},Filling->Axis,PlotRange->All]&/@{PowerSpectralDensity[proc,ω],PowerSpectralDensity[data,ω]}
(* Output *)
![image](img/image_009.png)
```

#### Smoothing

Obtain a smoothed estimate using a cutoff at 5:

```wolfram
data=Range[10];
```

```wolfram
sspec=PowerSpectralDensity[data,ω,5]
(* Output *)
(33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω])
```

Compare the smoothed spectrum to the original:

```wolfram
Plot[{Evaluate@PowerSpectralDensity[data,ω],sspec},{ω,-π,π},Filling->Axis,PlotRange->All]
(* Output *)
![image](img/image_011.png)
```

Compute the power spectral density using a [NuttallWindow](https://reference.wolfram.com/language/ref/NuttallWindow.html):

```wolfram
data=Range[10];
```

```wolfram
sspec=PowerSpectralDensity[data,ω,NuttallWindow]
(* Output *)
(33)/(4)+2 ((231 ((181035)/(2)+121849 Cos[(π)/(9)]+36058 Cos[(2 π)/(9)]) Cos[ω])/(10000000)+(12611277 Cos[3 ω])/(20000000)-(814649 Cos[6 ω])/(5000000)-(63 ((174733)/(2)-121849 Cos[(π)/(9)]+36058 Cos[(2 π)/(9)]) Cos[8 ω])/(5000000)-(17 Cos[5 ω] ((181035)/(2)-36058 Cos[(π)/(9)]-121849 Sin[(π)/(18)]))/(2000000)-(139 Cos[7 ω] ((181035)/(2)-121849 Cos[(2 π)/(9)]+36058 Sin[(π)/(18)]))/(10000000)+(17 Cos[2 ω] ((174733)/(2)+121849 Cos[(2 π)/(9)]+36058 Sin[(π)/(18)]))/(1250000)-(13 Cos[4 ω] ((174733)/(2)-36058 Cos[(π)/(9)]+121849 Sin[(π)/(18)]))/(5000000))
```

Compare the smoothed spectrum to the original:

```wolfram
Plot[{Evaluate@PowerSpectralDensity[data,ω],sspec},{ω,-π,π},Filling->Axis,PlotRange->All]
(* Output *)
![image](img/image_013.png)
```

Define a window using a pure function:

```wolfram
data=Range[10];
```

```wolfram
sspec=PowerSpectralDensity[data,ω,HannPoissonWindow[#,(1)/(2)]&]
(* Output *)
(33)/(4)+2 ((231 (1+Cos[(π)/(9)]) Cos[ω])/(80 ℯ^(1/18))+(17 (1+Cos[(2 π)/(9)]) Cos[2 ω])/(10 ℯ^(1/9))+(147 Cos[3 ω])/(160 ℯ^(1/6))-(31 Cos[6 ω])/(40 ℯ^(1/3))-(139 (1-Cos[(2 π)/(9)]) Cos[7 ω])/(80 ℯ^(7/18))-(63 (1-Cos[(π)/(9)]) Cos[8 ω])/(40 ℯ^(4/9))-(17 Cos[5 ω] (1-Sin[(π)/(18)]))/(16 ℯ^(5/18))-(13 Cos[4 ω] (1+Sin[(π)/(18)]))/(40 ℯ^(2/9)))
```

Compare the smoothed spectrum to the original:

```wolfram
Plot[{Evaluate@PowerSpectralDensity[data,ω],sspec},{ω,-π,π},Filling->Axis,PlotRange->All]
(* Output *)
![image](img/image_015.png)
```

Estimate the power spectral density using specified window function values:

```wolfram
n=10;
data=Range[n];
win=TukeyWindow[(#)/(2*(n-1)),(1)/(3)]&/@Range[0,n-1];
```

Compare to power spectral density with explicit [TukeyWindow](https://reference.wolfram.com/language/ref/TukeyWindow.html):

```wolfram
PowerSpectralDensity[data,ω,win]===PowerSpectralDensity[data,ω,TukeyWindow]
(* Output *)
True
```

Compare the smoothed spectrum to the original:

```wolfram
Plot[{Evaluate@PowerSpectralDensity[data,ω],Evaluate@PowerSpectralDensity[data,ω,TukeyWindow]},{ω,-π,π},Filling->Axis,PlotRange->All]
(* Output *)
![image](img/image_017.png)
```

Compute the power spectral density, given a cutoff and a window function:

```wolfram
data=Range[10];
```

```wolfram
sspec=PowerSpectralDensity[data,ω,{5,BartlettWindow}]
(* Output *)
(33)/(4)+2 ((77 Cos[ω])/(15)+(119)/(45) Cos[2 ω]+(49)/(60) Cos[3 ω]-(13)/(36) Cos[4 ω]-(17)/(18) Cos[5 ω])
```

Compare the smoothed spectrum to the original:

```wolfram
Plot[{Evaluate@PowerSpectralDensity[data,ω],sspec},{ω,-π,π},Filling->Axis,PlotRange->All]
(* Output *)
![image](img/image_019.png)
```

#### Random Processes

Power spectral density for an [ARProcess](https://reference.wolfram.com/language/ref/ARProcess.html):

```wolfram
proc[a_]=ARProcess[{a},1];
```

```wolfram
Plot[PowerSpectralDensity[#,w],{w,-π,π},Filling->Axis,PlotRange->All]&/@{proc[0.8],proc[-0.8]}
(* Output *)
{[Graphics],[Graphics]}
```

```wolfram
PowerSpectralDensity[proc[a],w]
(* Output *)
(1)/(1+a^2-2 a Cos[w])
```

Vector [ARProcess](https://reference.wolfram.com/language/ref/ARProcess.html):

```wolfram
a={{1/7,0},{1/3,0}};
Σ={{1,-1/3},{-1/3,1}};
proc=ARProcess[{a},Σ];
```

```wolfram
psd=PowerSpectralDensity[proc,ω];
psd//MatrixForm
(* Output *)
({{(49)/(50-14 Cos[ω]), -(7 (-7+8 Cos[ω]+8 ⅈ Sin[ω]))/(6 (-25+7 Cos[ω]))}, {(7 (7-8 Cos[ω]+8 ⅈ Sin[ω]))/(6 (-25+7 Cos[ω])), (513-224 Cos[ω])/(450-126 Cos[ω])}})
```

Cross spectral density:

```wolfram
ParametricPlot[{Re[#],Im[#]}&@psd[[1,2]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])]]
```

*([Graphics])*

Power spectral density for an [MAProcess](https://reference.wolfram.com/language/ref/MAProcess.html):

```wolfram
proc[b_]=MAProcess[{b},1];
```

```wolfram
Plot[PowerSpectralDensity[#,w],{w,-π,π},Filling->Axis]&/@{proc[0.9],proc[-0.9]}
(* Output *)
{[Graphics],[Graphics]}
```

```wolfram
PowerSpectralDensity[proc[b],w]
(* Output *)
1+b^2+2 b Cos[w]
```

Vector [MAProcess](https://reference.wolfram.com/language/ref/MAProcess.html):

```wolfram
a={{1/10,0},{1/3,0}};
Σ={{1,1/3},{1/3,1}};
proc=MAProcess[{a},Σ];
```

```wolfram
psd=PowerSpectralDensity[proc,ω];
psd//MatrixForm
(* Output *)
({{(101)/(100)+(Cos[ω])/(5), (1)/(30) (11+11 Cos[ω]+9 ⅈ Sin[ω])}, {(1)/(30) (11+11 Cos[ω]-9 ⅈ Sin[ω]), (2)/(9) (5+Cos[ω])}})
```

Cross spectral density:

```wolfram
ParametricPlot[{Re[#],Im[#]}&@psd[[1,2]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])],AspectRatio->1]
```

*([Graphics])*

Power spectral density for an [ARMAProcess](https://reference.wolfram.com/language/ref/ARMAProcess.html):

```wolfram
proc[a_,b_]=ARMAProcess[{a},{b},1];
```

```wolfram
Plot[PowerSpectralDensity[#,w],{w,-π,π},Filling->Axis,PlotRange->All]&/@{proc[0.3,0.2],proc[-0.3,-0.2]}
(* Output *)
{[Graphics],[Graphics]}
```

```wolfram
PowerSpectralDensity[proc[a,b],w]
(* Output *)
(1+b^2+2 b Cos[w])/(1+a^2-2 a Cos[w])
```

Vector [ARMAProcess](https://reference.wolfram.com/language/ref/ARMAProcess.html):

```wolfram
a={{1/9,0},{1/3,1/2}};
b={{1,3/4},{1/2,-1/4}};
Σ={{1,1/3},{1/3,1}};
proc=ARMAProcess[{a},{b},Σ];
```

```wolfram
psd=PowerSpectralDensity[proc,ω];
psd//FullSimplify//MatrixForm
(* Output *)
({{-(81 (49+40 Cos[ω]))/(32 (-41+9 Cos[ω])), (3 (468+ℯ^(2 ⅈ ω) (588+585 Cos[ω]-265 ⅈ Sin[ω])))/(8 (-9+ℯ^(ⅈ ω)) (-2+ℯ^(ⅈ ω)) (-1+9 ℯ^(ⅈ ω)))}, {-(3 (160+588 ℯ^(ⅈ ω)+425 ℯ^(2 ⅈ ω)+468 ℯ^(3 ⅈ ω)))/(8 (-9+100 ℯ^(ⅈ ω)-173 ℯ^(2 ⅈ ω)+18 ℯ^(3 ⅈ ω))), (7115+712 Cos[ω]+2880 Cos[2 ω])/(5352-5016 Cos[ω]+432 Cos[2 ω])}})
```

Cross spectral density:

```wolfram
ParametricPlot[{Re[#],Im[#]}&@psd[[1,2]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])],AspectRatio->1]
```

*([Graphics])*

Power spectral density for a fractionally integrated time series:

```wolfram
Plot[PowerSpectralDensity[FARIMAProcess[{.3,-.2},-.4,{1},1],w],{w,-π,π},Filling->Axis,PlotRange->All]
```

*([Graphics])*

```wolfram
PowerSpectralDensity[FARIMAProcess[{.3,.2},.4,{1},1],w]
(* Output *)
(2+2 Cos[w])/((2-2 Cos[w])^0.4 (1.1300000000000001-0.48 Cos[w]-0.4 Cos[2 w]))
```

Vector [FARIMAProcess](https://reference.wolfram.com/language/ref/FARIMAProcess.html):

```wolfram
a={{1/9,0},{1/3,1/2}};
b={{1,1/4},{1/2,0}};
Σ={{1,0},{0,1}};
proc=FARIMAProcess[{a},{1/3,-1/4},{b},Σ];
```

```wolfram
psd=PowerSpectralDensity[proc,ω];
psd//FullSimplify//MatrixForm
(* Output *)
({{-(81 (33+32 Cos[ω]))/(32 (1-ℯ^(ⅈ ω))^(1/3) (-41+9 Cos[ω]) (1-Cos[ω]+ⅈ Sin[ω])^(1/3)), (9 (1-ℯ^(ⅈ ω))^(1/4) (36+ℯ^(2 ⅈ ω) (163+156 Cos[ω]-76 ⅈ Sin[ω])))/(8 (1-ℯ^(-ⅈ ω))^(1/3) (-9+ℯ^(ⅈ ω)) (-2+ℯ^(ⅈ ω)) (-1+9 ℯ^(ⅈ ω)))}, {-(9 (1-ℯ^(-ⅈ ω))^(1/4) (40+163 ℯ^(ⅈ ω)+116 ℯ^(2 ⅈ ω)+36 ℯ^(3 ⅈ ω)))/(8 (1-ℯ^(ⅈ ω))^(1/3) (-9+100 ℯ^(ⅈ ω)-173 ℯ^(2 ⅈ ω)+18 ℯ^(3 ⅈ ω))), ((1-ℯ^(ⅈ ω))^(1/4) (2321+288 Cos[ω]+216 Cos[2 ω]) (1-Cos[ω]+ⅈ Sin[ω])^(1/4))/(8 (223-209 Cos[ω]+18 Cos[2 ω]))}})
```

Cross spectral density:

```wolfram
ParametricPlot[{Re[#],Im[#]}&@psd[[1,2]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])],AspectRatio->1]
(* Output *)
![image](img/image_021.png)
```

Power spectral density for a seasonal time series:

```wolfram
Plot[PowerSpectralDensity[SARMAProcess[{.2},{.4},{12,{.3},{1}},1],w],{w,-π,π},Filling->Axis]
(* Output *)
![image](img/image_023.png)
```

```wolfram
PowerSpectralDensity[SARMAProcess[{a},{b},{12,{g},{h}},σ^2],ω]
(* Output *)
(σ^2 (1+b^2+2 b Cos[ω]) (1+h^2+2 h Cos[12 ω]))/((1+a^2-2 a Cos[ω]) (1+g^2-2 g Cos[12 ω]))
```

Vector [SARMAProcess](https://reference.wolfram.com/language/ref/SARMAProcess.html):

```wolfram
a={{1/9,0},{1/3,1/2}};
b={{1/3,1/4},{1/2,0}};
Σ={{1,0},{0,1}};
proc=SARMAProcess[{},{b},{4,{b},{}},Σ];
```

```wolfram
psd=PowerSpectralDensity[proc,ω];
psd//FullSimplify//MatrixForm
(* Output *)
({{(721+384 Cos[ω]+72 Cos[3 ω]+48 Cos[4 ω]+144 Cos[5 ω])/(649-336 Cos[4 ω]-144 Cos[8 ω]), -(6 (8+81 Cos[ω]+8 Cos[3 ω]+81 Cos[4 ω]+21 ⅈ Sin[ω]+8 ⅈ Sin[3 ω]+21 ⅈ Sin[4 ω]))/(-649+336 Cos[4 ω]+144 Cos[8 ω])}, {-(6 (81 Cos[ω]+8 Cos[3 ω]+81 Cos[4 ω]-ⅈ (8 ⅈ+21 Sin[ω]+8 Sin[3 ω]+21 Sin[4 ω])))/(-649+336 Cos[4 ω]+144 Cos[8 ω]), (-937+48 Cos[ω]-288 Cos[3 ω]+384 Cos[4 ω]-144 Cos[5 ω])/(-649+336 Cos[4 ω]+144 Cos[8 ω])}})
```

Cross spectral density:

```wolfram
ParametricPlot[{Re[#],Im[#]}&@psd[[1,2]],{ω,-π,π},ColorFunction->Function[{x,y,ω},(ColorData["Rainbow"][ω])],AspectRatio->1]
(* Output *)
![image](img/image_025.png)
```

### Options

The default value of [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html):

```wolfram
PowerSpectralDensity[Range[10],ω,FourierParameters->{1,1}]
(* Output *)
(33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω]-(31)/(10) Cos[6 ω]-(139)/(40) Cos[7 ω]-(63)/(20) Cos[8 ω]-(81)/(40) Cos[9 ω])
```

```wolfram
PowerSpectralDensity[Range[10],ω]
(* Output *)
(33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω]-(31)/(10) Cos[6 ω]-(139)/(40) Cos[7 ω]-(63)/(20) Cos[8 ω]-(81)/(40) Cos[9 ω])
```

```wolfram
%-%%
(* Output *)
0
```

Change [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html):

```wolfram
PowerSpectralDensity[Range[10],ω,FourierParameters->{-1,1}]
(* Output *)
(1)/(2 π)((33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω]-(31)/(10) Cos[6 ω]-(139)/(40) Cos[7 ω]-(63)/(20) Cos[8 ω]-(81)/(40) Cos[9 ω]))
```

It is the default value scaled:

```wolfram
PowerSpectralDensity[Range[10],ω]
(* Output *)
(33)/(4)+2 ((231 Cos[ω])/(40)+(17)/(5) Cos[2 ω]+(49)/(40) Cos[3 ω]-(13)/(20) Cos[4 ω]-(17)/(8) Cos[5 ω]-(31)/(10) Cos[6 ω]-(139)/(40) Cos[7 ω]-(63)/(20) Cos[8 ω]-(81)/(40) Cos[9 ω])
```

```wolfram
%/%%//Simplify
(* Output *)
2 π
```

### Applications

Use power spectral density for estimating time series processes:

```wolfram
data=RandomFunction[ARProcess[{.3,.2},1],{10^3}];
```

```wolfram
EstimatedProcess[data,ARProcess[2],ProcessEstimator->"SpectralEstimator"]
(* Output *)
ARProcess[0.02202867099602999,{0.3501355911932765,0.19618551999681672},1.0394440247684484]
```

Use a smoothing window:

```wolfram
EstimatedProcess[data,ARProcess[2],ProcessEstimator->{"SpectralEstimator","Window"->5}]
(* Output *)
ARProcess[0.02202867099603011,{0.35013559119327176,0.19618551999681902},1.0394440247684418]
```

### Properties & Relations

Power spectral density of a time series is a transform of the [CovarianceFunction](https://reference.wolfram.com/language/ref/CovarianceFunction.html):

```wolfram
proc=ARProcess[c,{a},σ^2];
cf=CovarianceFunction[proc,h]
(* Output *)
(a^Abs[h] σ^2)/(1-a^2)
```

Use [FourierSequenceTransform](https://reference.wolfram.com/language/ref/FourierSequenceTransform.html):

```wolfram
ft=FourierSequenceTransform[cf,h,z]
(* Output *)
((-1+a^2) ℯ^(ⅈ z) σ^2)/((1-a^2) (-a+ℯ^(ⅈ z)) (-1+a ℯ^(ⅈ z)))
```

Compare to the power spectrum:

```wolfram
PowerSpectralDensity[proc,z]
(* Output *)
(σ^2)/(1+a^2-2 a Cos[z])
```

```wolfram
FullSimplify[%-ft]
(* Output *)
0
```

For a vector time series:

```wolfram
proc=ARProcess[{{{1/3,11/10},{1/5,1/10}}},{{1,0},{0,1}}];
```

```wolfram
cov=CovarianceFunction[proc,h]//PiecewiseExpand;
```

```wolfram
fst=FourierSequenceTransform[cov,h,w]
(* Output *)
{{(450 ℯ^(ⅈ w) (5-111 ℯ^(ⅈ w)+5 ℯ^(2 ⅈ w)))/((15+4 ℯ^(ⅈ w)) (-10+7 ℯ^(ⅈ w)) (-7+10 ℯ^(ⅈ w)) (4+15 ℯ^(ⅈ w))),-(150 ℯ^(ⅈ w) (165-58 ℯ^(ⅈ w)+30 ℯ^(2 ⅈ w)))/((15+4 ℯ^(ⅈ w)) (-10+7 ℯ^(ⅈ w)) (-7+10 ℯ^(ⅈ w)) (4+15 ℯ^(ⅈ w)))},{-(150 ℯ^(ⅈ w) (30-58 ℯ^(ⅈ w)+165 ℯ^(2 ⅈ w)))/((15+4 ℯ^(ⅈ w)) (-10+7 ℯ^(ⅈ w)) (-7+10 ℯ^(ⅈ w)) (4+15 ℯ^(ⅈ w))),(100 ℯ^(ⅈ w) (75-259 ℯ^(ⅈ w)+75 ℯ^(2 ⅈ w)))/(4200+7930 ℯ^(ⅈ w)-27509 ℯ^(2 ⅈ w)+7930 ℯ^(3 ⅈ w)+4200 ℯ^(4 ⅈ w))}}
```

```wolfram
fst-PowerSpectralDensity[proc,w]//Simplify
(* Output *)
{{0,0},{0,0}}
```

Power spectral density of data is a transform of the sample [CovarianceFunction](https://reference.wolfram.com/language/ref/CovarianceFunction.html):

```wolfram
sample=Range[20];
n=Length[sample];
cov=CovarianceFunction[sample,{-n+1,n-1}];
```

Apply [ListFourierSequenceTransform](https://reference.wolfram.com/language/ref/ListFourierSequenceTransform.html):

```wolfram
ListFourierSequenceTransform[cov,w,-n+1]//Simplify
(* Output *)
-(1)/(80) ℯ^(-19 ⅈ w) (19+17 ℯ^(ⅈ w)+15 ℯ^(2 ⅈ w)+13 ℯ^(3 ⅈ w)+11 ℯ^(4 ⅈ w)+9 ℯ^(5 ⅈ w)+7 ℯ^(6 ⅈ w)+5 ℯ^(7 ⅈ w)+3 ℯ^(8 ⅈ w)+ℯ^(9 ⅈ w)-ℯ^(10 ⅈ w)-3 ℯ^(11 ⅈ w)-5 ℯ^(12 ⅈ w)-7 ℯ^(13 ⅈ w)-9 ℯ^(14 ⅈ w)-11 ℯ^(15 ⅈ w)-13 ℯ^(16 ⅈ w)-15 ℯ^(17 ⅈ w)-17 ℯ^(18 ⅈ w)-19 ℯ^(19 ⅈ w))^2
```

Compare to `SamplePowerSpectralDensity`:

```wolfram
PowerSpectralDensity[sample,w]//Simplify
(* Output *)
(1)/(20) (Sin[(w)/(2)]+3 Sin[(3 w)/(2)]+5 Sin[(5 w)/(2)]+7 Sin[(7 w)/(2)]+9 Sin[(9 w)/(2)]+11 Sin[(11 w)/(2)]+13 Sin[(13 w)/(2)]+15 Sin[(15 w)/(2)]+17 Sin[(17 w)/(2)]+19 Sin[(19 w)/(2)])^2
```

```wolfram
%-%%//FullSimplify
(* Output *)
0
```

For a vector values time series:

```wolfram
data=RandomFunction[ARProcess[{{{.3,.1},{.2,.1}}},{{1,0},{0,1}}],{10}];
n=data["PathLength"];
```

```wolfram
cov=CovarianceFunction[data,{-n+1,n-1}]["Values"];
```

```wolfram
lfst=Table[ListFourierSequenceTransform[cov[[All,i,j]],w,-n+1],{i,1,2},{j,1,2}]
(* Output *)
{{0.8918096302806333+0.009390715396253035 ℯ^(-ⅈ w)+0.009390715396253035 ℯ^(ⅈ w)-0.09097683895192663 ℯ^(-2 ⅈ w)-0.09097683895192663 ℯ^(2 ⅈ w)-0.2211332720806821 ℯ^(-3 ⅈ w)-0.2211332720806821 ℯ^(3 ⅈ w)-0.46525080277078507 ℯ^(-4 ⅈ w)-0.46525080277078507 ℯ^(4 ⅈ w)+0.08555076653974467 ℯ^(-5 ⅈ w)+0.08555076653974467 ℯ^(5 ⅈ w)+0.1729281715686234 ℯ^(-6 ⅈ w)+0.1729281715686234 ℯ^(6 ⅈ w)+0.084705665942551 ℯ^(-7 ⅈ w)+0.084705665942551 ℯ^(7 ⅈ w)+0.061273666902055025 ℯ^(-8 ⅈ w)+0.061273666902055025 ℯ^(8 ⅈ w)-0.027926510206162804 ℯ^(-9 ⅈ w)-0.027926510206162804 ℯ^(9 ⅈ w)-0.05446637747998717 ℯ^(-10 ⅈ w)-0.05446637747998717 ℯ^(10 ⅈ w),-0.1188116377296755+0.3690487569121057 ℯ^(-ⅈ w)+0.22301142286650888 ℯ^(ⅈ w)+0.262250622564891 ℯ^(-2 ⅈ w)+0.017052758855269665 ℯ^(2 ⅈ w)-0.497213428191962 ℯ^(-3 ⅈ w)-0.13256527053457548 ℯ^(3 ⅈ w)-0.2988005106386115 ℯ^(-4 ⅈ w)+0.2753457769749431 ℯ^(4 ⅈ w)-0.2421972314099137 ℯ^(-5 ⅈ w)+0.026432793480719266 ℯ^(5 ⅈ w)-0.0145376543394306 ℯ^(-6 ⅈ w)-0.0673815249449312 ℯ^(6 ⅈ w)+0.33468389461497167 ℯ^(-7 ⅈ w)-0.04737180747968221 ℯ^(7 ⅈ w)+0.33556447324595395 ℯ^(-8 ⅈ w)-0.09098928638543984 ℯ^(8 ⅈ w)-0.08596607297528838 ℯ^(-9 ⅈ w)-0.0581186224806773 ℯ^(9 ⅈ w)-0.19751919992661424 ℯ^(-10 ⅈ w)+0.008081747521438653 ℯ^(10 ⅈ w)},{-0.1188116377296755+0.22301142286650888 ℯ^(-ⅈ w)+0.3690487569121057 ℯ^(ⅈ w)+0.017052758855269665 ℯ^(-2 ⅈ w)+0.262250622564891 ℯ^(2 ⅈ w)-0.13256527053457548 ℯ^(-3 ⅈ w)-0.497213428191962 ℯ^(3 ⅈ w)+0.2753457769749431 ℯ^(-4 ⅈ w)-0.2988005106386115 ℯ^(4 ⅈ w)+0.026432793480719266 ℯ^(-5 ⅈ w)-0.2421972314099137 ℯ^(5 ⅈ w)-0.0673815249449312 ℯ^(-6 ⅈ w)-0.0145376543394306 ℯ^(6 ⅈ w)-0.04737180747968221 ℯ^(-7 ⅈ w)+0.33468389461497167 ℯ^(7 ⅈ w)-0.09098928638543984 ℯ^(-8 ⅈ w)+0.33556447324595395 ℯ^(8 ⅈ w)-0.0581186224806773 ℯ^(-9 ⅈ w)-0.08596607297528838 ℯ^(9 ⅈ w)+0.008081747521438653 ℯ^(-10 ⅈ w)-0.19751919992661424 ℯ^(10 ⅈ w),1.339459464095799+0.25018364003858845 ℯ^(-ⅈ w)+0.25018364003858845 ℯ^(ⅈ w)-0.3368822543823049 ℯ^(-2 ⅈ w)-0.3368822543823049 ℯ^(2 ⅈ w)+0.2735307719045323 ℯ^(-3 ⅈ w)+0.2735307719045323 ℯ^(3 ⅈ w)+0.002347697254395113 ℯ^(-4 ⅈ w)+0.002347697254395113 ℯ^(4 ⅈ w)-0.32221943965368105 ℯ^(-5 ⅈ w)-0.32221943965368105 ℯ^(5 ⅈ w)+0.008568346462587521 ℯ^(-6 ⅈ w)+0.008568346462587521 ℯ^(6 ⅈ w)-0.032242431608664475 ℯ^(-7 ⅈ w)-0.032242431608664475 ℯ^(7 ⅈ w)-0.32928881805194193 ℯ^(-8 ⅈ w)-0.32928881805194193 ℯ^(8 ⅈ w)-0.21303523874687838 ℯ^(-9 ⅈ w)-0.21303523874687838 ℯ^(9 ⅈ w)+0.029307994735467694 ℯ^(-10 ⅈ w)+0.029307994735467694 ℯ^(10 ⅈ w)}}
```

```wolfram
PowerSpectralDensity[data,w]-lfst//Simplify//Chop
(* Output *)
{{0,0},{0,0}}
```

Power spectrum of white noise:

```wolfram
FourierSequenceTransform[Piecewise[{{1,x==0}},0],x,w]
(* Output *)
1
```

Compare to special case of an [MAProcess](https://reference.wolfram.com/language/ref/MAProcess.html):

```wolfram
AbsoluteCorrelationFunction[MAProcess[{},1],h]
(* Output *)
{, {{1, Abs[h]==0}, {0, True}}}
```

```wolfram
PowerSpectralDensity[MAProcess[{},1],w]
(* Output *)
1
```

Integrate to find the variance:

```wolfram
Integrate[PowerSpectralDensity[ARProcess[{a},σ^2],w,FourierParameters->{-1,1}],{w,-π,π},Assumptions->σ>0&&-1<=a<=1]
(* Output *)
(σ^2)/(1-a^2)
```

Compare to the variance of the time series:

```wolfram
CovarianceFunction[ARProcess[{a},σ^2],0]
(* Output *)
-(σ^2)/(-1+a^2)
```

```wolfram
%-%%//Simplify
(* Output *)
0
```

Integrate to find the sample second moment:

```wolfram
data=TimeSeries[...];
```

```wolfram
NIntegrate[PowerSpectralDensity[data,w,FourierParameters->{-1,1}],{w,-π,π}]
(* Output *)
0.9464571483850154
```

Compare to the sample second moment:

```wolfram
CovarianceFunction[data,0]
(* Output *)
0.9464571483849886
```

```wolfram
%-%%//Chop
(* Output *)
0
```

Power spectral density for harmonic frequencies is related to [PeriodogramArray](https://reference.wolfram.com/language/ref/PeriodogramArray.html):

```wolfram
data=Range[20];
n=Length[data];
spec=Table[PowerSpectralDensity[data,2π j/n],{j,0,n-1}]//N
(* Output *)
{0.,204.317290945307,52.3606797749979,24.259199981595927,14.47213595499958,10.,7.639320225002102,6.2980809184124915,5.52786404500042,5.125428154684592,5.,5.125428154684592,5.52786404500042,6.2980809184124915,7.639320225002102,10.,14.47213595499958,24.259199981595927,52.3606797749979,204.317290945307}
```

Compare with [PeriodogramArray](https://reference.wolfram.com/language/ref/PeriodogramArray.html):

```wolfram
periodogram=PeriodogramArray[data]
(* Output *)
{2205.,204.317290945307,52.3606797749979,24.25919998159593,14.47213595499957,10.000000000000005,7.639320225002107,6.2980809184124995,5.527864045000413,5.125428154684584,5.000000000000001,5.125428154684584,5.527864045000413,6.298080918412496,7.639320225002105,10.00000000000001,14.47213595499957,24.259199981595923,52.360679774997905,204.317290945307}
```

For zero frequency:

```wolfram
First[periodogram]-n*Mean[data]^2
(* Output *)
0.
```

For nonzero frequencies:

```wolfram
Rest[spec]-Rest[periodogram]//Chop
(* Output *)
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
```

Diagonal elements of the power spectral density for vector data:

```wolfram
data=TimeSeries[...];
```

```wolfram
data["ValueDimensions"]
(* Output *)
2
```

```wolfram
psd=PowerSpectralDensity[data,ω];
```

Compare to univariate power spectral density for each data component:

```wolfram
dens=Table[PowerSpectralDensity[data["Values"][[All,i]],ω],{i,2}];
```

```wolfram
Table[Plot[{psd[[i,i]],dens[[i]]},{ω,-π,π},PlotStyle->{Automatic, Dotted},PlotLegends->{"Vector Component","Univariate"}],{i,2}]
(* Output *)
![image](img/image_027.png)
```

Power spectral density of a vector process is conjugate symmetric about zero:

```wolfram
proc=ARProcess[{{{1/3,1/10},{1/5,1/10}}},{{1,0},{0,1}}];
```

```wolfram
Conjugate/@PowerSpectralDensity[proc,-2]
(* Output *)
{{-(450 (-51+10 Cos[2]))/(26729-19760 Cos[2]+600 Cos[4]),(150 (-8+45 Cos[2]+15 ⅈ Sin[2]))/(26729-19760 Cos[2]+600 Cos[4])},{(150 (-8+45 Cos[2]-15 ⅈ Sin[2]))/(26729-19760 Cos[2]+600 Cos[4]),-(100 (-259+150 Cos[2]))/(26729-19760 Cos[2]+600 Cos[4])}}
```

```wolfram
PowerSpectralDensity[proc,2]
(* Output *)
{{-(450 (-51+10 Cos[2]))/(26729-19760 Cos[2]+600 Cos[4]),(150 (-8+45 Cos[2]+15 ⅈ Sin[2]))/(26729-19760 Cos[2]+600 Cos[4])},{(150 (-8+45 Cos[2]-15 ⅈ Sin[2]))/(26729-19760 Cos[2]+600 Cos[4]),-(100 (-259+150 Cos[2]))/(26729-19760 Cos[2]+600 Cos[4])}}
```

```wolfram
%-%%
(* Output *)
{{0,0},{0,0}}
```

Power spectral density of a univariate process is symmetric about zero:

```wolfram
proc=ARMAProcess[{a},{b},σ^2];
```

```wolfram
PowerSpectralDensity[proc,w]==PowerSpectralDensity[proc,-w]
(* Output *)
True
```

Power spectral density of a vector process is Hermitian:

```wolfram
proc=ARMAProcess[{{{1/3,1/10},{1/5,1/10}}},{{{1/8,-3/10},{3/5,1/6}}},{{1,-.3},{-.3,1}}];
```

```wolfram
HermitianMatrixQ[PowerSpectralDensity[proc,RandomReal[{-π,π}]]]
(* Output *)
True
```

Also non-negative definite:

```wolfram
PositiveSemidefiniteMatrixQ[PowerSpectralDensity[proc,RandomReal[{-π,π}]]]
(* Output *)
True
```

The magnitude of the sample cross spectral density is given by each component:

```wolfram
data=EventSeries[...];
```

```wolfram
data["ValueDimensions"]
(* Output *)
2
```

```wolfram
psd=PowerSpectralDensity[data,3]//Chop
(* Output *)
{{0.1413133945054320042453303578590854855730696622915209735265,-0.1400831942442592704058233955839800042587176702836541704497-0.3643350552757414758324529308343085967847722985298583615981 ⅈ},{-0.1400831942442592704058233955839800042587176702836541704497+0.3643350552757414758324529308343085967847722985298583615981 ⅈ,1.0781945642569339931670159228357343142408265440875456395064}}
```

```wolfram
Abs[psd[[1,2]]]^2-Times@@Diagonal[psd]
(* Output *)
0`23.64791712097895
```

The determinant of the sample power spectral density is constant equal to zero:

```wolfram
FindMaximum[Abs@Det[PowerSpectralDensity[data,w]],{w,-3,3}]
(* Output *)
{9.302259532034058×10^-16,{w->-3.0716186271094075}}
```

Use [TransferFunctionModel](https://reference.wolfram.com/language/ref/TransferFunctionModel.html) to calculate [PowerSpectralDensity](https://reference.wolfram.com/language/ref/PowerSpectralDensity.html) of a time series:

```wolfram
proc=ARMAProcess[{a},{b},σ^2];
psd[ω_]=PowerSpectralDensity[proc,ω]
(* Output *)
(σ^2 (1+b^2+2 b Cos[ω]))/(1+a^2-2 a Cos[ω])
```

Define transfer function:

```wolfram
g[z_]=TransferFunctionModel[proc][z][[1,1]]
(* Output *)
(b+z)/(-a+z)
```

Calculate spectral density:

```wolfram
σ^2g[Exp[I ω]]g[Exp[-I ω]]//ComplexExpand//TrigExpand//Simplify
(* Output *)
(σ^2 (1+b^2+2 b Cos[ω]))/(1+a^2-2 a Cos[ω])
```

Check:

```wolfram
%-psd[ω]//Simplify
(* Output *)
0
```

### Neat Examples

Plot a product of two power spectral densities in 3D:

```wolfram
proc=SARMAProcess[{.2},{.4},{24,{.3},{.1}},1];
Plot3D[Evaluate[PowerSpectralDensity[proc,w]*PowerSpectralDensity[proc,z]],{w,-π,π},{z,-π,π},PlotRange->All,ColorFunction->"IslandColors",ViewPoint->{-3,3,.5}]
(* Output *)
![image](img/image_029.png)
```

## Related Guides ▪Fourier Analysis ▪Time Series Processes

## Related Links ▪ Window Functions ▪ Time Series Processes

## History Introduced in 2012 (9.0)
