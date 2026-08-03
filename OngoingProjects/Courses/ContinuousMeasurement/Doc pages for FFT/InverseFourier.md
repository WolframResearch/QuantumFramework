# InverseFourier | [SpanFromLeft]

> [InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html)[*list*] — finds the discrete inverse Fourier transform of a list of complex numbers.
> [InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html)[*list*,{*p*_1,*p*_2,…}] — returns the specified positions of the discrete inverse Fourier transform.

## Details and Options

The inverse Fourier transform $u_{r}$ of a list $v_{s}$ of length $n$ is defined to be $\frac{1}{\sqrt{n}}\sum_{s=1}^{n}v_{s}e^{-2 \pi i(r-1)(s-1)/n}$.

Note that the zero frequency term must appear at position 1 in the input list.

Other definitions are used in some scientific and technical fields.

Different choices of definitions can be specified using the option [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html).

With the setting [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html)->{*a*,*b*} the discrete Fourier transform computed by [InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html) is $\frac{1}{n^{(1+a)/2}}\sum_{s=1}^{n}v_{s}e^{-2 \pi i b(r-1)(s-1)/n}$.

Some common choices for `{*a*,*b*}` are `{0,1}` (default), `{-1,1}` (data analysis), `{1,-1}` (signal processing).

The setting `*b*=-1 effectively corresponds to conjugating both input and output lists.

To ensure a unique discrete Fourier transform, [Abs](https://reference.wolfram.com/language/ref/Abs.html)[*b*] must be relatively prime to $n$.

The list of data need not have a length equal to a power of two.

The `*list*` given in [InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html)[*list*] can be nested to represent an array of data in any number of dimensions.

The array of data must be rectangular.

[InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html)[*list*,{*p*_1,*p*_2,…}] is typically equivalent to [Extract](https://reference.wolfram.com/language/ref/Extract.html)[InverseFourier[*list**]*,{*p*_1,*p*_2,…}]. Cases with just a few positions `*p*` are computed using an algorithm that takes less time and memory but is more subject to numerical error, particularly when the length of `*list*` is long.

If the elements of `*list*` are exact numbers, [InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html) begins by applying [N](https://reference.wolfram.com/language/ref/N.html) to them.

## Examples

### Basic Examples

Inverse Fourier transform of a real list:

```wolfram
InverseFourier[{1,0,1,1,0,1,1,1,0,1,1,1,1}]
(* Output *)
{2.773500981126146+0. ⅈ,0.0603679140667378+0.20430072968345037 ⅈ,0.016615856108415265+0.3196654252552805 ⅈ,-0.43656498892450984-0.08181776968389065 ⅈ,0.33420964844223305+0.6010289732440861 ⅈ,0.23130902217903387+0.24643392871528733 ⅈ,0.2100876952970121-0.4212072796263528 ⅈ,0.2100876952970121+0.4212072796263528 ⅈ,0.23130902217903387-0.24643392871528733 ⅈ,0.33420964844223305-0.6010289732440861 ⅈ,-0.43656498892450984+0.08181776968389065 ⅈ,0.016615856108415265-0.3196654252552805 ⅈ,0.0603679140667378-0.20430072968345037 ⅈ}
```

Inverse Fourier transform of a complex list:

```wolfram
InverseFourier[{1+2I,3+4I}]
(* Output *)
{2.82842712474619+4.242640687119285 ⅈ,-1.414213562373095-1.414213562373095 ⅈ}
```

### Scope

`*x*` is a list of real values:

```wolfram
x={0,0,0,1,0,0,0};
```

Compute the inverse Fourier transform with machine arithmetic:

```wolfram
InverseFourier[x]
(* Output *)
{0.3779644730092272+0. ⅈ,-0.3405342233544578-0.1639926388028409 ⅈ,0.23565699438616372+0.29550452425305174 ⅈ,-0.0841050075363194-0.3684881145497891 ⅈ,-0.0841050075363194+0.3684881145497891 ⅈ,0.23565699438616372-0.29550452425305174 ⅈ,-0.3405342233544578+0.1639926388028409 ⅈ}
```

Compute using 24-digit precision arithmetic:

```wolfram
InverseFourier[N[x,24]]
(* Output *)
{0.3779644730092272272145165362341800608157513118689214543384,-0.3405342233544579049302902426908887412965089041467792502342-0.1639926388028408838980160125022749532043526605606659087136 ⅈ,0.2356569943861637213799973131958250370646165677548010403279+0.2955045242530517627289728589996899042235694097750702879501 ⅈ,-0.0841050075363194300569653386220263261759833195424825172038-0.3684881145497891211690431535025850489807832507855956206812 ⅈ,-0.0841050075363194300569653386220263261759833195424825172038+0.3684881145497891211690431535025850489807832507855956206812 ⅈ,0.2356569943861637213799973131958250370646165677548010403279-0.2955045242530517627289728589996899042235694097750702879501 ⅈ,-0.3405342233544579049302902426908887412965089041467792502342+0.1639926388028408838980160125022749532043526605606659087136 ⅈ}
```

Compute a 2D inverse Fourier transform:

```wolfram
InverseFourier[RandomComplex[1+I,{4,5}]]
(* Output *)
{{2.4373544008954253+1.5944579818778712 ⅈ,0.30216833793951153+0.11735100210594394 ⅈ,0.40146138154230426-0.19855388277580868 ⅈ,-0.023893301398035912+0.025676837851804062 ⅈ,-0.05629173396766337-0.09993379323987886 ⅈ},{-0.20715518081966505+0.3344504771738099 ⅈ,0.03222343805516017-0.30542096169425986 ⅈ,-0.5153618503742753-0.07279009649057855 ⅈ,0.20350890377196418-0.17769111321888792 ⅈ,0.37095287378514785-0.19226056532005115 ⅈ},{0.08170862176577313+0.32083110509019797 ⅈ,0.05147656665336799-0.6553513945951109 ⅈ,-0.3098428516349509+0.44891213543375613 ⅈ,-0.07459469734530826+0.18381120421061184 ⅈ,0.7395195134801859+0.12140153930003476 ⅈ},{0.04455268168131593+0.05848245308142783 ⅈ,0.006672541955659615-0.020950397845291072 ⅈ,-0.24936375299285726+0.24117074304579733 ⅈ,0.06382333518248073-0.006839518154346554 ⅈ,0.3985375752694679-0.16486803507804645 ⅈ}}
```

`*x*` is a rank-4 tensor with a single nonzero entry:

```wolfram
x=ConstantArray[0,{2,3,3,2}];x[[1,1,1,1]]=1;
```

Compute the 4D inverse Fourier transform:

```wolfram
InverseFourier[x]
(* Output *)
{{{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}},{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}},{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}}},{{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}},{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}},{{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666},{0.16666666666666666,0.16666666666666666}}}}
```

### Options

#### FourierParameters

No normalization:

```wolfram
a=InverseFourier[{1,0,1,0,0,1,0,0,0,1},FourierParameters->{1,1}]
(* Output *)
{0.4+0. ⅈ,0.1118033988749895-0.03632712640026803 ⅈ,0.15000000000000002+0.03632712640026803 ⅈ,-0.1118033988749895+0.15388417685876268 ⅈ,0.15000000000000002+0.15388417685876268 ⅈ,0.+0. ⅈ,0.15000000000000002-0.15388417685876268 ⅈ,-0.1118033988749895-0.15388417685876268 ⅈ,0.15000000000000002-0.03632712640026803 ⅈ,0.1118033988749895+0.03632712640026803 ⅈ}
```

Normalization by $1/\sqrt{n}$:

```wolfram
InverseFourier[{1,0,1,0,0,1,0,0,0,1}]
(* Output *)
{1.2649110640673518+0. ⅈ,0.3535533905932738-0.11487646027368055 ⅈ,0.4743416490252569+0.11487646027368055 ⅈ,-0.3535533905932738+0.4866244947338651 ⅈ,0.4743416490252569+0.4866244947338651 ⅈ,0.+0. ⅈ,0.4743416490252569-0.4866244947338651 ⅈ,-0.3535533905932738-0.4866244947338651 ⅈ,0.4743416490252569-0.11487646027368055 ⅈ,0.3535533905932738+0.11487646027368055 ⅈ}
```

```wolfram
%*Sqrt[10]
(* Output *)
{4.+0. ⅈ,1.118033988749895-0.36327126400268034 ⅈ,1.5+0.36327126400268034 ⅈ,-1.118033988749895+1.5388417685876268 ⅈ,1.5+1.5388417685876268 ⅈ,0.+0. ⅈ,1.5-1.5388417685876268 ⅈ,-1.118033988749895-1.5388417685876268 ⅈ,1.5-0.36327126400268034 ⅈ,1.118033988749895+0.36327126400268034 ⅈ}
```

Normalization by $1/n$:

```wolfram
InverseFourier[{1,0,1,0,0,1,0,0,0,1},FourierParameters->{-1,1}]
(* Output *)
{4.+0. ⅈ,1.118033988749895-0.3632712640026803 ⅈ,1.5+0.3632712640026803 ⅈ,-1.118033988749895+1.5388417685876268 ⅈ,1.5+1.5388417685876268 ⅈ,0.+0. ⅈ,1.5-1.5388417685876268 ⅈ,-1.118033988749895-1.5388417685876268 ⅈ,1.5-0.3632712640026803 ⅈ,1.118033988749895+0.3632712640026803 ⅈ}
```

[InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html) is the same as [Fourier](https://reference.wolfram.com/language/ref/Fourier.html) with parameter $b=-1$:

```wolfram
v=RandomVariate[NormalDistribution[],{8,2}].{1,ⅈ}
(* Output *)
{-2.1560769774957254-0.42195617231169325 ⅈ,0.3868236322087926-0.1333710020375242 ⅈ,-1.6616604651821527+1.5742189399538062 ⅈ,0.8192593972588775-0.7085064680337076 ⅈ,-0.21319883410116655+0.9763601433534722 ⅈ,1.2013847967970253-0.5104104566592209 ⅈ,0.2904482043494613-1.3178023929173275 ⅈ,-0.5532019465972225+0.1158714158656435 ⅈ}
```

```wolfram
InverseFourier[v]
(* Output *)
{-0.6668802516633237-0.15047090627259577 ⅈ,-0.3230174416325396+0.35667445420740823 ⅈ,-0.37095173784512225-0.3620964971585634 ⅈ,-1.2744740832197095-1.6243834928919691 ⅈ,-1.9780442292990323+0.7238075929309215 ⅈ,0.9941630272284538+0.03491587458081778 ⅈ,-0.334785766955368+0.5728054253810222 ⅈ,-2.1443161228030485-0.7447247339975115 ⅈ}
```

```wolfram
Fourier[v,FourierParameters->{0,-1}]
(* Output *)
{-0.6668802516633237-0.15047090627259577 ⅈ,-0.3230174416325396+0.35667445420740823 ⅈ,-0.37095173784512225-0.3620964971585634 ⅈ,-1.2744740832197095-1.6243834928919691 ⅈ,-1.9780442292990323+0.7238075929309215 ⅈ,0.9941630272284538+0.03491587458081778 ⅈ,-0.334785766955368+0.5728054253810222 ⅈ,-2.1443161228030485-0.7447247339975115 ⅈ}
```

Data from a sinc function with noise:

```wolfram
x=Table[Sinc[x-10],{x,100}]+RandomReal[{-.05,.05},{100}];
```

```wolfram
ListPlot[x,PlotRange->All]
```

*([Graphics])*

Get the Fourier transform:

```wolfram
f=Fourier[x,FourierParameters->{-1,1}];
```

```wolfram
ListPlot[Abs[f]]
```

*([Graphics])*

Reconstruct the signal from part of the spectrum:

```wolfram
s=2 InverseFourier[Take[f,20],FourierParameters->{-1,.2}];
(* Output *)
InverseFourier
```

```wolfram
Show[ListPlot[Re[s],PlotStyle->Red],ListPlot[x]]
```

*([Graphics])*

### Applications

Some Gaussian data:

```wolfram
data=Table[Exp[-((x-50)/10)^2],{x,100}];
ListPlot[data,PlotRange->All]
```

*([Graphics])*

The multiplication of each mode to get the first derivative:

```wolfram
k=-I*(2*Pi/100)*Join[{0},Range[49],{0},-Reverse[Range[49]]];
ListPlot[Abs[k]]
```

*([Graphics])*

Approximate the first derivative of the data:

```wolfram
d=InverseFourier[Fourier[data]*k];
```

```wolfram
Show[ListPlot[Re[d],PlotRange->All],Plot[Evaluate[D[Exp[-((x-50)/10)^2],x]],{x,0,100}]]
```

*([Graphics])*

Note the derivative approximation implicitly assumes periodicity:

```wolfram
ListLinePlot[Re[InverseFourier[Fourier[Range[100]]*k]]]
```

*([Graphics])*

### Properties & Relations

$v_{s}$ is given by $\frac{1}{\sqrt{\mathit{n}}}\sum_{\mathit{r}=1}^{\mathit{n}}\mathit{u}_{\mathit{r}}\mathit{e}^{-2 \pi \mathit{i}(\mathit{r}-1)(\mathit{s}-1)/\mathit{n}}$:

```wolfram
x=N[{1,2,3,4,5,6,7}];
n=Length[x];
u/:u_r_:=x[[r]]
```

```wolfram
v=Table[(1)/(Sqrt[n])∑_{r=1}^{n}u_rℯ^(-2π ⅈ(r-1)(s-1)/n),{s,1,n}]
(* Output *)
{10.583005244258361,-1.3228756555322945+2.746979603717467 ⅈ,-1.3228756555322945+1.054958132087371 ⅈ,-1.3228756555322945+0.3019377358048386 ⅈ,-1.3228756555322945-0.3019377358048386 ⅈ,-1.3228756555322945-1.054958132087371 ⅈ,-1.3228756555322945-2.746979603717467 ⅈ}
```

```wolfram
Chop[v-InverseFourier[x]]
(* Output *)
{0,0,0,0,0,0,0}
```

[InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html) is equivalent to multiplication with [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) with specified parameters:

```wolfram
m[n_]:=FourierMatrix[n,FourierParameters->{0,-1}]
```

```wolfram
v=m[6].N[{1,2,3,4,5,6}]
(* Output *)
{8.573214099741124+0. ⅈ,-1.2247448713915883+2.1213203435596433 ⅈ,-1.224744871391588+0.7071067811865472 ⅈ,-1.2247448713915894+0. ⅈ,-1.224744871391588-0.7071067811865472 ⅈ,-1.2247448713915883-2.1213203435596433 ⅈ}
```

```wolfram
Chop[v-InverseFourier[{1,2,3,4,5,6}]]
(* Output *)
{0,0,0,0,0,0}
```

The conjugate transpose of the matrix is equivalent to [Fourier](https://reference.wolfram.com/language/ref/Fourier.html):

```wolfram
Chop[ConjugateTranspose[m[6]].N[{1,2,3,4,5,6}]-Fourier[{1,2,3,4,5,6}]]
(* Output *)
{0,0,0,0,0,0}
```

### Possible Issues

[InverseFourier](https://reference.wolfram.com/language/ref/InverseFourier.html) uses an efficient algorithm when only a small number of coefficients is needed:

```wolfram
bdata=RandomReal[1,2^20,WorkingPrecision->50];data=N[bdata];pos={{1},{3},{2^20-5}};
AbsoluteTiming[res=Extract[InverseFourier[bdata],pos]]
(* Output *)
{17.328942,{512.05267305662334734066844336847530935455442690673350319446805671295000126475233,-0.18383055759370937607861070514656063554301082603087081642057297007146349990304+0.38104314919085593014909377337612901880306323063122975541084107848815619989935 ⅈ,-0.21919036395660422801475519913156810176094615918809252251784006420948796575385-0.19436914286018982082259108345185086915651328565645052275338589595580593068222 ⅈ}}
```

The fast and efficient implementation may result in significant numerical error:

```wolfram
AbsoluteTiming[Extract[InverseFourier[data],pos]-res]
(* Output *)
{0.028899,{0.-3.760334171782637×10^-17 ⅈ,5.551115123125783×10^-17+1.1102230246251565×10^-16 ⅈ,2.7755575615628914×10^-17-1.3877787807814457×10^-16 ⅈ}}
```

```wolfram
AbsoluteTiming[InverseFourier[data,pos]-res]
(* Output *)
{0.019503,{-0.00001107146442791418+0. ⅈ,-0.00015599512077724143+9.991203808734639×10^-9 ⅈ,0.000019868207109363656+7.346321317935889×10^-8 ⅈ}}
```

## Tech Notes ▪Manipulating Numerical Data ▪Discrete Fourier Transforms

## Related Guides ▪Summation Transforms ▪Fourier Analysis ▪GPU Computing ▪Data Transforms and Smoothing ▪GPU Computing with NVIDIA ▪GPU Computing with Apple ▪Signal Transforms ▪Signal Visualization & Analysis

## History Introduced in 1988 (1.0) | Updated in 1996 (3.0) ▪ 1999 (4.0) ▪ 2000 (4.1) ▪ 2002 (4.2) ▪ 2012 (9.0)
