# FourierDCT | [SpanFromLeft]

> [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*]  — finds the Fourier discrete cosine transform of a list of real numbers.
> [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*,*m*] — finds the Fourier discrete cosine transform of type `*m*`.

## Details

Possible types `*m*` of discrete cosine transform for a list $u_{r}$ of length $n$ giving a result $v_{s}$ are:

1 (DCT-I) | $\sqrt{\frac{2}{n-1}}(\frac{u_{1}}{2}+\sum_{r=2}^{n-1}u_{r} cos(\frac{\pi}{n-1} (r-1) (s-1))+ (-1)^{s-1}\frac{u_{n}}{2} )$
2 (DCT-II) | $\frac{1}{\sqrt{n}}\sum_{r=1}^{n}u_{r} cos(\frac{\pi}{n} (r-\frac{1}{2})(s-1) )$
3 (DCT-III) | $\frac{1}{\sqrt{n}}(u_{1}+2 \sum_{r=2}^{n}u_{r} cos(\frac{\pi}{n} (r-1)(s-\frac{1}{2}) ))$
4 (DCT-IV) | $\sqrt{\frac{2}{n}}(\sum_{r=1}^{n}u_{r} cos(\frac{\pi}{n} (r-\frac{1}{2}) (s-\frac{1}{2})))$

[FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*] is equivalent to [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*,2].

The inverse discrete cosine transforms for types 1, 2, 3, and 4 are types 1, 3, 2, and 4, respectively.

The `*list*` given in [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*] can be nested to represent an array of data in any number of dimensions.

The array of data must be rectangular.

If the elements of `*list*` are exact numbers, [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html) begins by applying [N](https://reference.wolfram.com/language/ref/N.html) to them.

[FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html) can be used on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) objects.

## Examples

### Basic Examples

Find a discrete cosine transform:

```wolfram
FourierDCT[{0,0,1,0,1}]
(* Output *)
{0.8944271909999159,-0.4253254041760199,-0.08541019662496845,-0.2628655560595667,0.5854101966249684}
```

Find the inverse discrete cosine transform:

```wolfram
FourierDCT[%,3]
(* Output *)
{1.1102230246251565×10^-16,-1.1102230246251565×10^-16,0.9999999999999999,-1.1102230246251565×10^-16,0.9999999999999998}
```

Find a discrete cosine transform of type 1 (DCT-I):

```wolfram
FourierDCT[{1,0,0,1,2},1]
(* Output *)
{1.7677669529663687,-0.8535533905932737,1.0606601717798212,0.14644660940672627,0.35355339059327373}
```

Find the inverse discrete cosine transform:

```wolfram
FourierDCT[%,1]
(* Output *)
{0.9999999999999998,-1.1644331332494318×10^-16,-1.3877787807814457×10^-17,0.9999999999999999,1.9999999999999998}
```

### Scope

Use machine arithmetic to compute the discrete cosine transform:

```wolfram
v = {0,1,2,3,4,3,2,1,0};
```

```wolfram
FourierDCT[v]
(* Output *)
{5.333333333333333,-1.0779341286758226×10^-17,-2.7636197897938626,-3.2051726724102725×10^-17,0.09437286095264925,-1.4990789100506607×10^-16,-0.3333333333333333,-3.312237636943216×10^-16,0.14200734925348715}
```

Use 24-digit precision arithmetic:

```wolfram
FourierDCT[N[v, 24]]
(* Output *)
{5.33333333333333333333333333333333333332,0`22.97063211792824,-2.76361978979386320220707413739488395336,0`22.889215187534685,0.09437286095264951896462088339994358711,0`22.87924626186094,-0.33333333333333333333333333333277631253,0`22.686926201294767,0.14200734925348727882830497944843777253}
```

A two-dimensional discrete cosine transform:

```wolfram
m = RandomReal[1, {3,4}];
```

```wolfram
FourierDCT[m]
(* Output *)
{{1.663592499554336,0.285490095667454,-0.18215854158460085,0.05197600944631247},{0.48088253516152596,0.17495467764446823,-0.03086279893139374,0.12274786033971653},{-0.029090664927172098,0.08104950367167936,-0.28603332357886696,0.08150951474064427}}
```

A five-dimensional discrete cosine transform:

```wolfram
m = RandomReal[1, {5,5,5,5,5}];
```

```wolfram
Max[FourierDCT[m]]
(* Output *)
28.218282567805694
```

### Generalizations & Extensions

The list may have complex values:

```wolfram
FourierDCT[{1,2I,3,4I}]
(* Output *)
{2.+3. ⅈ,-0.1120853822919913-1.4650756326574836 ⅈ,-0.7071067811865476+0.7071067811865476 ⅈ,1.577161014949475-1.6892463972414664 ⅈ}
```

You can use "I", "II", "III", or "IV" for the types 1, 2, 3, and 4, respectively:

```wolfram
FourierDCT[{0,0,1,0,0},"IV"]
(* Output *)
{0.44721359549995804,-0.44721359549995804,-0.44721359549995804,0.44721359549995815,0.4472135954999578}
```

### Applications

#### Compressing Image Data

Import some image data:

```wolfram
ArrayPlot[data=256-Import["ExampleData/ocelot.jpg","Data"]]
```

*([Graphics])*

The two-dimensional DCT:

```wolfram
t=FourierDCT[data];
```

The diagonal spectra shows exponential decay:

```wolfram
ListLogPlot[Abs@Diagonal[t],Joined->True,PlotRange->All]
```

*([Graphics])*

Truncate modes in each axis, effectively compressing by a factor of $f$:

```wolfram
truncate[data_,f_]:=
Module[{i,j},
{i,j}=Floor[Dimensions[data]/Sqrt[f]];
PadRight[Take[data,i,j],Dimensions[data],0.]
];
```

Invert the DCT:

```wolfram
{ArrayPlot[FourierDCT[truncate[t,4],3]],ArrayPlot[FourierDCT[truncate[t,9],3]],ArrayPlot[FourierDCT[truncate[t,16],3]]}
(* Output *)
![image](img/image_001.png)
```

#### Cosine Series Expansion

Get an expansion for an even function as a sum of cosines:

```wolfram
f[x_]:=Exp[-100(x-1/2)^2];
```

The function values on a uniformly spaced grid with $n$ points on $[-L,L)$:

```wolfram
n=25;
xg=N[Range[0,n-1]]/n;
fg=Map[f,xg];
fp=ListPlot[Transpose[{xg,fg}],PlotRange->All]
```

*([Graphics])*

Compute the DCT-III and renormalize:

```wolfram
cc=FourierDCT[fg,3]/Sqrt[n];
```

The function has, in effect, been periodized with a particular symmetry:

```wolfram
Show[fp,Plot[Sum[cc[[r]]*Cos[Pi(r-1/2)x],{r,Length[cc]}],{x,-1,3}, PlotRange->All],PlotRange->All]
(* Output *)
![image](img/image_003.png)
```

Plot the expansion error where the points are defined:

```wolfram
Plot[f[x]-Sum[cc[[r]]*Cos[Pi(r-1/2)x],{r,Length[cc]}],{x,0,1-1/n},PlotRange->All]
```

*([Graphics])*

#### Chebyshev Basis Expansion

Get an expansion for a function in the Chebyshev polynomials:

```wolfram
f[x_]:=1/(1+(5x)^2);
```

The values of the function at the Chebyshev nodes:

```wolfram
n=47;
cnodes=N[Cos[Pi Range[0,n]/n]];
fc=Map[f,cnodes];
pf=ListPlot[Transpose[{cnodes,fc}],PlotRange->All]
```

*([Graphics])*

Find the Chebyshev coefficients:

```wolfram
cc=FourierDCT[fc,1]*Sqrt[2/n];
cc[[{1,-1}]]/=2;
```

Show the error:

```wolfram
Plot[f[x]-Sum[ChebyshevT[i-1,x]*cc[[i]],{i,Length[cc]}],{x,-1,1},PlotRange->All]
```

*([Graphics])*

### Properties & Relations

DCT-I and DCT-IV are their own inverses:

```wolfram
data=RandomReal[1,7];
```

```wolfram
Chop[FourierDCT[FourierDCT[data,1],1]-data]
(* Output *)
{0,0,0,0,0,0,0}
```

```wolfram
Chop[FourierDCT[FourierDCT[data,4],4]-data]
(* Output *)
{0,0,0,0,0,0,0}
```

DCT-II and DCT-III are inverses of each other:

```wolfram
data=RandomReal[1,10];
```

```wolfram
Chop[FourierDCT[FourierDCT[data,2],3]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0,0}
```

```wolfram
Chop[FourierDCT[FourierDCT[data,3],2]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0,0}
```

The DCT is equivalent to matrix multiplication:

```wolfram
dctII[n_]:=(1)/(Sqrt[n])Table[Cos[Pi(r-1/2)(s-1)/n],{s,n},{r,n}]
```

```wolfram
MatrixForm[dctII[7]]
(* Output *)
({{(1)/(Sqrt[7]), (1)/(Sqrt[7]), (1)/(Sqrt[7]), (1)/(Sqrt[7]), (1)/(Sqrt[7]), (1)/(Sqrt[7]), (1)/(Sqrt[7])}, {(Cos[(π)/(14)])/(Sqrt[7]), (Cos[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(7)])/(Sqrt[7]), 0, -(Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(3 π)/(14)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7])}, {(Cos[(π)/(7)])/(Sqrt[7]), (Sin[(π)/(14)])/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), -(1)/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7])}, {(Cos[(3 π)/(14)])/(Sqrt[7]), -(Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7]), 0, (Cos[(π)/(14)])/(Sqrt[7]), (Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(3 π)/(14)])/(Sqrt[7])}, {(Sin[(3 π)/(14)])/(Sqrt[7]), -(Cos[(π)/(7)])/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), (1)/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), -(Cos[(π)/(7)])/(Sqrt[7]), (Sin[(3 π)/(14)])/(Sqrt[7])}, {(Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7]), (Cos[(3 π)/(14)])/(Sqrt[7]), 0, -(Cos[(3 π)/(14)])/(Sqrt[7]), (Cos[(π)/(14)])/(Sqrt[7]), -(Sin[(π)/(7)])/(Sqrt[7])}, {(Sin[(π)/(14)])/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), -(1)/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(14)])/(Sqrt[7])}})
```

```wolfram
data=RandomReal[1,7];
```

```wolfram
Chop[dctII[7].data-FourierDCT[data]]
(* Output *)
{0,0,0,0,0,0,0}
```

### Possible Issues

[FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html) always returns normalized results:

```wolfram
TableForm[Table[FourierDCT[ConstantArray[1,5],type],{type,1,4}]]
(* Output *)
{{2.82842712474619, 0., 0., 0., 0.}, {2.23606797749979, 9.444121133484361×10^-17, 8.033649276373055×10^-17, 5.836787854364518×10^-17, 3.06858096987851×10^-17}, {2.8235955159711317, -0.8777061007329484, 0.44721359549995776, -0.22786670826713562, 0.07083167502878462}, {2.021471201601977, -0.6965515053690705, 0.4472135954999583, -0.3549107188691969, 0.32016958489789743}}
```

To get unnormalized results, you can multiply by the normalization:

```wolfram
nc[n_,1]=Sqrt[(n-1)/2];
nc[n_,2|3]=Sqrt[n];
nc[n_,4]=Sqrt[n/2];
```

```wolfram
unnormalizedDCT[data_,type_]:=FourierDCT[data,type]*nc[Length[data],type]
```

```wolfram
TableForm[Table[unnormalizedDCT[ConstantArray[1,5],type],{type,1,4}]]
(* Output *)
{{4., 0., 0., 0., 0.}, {5.000000000000001, 2.1117696842213397×10^-16, 1.7963785889362148×10^-16, 1.3051454412604205×10^-16, 6.861555643110583×10^-17}, {6.313751514675044, -1.9626105055051508, 0.9999999999999997, -0.5095254494944286, 0.1583844403245368}, {3.196226610749831, -1.1013446322926335, 0.7071067811865481, -0.5611631188171807, 0.5062325628940022}}
```

## Tech Notes ▪Discrete Fourier Transforms

## Related Guides ▪Fourier Analysis ▪Image Filtering & Neighborhood Processing ▪Summation Transforms

## History Introduced in 2007 (6.0)
