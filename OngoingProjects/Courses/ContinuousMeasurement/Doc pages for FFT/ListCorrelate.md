# ListCorrelate | [SpanFromLeft]

> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*] — forms the correlation of the kernel `*ker*` with `*list*`.
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,*k*] — forms the cyclic correlation in which the `*k*`$^{th}$ element of `*ker*` is aligned with each element in `*list*`.
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,{*k*_*L*,*k*_*R*}] — forms the cyclic correlation whose first element contains `*list*[[1]]*ker*[[*k*_*L*]]` and whose last element contains `*list*[[-1]]*ker*[[*k*_*R*]]`.
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,*klist*,*p*] — forms the correlation in which `*list*` is padded at each end with repetitions of the element `*p*`.
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,*klist*,{*p*_1,*p*_2,…}] — forms the correlation in which `*list*` is padded at each end with cyclic repetitions of the `*p*_*i*`.
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,*klist*,*padding*,*g*,*h*] — forms a generalized correlation in which `*g*` is used in place of [Times](https://reference.wolfram.com/language/ref/Times.html) and `*h*` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html).
> [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*,*klist*,*padding*,*g*,*h*,*lev*] — forms a correlation using elements at level `*lev*` in `*ker*` and `*list*`.

## Details

With kernel `*K*_*r*` and list `*a*_*s*`, [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*] computes $\sum_{r}K_{r}a_{s+r}$, where the limits of the sum are such that the kernel never overhangs either end of the list.

For a one[Hyphen]dimensional list [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html)[*ker*,*list*] is equivalent to [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[[Reverse](https://reference.wolfram.com/language/ref/Reverse.html)[*ker*],*list*].

For higher-dimensional lists, `*ker*` must be reversed at every level.

Settings for `*k*_*L*` and `*k*_*R*` are negated in [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) relative to [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html).

Common settings for `{*k*_*L*,*k*_*R*}` in [ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html) are:

{1,-1} | no overhangs (default)
{1,1} | maximal overhang at the right[Hyphen]hand end
{-1,-1} | maximal overhang at the left[Hyphen]hand end
{-1,1} | maximal overhangs at both beginning and end

## Examples

### Basic Examples

Correlate a kernel `{x,y}` with a list of data:

```wolfram
ListCorrelate[{x,y},{a,b,c,d,e,f}]
(* Output *)
{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y}
```

Make a cyclic correlation the same length as the original data:

```wolfram
ListCorrelate[{x,y},{a,b,c,d,e,f},1]
(* Output *)
{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y,f x+a y}
```

Align element 2 in the kernel with successive elements in the data:

```wolfram
ListCorrelate[{x,y},{a,b,c,d,e,f},2]
(* Output *)
{f x+a y,a x+b y,b x+c y,c x+d y,d x+e y,e x+f y}
```

Pad with `zzz` instead of using the data cyclically:

```wolfram
ListCorrelate[{x,y},{a,b,c,d,e,f},1,zzz]
(* Output *)
{a x+b y,b x+c y,c x+d y,d x+e y,e x+f y,f x+y zzz}
```

Two-dimensional correlation:

```wolfram
ListCorrelate[{{1,1},{1,1}},{{a,b,c},{d,e,f},{g,h,i}}]
(* Output *)
{{a+b+d+e,b+c+e+f},{d+e+g+h,e+f+h+i}}
```

### Scope

Use exact arithmetic to compute the correlation:

```wolfram
ker = Differences[Range[10]^2];
list = 1/Range[20];
```

```wolfram
ListCorrelate[ker, list]
(* Output *)
{(52489)/(2520),(40499)/(2520),(124189)/(9240),(64591)/(5544),(531397)/(51480),(371809)/(40040),(55361)/(6552),(86017)/(11088),(5860303)/(816816),(962251)/(144144),(76494941)/(12252240),(5910403)/(1007760)}
```

Use machine arithmetic:

```wolfram
ListCorrelate[N[ker], N[list]]
(* Output *)
{20.828968253968256,16.071031746031746,13.440367965367965,11.650613275613276,10.322397047397049,9.28593906093906,8.449481074481076,7.757665945165945,7.174569303245774,6.675622988122988,6.243343339666869,5.864891442406922}
```

Use 24-digit precision arithmetic:

```wolfram
ListCorrelate[N[ker,24], N[list,24]]
(* Output *)
{20.82896825396825396825396825396825396825,16.07103174603174603174603174603174603176,13.44036796536796536796536796536796536797,11.65061327561327561327561327561327561326,10.32239704739704739704739704739704739705,9.28593906093906093906093906093906093905,8.44948107448107448107448107448107448105,7.75766594516594516594516594516594516592,7.1745693032457738340091281267751855987,6.67562298812298812298812298812298812296,6.24334333966686907863378451613745731389,5.86489144240692228308327379534809875365}
```

Correlation of complex data:

```wolfram
ListCorrelate[{1, I}, RandomComplex[1 + I, 5]]
(* Output *)
{0.7320679396401226+1.1094562546505216 ⅈ,-0.3843498997067414+1.000469396149658 ⅈ,-0.005743311997792766+1.4992688269259575 ⅈ,-0.40645910761714177+1.073586651240094 ⅈ}
```

Two-dimensional correlation:

```wolfram
ListCorrelate[{{1,0,1},{0,-4,0},{1,0,1}}, Array[a_#&, {4,4}]]
(* Output *)
{{2 a_1-4 a_2+2 a_3,2 a_1-4 a_2+2 a_3},{2 a_2-4 a_3+2 a_4,2 a_2-4 a_3+2 a_4}}
```

Cyclic two-dimensional correlation:

```wolfram
c = ListCorrelate[{{1,0,1},{0,-4,0},{1,0,1}}, Array[a_#&, {4,4}],{1,1}];
```

```wolfram
MatrixForm[c]
(* Output *)
({{2 a_1-4 a_2+2 a_3, 2 a_1-4 a_2+2 a_3, 2 a_1-4 a_2+2 a_3, 2 a_1-4 a_2+2 a_3}, {2 a_2-4 a_3+2 a_4, 2 a_2-4 a_3+2 a_4, 2 a_2-4 a_3+2 a_4, 2 a_2-4 a_3+2 a_4}, {2 a_1+2 a_3-4 a_4, 2 a_1+2 a_3-4 a_4, 2 a_1+2 a_3-4 a_4, 2 a_1+2 a_3-4 a_4}, {-4 a_1+2 a_2+2 a_4, -4 a_1+2 a_2+2 a_4, -4 a_1+2 a_2+2 a_4, -4 a_1+2 a_2+2 a_4}})
```

Two-dimensional correlation with maximal overhangs and zero padding:

```wolfram
c = ListCorrelate[{{1,0,1},{0,-4,0},{1,0,1}}, Array[a_#&, {4,4}],{-1,1}, 0];
```

```wolfram
MatrixForm[c]
(* Output *)
({{a_1, a_1, 2 a_1, 2 a_1, a_1, a_1}, {a_2, -4 a_1+a_2, -4 a_1+2 a_2, -4 a_1+2 a_2, -4 a_1+a_2, a_2}, {a_1+a_3, a_1-4 a_2+a_3, 2 a_1-4 a_2+2 a_3, 2 a_1-4 a_2+2 a_3, a_1-4 a_2+a_3, a_1+a_3}, {a_2+a_4, a_2-4 a_3+a_4, 2 a_2-4 a_3+2 a_4, 2 a_2-4 a_3+2 a_4, a_2-4 a_3+a_4, a_2+a_4}, {a_3, a_3-4 a_4, 2 a_3-4 a_4, 2 a_3-4 a_4, a_3-4 a_4, a_3}, {a_4, a_4, 2 a_4, 2 a_4, a_4, a_4}})
```

### Generalizations & Extensions

Use functions `f` and `g` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html) and [Times](https://reference.wolfram.com/language/ref/Times.html):

```wolfram
ListCorrelate[{x,y,z},{1,2,3,4,5},{1,-1},{x,y,z}, f,g]
(* Output *)
{g[f[x,1],f[y,2],f[z,3]],g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]]}
```

Use functions `f` and `g` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html) and [Times](https://reference.wolfram.com/language/ref/Times.html) with maximal overhangs and zero padding:

```wolfram
ListCorrelate[{x,y,z},{1,2,3,4,5},{-1,1},0, f,g]
(* Output *)
{g[f[x,0],f[y,0],f[z,1]],g[f[x,0],f[y,1],f[z,2]],g[f[x,1],f[y,2],f[z,3]],g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]],g[f[x,4],f[y,5],f[z,0]],g[f[x,5],f[y,0],f[z,0]]}
```

Use functions `f` and `g` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html) and [Times](https://reference.wolfram.com/language/ref/Times.html) with maximal overhangs and empty padding:

```wolfram
ListCorrelate[{x,y,z},{1,2,3,4,5},{-1,1},{}, f,g]
(* Output *)
{g[f[x],f[y],f[z,1]],g[f[x],f[y,1],f[z,2]],g[f[x,1],f[y,2],f[z,3]],g[f[x,2],f[y,3],f[z,4]],g[f[x,3],f[y,4],f[z,5]],g[f[x,4],f[y,5],f[z]],g[f[x,5],f[y],f[z]]}
```

[ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html) works with [TimeSeries](https://reference.wolfram.com/language/ref/TimeSeries.html):

```wolfram
ts=TimeSeries[Table[x_t,{t,12}]]
(* Output *)
TimeSeries[...]
```

```wolfram
ts["Path"]
(* Output *)
{{0,x_1},{1,x_2},{2,x_3},{3,x_4},{4,x_5},{5,x_6},{6,x_7},{7,x_8},{8,x_9},{9,x_10},{10,x_11},{11,x_12}}
```

```wolfram
ListCorrelate[{1/2,1/2},ts,{1}]
(* Output *)
TimeSeries[...]
```

```wolfram
%["Path"]
(* Output *)
{{1,(x_2)/(2)+(x_3)/(2)},{2,(x_3)/(2)+(x_4)/(2)},{3,(x_4)/(2)+(x_5)/(2)},{4,(x_5)/(2)+(x_6)/(2)},{5,(x_6)/(2)+(x_7)/(2)},{6,(x_7)/(2)+(x_8)/(2)},{7,(x_8)/(2)+(x_9)/(2)},{8,(x_9)/(2)+(x_10)/(2)},{9,(x_10)/(2)+(x_11)/(2)},{10,(x_11)/(2)+(x_12)/(2)},{11,(1)/(2)+(x_12)/(2)}}
```

### Applications

Smooth data with a weighted running average:

```wolfram
x = N[Range[1000]/1000];
data = RandomReal[{-1,1}, 1000] +2  Sin[7 x];
```

Normalized Gaussian profile for averaging weights:

```wolfram
w = Exp[-.01 N[Range[-24,24]^2]];
w /= Total[w];
```

```wolfram
smoother = ListCorrelate[w, data];
```

```wolfram
Show[ListLinePlot[Transpose[{x,data}]], ListLinePlot[Transpose[{Take[x,{25,-25}],smoother}], PlotStyle->Red]]
```

*([Graphics])*

Gaussian smoothing of an image:

```wolfram
img = 256-Import["ExampleData/turtle.jpg", "Data"];
ArrayPlot[img]
```

*([Graphics])*

Gaussian kernel with a 5×5 pixel stencil:

```wolfram
g[σ_][x_, y_] := (1)/(2 π σ^2)Exp[-(1)/(2σ^2)(x^2 + y^2)];
gk = N[Outer[g[1], Range[-2,2], Range[-2,2]]];
```

Smooth the image:

```wolfram
smooth = ListCorrelate[gk, img];
```

```wolfram
ArrayPlot[smooth]
(* Output *)
![image](img/image_001.png)
```

Edge detection in an image:

```wolfram
img = 256-Import["ExampleData/ocelot.jpg", "Data"];
ArrayPlot[img]
```

*([Graphics])*

Correlate with a Laplacian filter kernel:

```wolfram
dimg = ListCorrelate[{{1,0,1},{0,-4,0},{1,0,1}},img];
```

```wolfram
ArrayPlot[dimg]
```

*([Graphics])*

Use a Laplacian of a Gaussian filter kernel:

```wolfram
g[σ_][x_, y_] := (1)/(2 π σ^2)Exp[-(1)/(2σ^2)(x^2 + y^2)];
f[σ_][x_, y_] = D[g[σ][x,y],{x,2}] + D[g[σ][x,y],{y,2}];
```

```wolfram
ker = Outer[f[1.5], N[Range[-5,5]], N[Range[-5,5]]];
ListPlot3D[ker, PlotRange->All]
```

*([Graphics3D])*

```wolfram
dimg = ListCorrelate[ker,img];
```

```wolfram
ArrayPlot[dimg]
(* Output *)
![image](img/image_003.png)
```

Generate Pascal's triangle:

```wolfram
NestList[ListCorrelate[{1,1},#,{-1,1},0]&,{1},10]
(* Output *)
{{1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1},{1,6,15,20,15,6,1},{1,7,21,35,35,21,7,1},{1,8,28,56,70,56,28,8,1},{1,9,36,84,126,126,84,36,9,1},{1,10,45,120,210,252,210,120,45,10,1}}
```

```wolfram
Column[%,Center]
(* Output *)
{{{1}}, {{1,1}}, {{1,2,1}}, {{1,3,3,1}}, {{1,4,6,4,1}}, {{1,5,10,10,5,1}}, {{1,6,15,20,15,6,1}}, {{1,7,21,35,35,21,7,1}}, {{1,8,28,56,70,56,28,8,1}}, {{1,9,36,84,126,126,84,36,9,1}}, {{1,10,45,120,210,252,210,120,45,10,1}}}
```

Additive cellular automata:

```wolfram
{ArrayPlot[PadRight[Mod[NestList[ListCorrelate[{1,1},#,{-1,1},0]&,{1},100],2]]],ArrayPlot[PadRight[Mod[NestList[ListCorrelate[{1,2},#,{-1,1},0]&,{1},100],5]]]}
(* Output *)
{[Graphics],[Graphics]}
```

Apply a finite difference formula to a uniformly sampled function:

```wolfram
n = 100;h = 1./n;
x = N[h Range[n]];
f =Sin[2 π x];
```

```wolfram
ListPlot[{f, ListCorrelate[({-1,0,1})/(2 h),f,2]}, DataRange->{0,1}]
```

*([Graphics])*

Show the error for different numbers of grid points:

```wolfram
Table[n = 10^k; h = 1./n;x = N[h Range[n]];
ListPlot[2π Cos[2 π x] - ListCorrelate[({-1,0,1})/(2 h),Sin[2 π x],2], DataRange->{h,1}] ,
{k,1,3}]
(* Output *)
{[Graphics],[Graphics],[Graphics]}
```

Show the error for different numbers of grid points for a second derivative approximation:

```wolfram
Table[n = 10^k; h = 1./n;x = N[h Range[n]];
ListPlot[-4 π^2Sin[2 π x] - ListCorrelate[({-1,16,-30,16,-1})/(12 h^2),Sin[2 π x],3], DataRange->{h,1}] ,
{k,1,3}]
(* Output *)
{[Graphics],[Graphics],[Graphics]}
```

### Properties & Relations

[ListCorrelate](https://reference.wolfram.com/language/ref/ListCorrelate.html) is equivalent to [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) with the kernel reversed:

```wolfram
ker = RandomReal[1, {123}];
list = RandomReal[1, {10000}];
```

```wolfram
ListConvolve[ker, list] == ListCorrelate[Reverse[ker], list]
(* Output *)
True
```

Generate two random vectors:

```wolfram
{a, b} = RandomReal[1, {2,31}];
```

A function for constructing a circulant matrix from a vector:

```wolfram
circulant[v_]:=ToeplitzMatrix[v,RotateRight[Reverse[v]]]
```

Cyclic correlation is equivalent to multiplication with a circulant matrix:

```wolfram
ListCorrelate[a, b, {-1,-1}]==circulant[Reverse[a]].b
(* Output *)
True
```

Cyclic correlation is also equivalent to multiplication in the discrete Fourier transform domain:

```wolfram
ListCorrelate[a, b, {-1,-1}] ==InverseFourier[Fourier[Reverse[a]] Fourier[b]] Sqrt[31]
(* Output *)
True
```

Generate two random vectors:

```wolfram
{a, b} = RandomReal[1, {2,31}];
```

A function for constructing an upper triangular Toeplitz matrix from a vector:

```wolfram
uttoep[v_]:=With[{n=Length[v]},ToeplitzMatrix[First[v]UnitVector[n,1],v]]
```

Cyclic correlation with zero-padding is equivalent to multiplication with an upper triangular Toeplitz matrix:

```wolfram
ListCorrelate[a,b,1,0]==uttoep[a].b
(* Output *)
True
```

## Tech Notes ▪Convolutions and Correlations ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Data Transforms and Smoothing ▪Linear and Nonlinear Filters ▪Statistical Data Analysis ▪Signal Filtering & Filter Design

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=ListCorrelate) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1999 (4.0)
