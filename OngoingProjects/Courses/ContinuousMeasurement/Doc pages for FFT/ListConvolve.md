# ListConvolve | [SpanFromLeft]

> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*] — forms the convolution of the kernel `*ker*` with `*list*`.
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*k*] — forms the cyclic convolution in which the `*k*`$^{th}$ element of `*ker*` is aligned with each element in `*list*`.
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{*k*_*L*,*k*_*R*}] — forms the cyclic convolution whose first element contains `*list*[[1]]*ker*[[*k*_*L*]]` and whose last element contains `*list*[[-1]]*ker*[[*k*_*R*]]`.
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*klist*,*p*] — forms the convolution in which `*list*` is padded at each end with repetitions of the element `*p*`.
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*klist*,{*p*_1,*p*_2,…}] — forms the convolution in which `*list*` is padded at each end with cyclic repetitions of the `*p*_*i*`.
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*klist*,*padding*,*g*,*h*] — forms a generalized convolution in which `*g*` is used in place of [Times](https://reference.wolfram.com/language/ref/Times.html) and `*h*` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html).
> [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*klist*,*padding*,*g*,*h*,*lev*] — forms a convolution using elements at level `*lev*` in `*ker*` and `*list*`.

## Details

With kernel `*K*_*r*` and list `*a*_*s*`, [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*] computes $\sum_{r}K_{r}a_{s-r}$, where the limits of the sum are such that the kernel never overhangs either end of the list.

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*] gives a result of length [Length](https://reference.wolfram.com/language/ref/Length.html)[*list*]-[Length](https://reference.wolfram.com/language/ref/Length.html)[*ker*]+1.

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*] allows no overhangs and is equivalent to [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{-1,1}].

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,*k*] is equivalent to [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{*k*,*k*}].

The values of `*k*_*L*` and `*k*_*R*` in [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{*k*_*L*,*k*_*R*}] determine the amount of overhang to allow at each end of `*list*`.

Common settings for `{*k*_*L*,*k*_*R*}` are:

{-1,1} | no overhangs (default)
{-1,-1} | maximal overhang at the right[Hyphen]hand end
{1,1} | maximal overhang at the left[Hyphen]hand end
{1,-1} | maximal overhangs at both beginning and end

With maximal overhang at one end only, the result from [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) is the same length as `*list*`.

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{*k*_*L*,*k*_*R*},*padlist*] effectively lays down repeated copies of `*padlist*`, then superimposes one copy of `*list*` on them and forms a convolution of the result.

Common settings for `*padlist*` are:

*p* | pad with repetitions of a single element
{*p*_1,*p*_2,…} | pad with cyclic repetitions of a sequence of elements
*list* | pad by treating `*list*` as cyclic (default)
{} | do no padding

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) works with multidimensional kernels and lists of data.

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html)[*ker*,*list*,{{*k*_*L***1,*k*_*L***2,…},{*k*_*R***1,*k*_*R***2,…}}] forms the cyclic convolution whose `{1,1,…}` element contains `*ker*[[*k*_*L***1,*k*_*L***2,…]]*list*[[1,1,…]]`, and whose `{-1,-1,…}` element contains `*ker*[[*k*_*R***1,*k*_*R***2,…]]*list*[[-1,-1,…]]`.

`{*k*_*L*,*k*_*R*}` is taken to be equivalent to `{{*k*_*L*,*k*_*L*,…},{*k*_*R*,*k*_*R*,…}}`.

When a function `*h*` is specified to use in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html), explicit nested `*h*` expressions are generated with a depth equal to the depth of `*ker*`.

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) works with exact numbers and symbolic data as well as approximate numbers.

## Examples

### Basic Examples

Convolve a kernel `{x,y}` with a list of data:

```wolfram
ListConvolve[{x,y},{a,b,c,d,e,f}]
(* Output *)
{b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}
```

Make a cyclic convolution the same length as the original data:

```wolfram
ListConvolve[{x,y},{a,b,c,d,e,f},1]
(* Output *)
{a x+f y,b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}
```

Align element 2 in the kernel with successive elements in the data:

```wolfram
ListConvolve[{x,y},{a,b,c,d,e,f},2]
(* Output *)
{b x+a y,c x+b y,d x+c y,e x+d y,f x+e y,a x+f y}
```

Pad with `zzz` instead of using the data cyclically:

```wolfram
ListConvolve[{x,y},{a,b,c,d,e,f},1,zzz]
(* Output *)
{a x+y zzz,b x+a y,c x+b y,d x+c y,e x+d y,f x+e y}
```

Two-dimensional convolution:

```wolfram
ListConvolve[{{1,1},{1,1}},{{a,b,c},{d,e,f},{g,h,i}}]
(* Output *)
{{a+b+d+e,b+c+e+f},{d+e+g+h,e+f+h+i}}
```

### Scope

#### Overhangs and Alignments

"Slide" the kernel along the data, allowing no overhangs:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6}]
(* Output *)
{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}
```

Maximal overhang at the beginning; none at the end:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},1]
(* Output *)
{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}
```

Maximal overhang at the end; none at the beginning:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},-1]
(* Output *)
{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}
```

Maximal overhangs at both beginning and end:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},{1,-1}]
(* Output *)
{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}
```

Align element 1 of the kernel with the first element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},1]
(* Output *)
{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}
```

Align element 2 of the kernel with the first element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},2]
(* Output *)
{2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z}
```

Align element 3 of the kernel with the first element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},3]
(* Output *)
{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}
```

Align the last element of the kernel with the first element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},-1]
(* Output *)
{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}
```

Align the first element of the kernel with both the first and last elements of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},{1,1}]
(* Output *)
{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}
```

Align element 2 of the kernel with the first element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},{2,1}]
(* Output *)
{2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}
```

Align element 2 of the kernel with the last element of the data:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5,6},{1,2}]
(* Output *)
{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z}
```

#### Data Padding

Use padding `aa`:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5},1,aa]
(* Output *)
{x+aa y+aa z,2 x+y+aa z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}
```

Cyclically use a list of padding elements:

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5},1,{aa,bb}]
(* Output *)
{x+bb y+aa z,2 x+y+bb z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}
```

```wolfram
ListConvolve[{x,y,z},{1,2,3,4,5},1,{aa,bb,cc}]
(* Output *)
{x+cc y+bb z,2 x+y+cc z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z}
```

```wolfram
ListConvolve[{x,y,z},{1,2,3},{1,-1},{aa,bb,cc}]
(* Output *)
{x+cc y+bb z,2 x+y+cc z,3 x+2 y+z,aa x+3 y+2 z,bb x+aa y+3 z}
```

#### Higher Dimensions

Allow no overhangs:

```wolfram
ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}}]
(* Output *)
{{5 a+4 b+2 c+d,6 a+5 b+3 c+2 d},{8 a+7 b+5 c+4 d,9 a+8 b+6 c+5 d}}
```

Align with the `{1,1}` elements of the kernel and data:

```wolfram
ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},1]
(* Output *)
{{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d},{4 a+6 b+c+3 d,5 a+4 b+2 c+d,6 a+5 b+3 c+2 d},{7 a+9 b+4 c+6 d,8 a+7 b+5 c+4 d,9 a+8 b+6 c+5 d}}
```

The result has the same dimensions as the input data:

```wolfram
Dimensions[%]
(* Output *)
{3,3}
```

Give a different overhang in each dimension:

```wolfram
ListConvolve[{{a,b},{c,d}},{{1,2,3},{4,5,6},{7,8,9}},{{1,1},{2,2}}]
(* Output *)
{{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d,a+3 b+7 c+9 d},{4 a+6 b+c+3 d,5 a+4 b+2 c+d,6 a+5 b+3 c+2 d,4 a+6 b+c+3 d},{7 a+9 b+4 c+6 d,8 a+7 b+5 c+4 d,9 a+8 b+6 c+5 d,7 a+9 b+4 c+6 d},{a+3 b+7 c+9 d,2 a+b+8 c+7 d,3 a+2 b+9 c+8 d,a+3 b+7 c+9 d}}
```

### Generalizations & Extensions

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) works with sparse arrays:

```wolfram
ListConvolve[{a,b,c},SparseArray[6->x,20]]
(* Output *)
{0,0,0,a x,b x,c x,0,0,0,0,0,0,0,0,0,0,0,0}
```

Use `f` in place of [Times](https://reference.wolfram.com/language/ref/Times.html):

```wolfram
ListConvolve[{x,y},{1,2,3,4},1,0,f]
(* Output *)
{f[x,1]+f[y,0],f[x,2]+f[y,1],f[x,3]+f[y,2],f[x,4]+f[y,3]}
```

Use `g` in place of [Plus](https://reference.wolfram.com/language/ref/Plus.html):

```wolfram
ListConvolve[{x,y},{1,2,3,4},1,0,f,g]
(* Output *)
{g[f[y,0],f[x,1]],g[f[y,1],f[x,2]],g[f[y,2],f[x,3]],g[f[y,3],f[x,4]]}
```

```wolfram
ListConvolve[{x,y},{1,2,3,4},1,0,f,List]
(* Output *)
{{f[y,0],f[x,1]},{f[y,1],f[x,2]},{f[y,2],f[x,3]},{f[y,3],f[x,4]}}
```

Use `f` and `g` in place of [Times](https://reference.wolfram.com/language/ref/Times.html) and [Plus](https://reference.wolfram.com/language/ref/Plus.html), with empty data padding:

```wolfram
ListConvolve[{x,y,z},{1,2,3},1,{},f,g]
(* Output *)
{g[f[z],f[y],f[x,1]],g[f[z],f[y,1],f[x,2]],g[f[z,1],f[y,2],f[x,3]]}
```

```wolfram
ListConvolve[{x,y,z},{1,2,3},-1,{},f,g]
(* Output *)
{g[f[z,1],f[y,2],f[x,3]],g[f[z,2],f[y,3],f[x]],g[f[z,3],f[y],f[x]]}
```

```wolfram
ListConvolve[{x,y,z},{1,2,3},{1,-1},{},f,g]
(* Output *)
{g[f[z],f[y],f[x,1]],g[f[z],f[y,1],f[x,2]],g[f[z,1],f[y,2],f[x,3]],g[f[z,2],f[y,3],f[x]],g[f[z,3],f[y],f[x]]}
```

[ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) works with [TimeSeries](https://reference.wolfram.com/language/ref/TimeSeries.html):

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
ListConvolve[{1/2,1/2},ts,{1}]
(* Output *)
TimeSeries[...]
```

```wolfram
%["Path"]
(* Output *)
{{0,(x_1)/(2)+(x_12)/(2)},{1,(x_1)/(2)+(x_2)/(2)},{2,(x_2)/(2)+(x_3)/(2)},{3,(x_3)/(2)+(x_4)/(2)},{4,(x_4)/(2)+(x_5)/(2)},{5,(x_5)/(2)+(x_6)/(2)},{6,(x_6)/(2)+(x_7)/(2)},{7,(x_7)/(2)+(x_8)/(2)},{8,(x_8)/(2)+(x_9)/(2)},{9,(x_9)/(2)+(x_10)/(2)},{10,(x_10)/(2)+(x_11)/(2)},{11,(x_11)/(2)+(x_12)/(2)}}
```

### Applications

Find a moving average:

```wolfram
ListConvolve[{1,1}/2,{a,b,c,d,e}]
(* Output *)
{(a)/(2)+(b)/(2),(b)/(2)+(c)/(2),(c)/(2)+(d)/(2),(d)/(2)+(e)/(2)}
```

Or use the [MovingAverage](https://reference.wolfram.com/language/ref/MovingAverage.html) function:

```wolfram
MovingAverage[{a,b,c,d,e},2]
(* Output *)
{(a+b)/(2),(b+c)/(2),(c+d)/(2),(d+e)/(2)}
```

Smooth noisy data:

```wolfram
data=Table[Sqrt[x]+RandomReal[],{x,30}];
```

```wolfram
ListLinePlot[data]
```

*([Graphics])*

```wolfram
ListLinePlot[ListConvolve[{1,1,2,1,1}/6,data]]
```

*([Graphics])*

Find the autocorrelation of a list:

```wolfram
data=Table[Mod[i^2,17],{i,100}];
```

```wolfram
ListLinePlot[ListConvolve[data,data,{1,1}]]
```

*([Graphics])*

Apply a simple image processing filter:

```wolfram
ArrayPlot[CellularAutomaton[150,{{1},0},30]]
```

*([Graphics])*

```wolfram
ArrayPlot[ListConvolve[{{-1,1},{1,-1}},CellularAutomaton[150,{{1},0},30]]]
```

*([Graphics])*

Multiply polynomials by convolving coefficient lists:

```wolfram
CoefficientList[(a+b x+c x^2)(1+2x+3x^2+4x^3),x]
(* Output *)
{a,2 a+b,3 a+2 b+c,4 a+3 b+2 c,4 b+3 c,4 c}
```

```wolfram
ListConvolve[{a,b,c},{1,2,3,4},{1,-1},0]
(* Output *)
{a,2 a+b,3 a+2 b+c,4 a+3 b+2 c,4 b+3 c,4 c}
```

Multiply numbers by convolving digit lists:

```wolfram
ListConvolve[IntegerDigits[24561],IntegerDigits[1247],{1,-1},0]
(* Output *)
{2,8,21,46,61,61,46,7}
```

```wolfram
FromDigits[%,10]
(* Output *)
30627567
```

```wolfram
24561 1247
(* Output *)
30627567
```

Pascal's triangle:

```wolfram
Column[NestList[ListConvolve[{1,1},#,{1,-1},0]&,{1},5],Center]
(* Output *)
{{{1}}, {{1,1}}, {{1,2,1}}, {{1,3,3,1}}, {{1,4,6,4,1}}, {{1,5,10,10,5,1}}}
```

Additive cellular automaton in base 5:

```wolfram
ArrayPlot[NestList[Mod[ListConvolve[{2,3,1},#,2],5]&,SparseArray[50->1,101],50]]
```

*([Graphics])*

Fast multiplication of high[Hyphen]degree polynomials:

```wolfram
p1[x_] = RandomInteger[9, 101].x^Range[0,100];
p2[x_] = RandomInteger[9, 101].x^Range[0,100];
```

Use [ListConvolve](https://reference.wolfram.com/language/ref/ListConvolve.html) with maximal overhangs and zero padding:

```wolfram
c3 = ListConvolve[CoefficientList[p1[x],x], CoefficientList[p2[x],x],{1,-1}, 0] ;
```

```wolfram
p3[x_] = c3.x^Range[0, Length[c3] - 1];
```

```wolfram
p3[x] == Expand[p1[x] p2[x]]
(* Output *)
True
```

### Properties & Relations

Generate two random vectors:

```wolfram
{a, b} = RandomReal[1, {2,31}];
```

A function for constructing a circulant matrix from a vector:

```wolfram
circulant[v_]:=ToeplitzMatrix[v,RotateRight[Reverse[v]]]
```

Cyclic convolution is equivalent to multiplication with a circulant matrix:

```wolfram
ListConvolve[a, b, {1,1}]==circulant[a].b
(* Output *)
True
```

Cyclic convolution is also equivalent to multiplication in the discrete Fourier transform domain:

```wolfram
ListConvolve[a, b, {1,1}] ==InverseFourier[Fourier[a] Fourier[b]] Sqrt[31]
(* Output *)
True
```

Generate two random vectors:

```wolfram
{a, b} = RandomReal[1, {2,31}];
```

A function for constructing a lower triangular Toeplitz matrix from a vector:

```wolfram
lttoep[v_]:=With[{n=Length[v]},ToeplitzMatrix[v,First[v]UnitVector[n,1]]]
```

Cyclic convolution with zero-padding is equivalent to multiplication with a lower triangular Toeplitz matrix:

```wolfram
ListConvolve[a,b,1,0]==lttoep[a].b
(* Output *)
True
```

Convolve with a single element:

```wolfram
ListConvolve[{1},{1,2,3,4,5,6}]
(* Output *)
{1,2,3,4,5,6}
```

```wolfram
ListConvolve[{1},{1,2,3,4,5,6},2]
(* Output *)
{2,3,4,5,6,1}
```

A kernel of the same length as the data, with no overhangs, is like a reversed dot product:

```wolfram
ListConvolve[{x,y,z},{a,b,c}]
(* Output *)
{c x+b y+a z}
```

Successive differences:

```wolfram
ListConvolve[{1,-1},{a,b,c,d,e}]
(* Output *)
{-a+b,-b+c,-c+d,-d+e}
```

```wolfram
Differences[{a,b,c,d,e}]
(* Output *)
{-a+b,-b+c,-c+d,-d+e}
```

Align with successive elements:

```wolfram
Table[ListConvolve[{x,y,z},{1,2,3,4,5,6},i],{i,3}]//Column
(* Output *)
{{{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}}, {{2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z}}, {{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}}}
```

```wolfram
Table[ListConvolve[{x,y,z},{1,2,3,4,5,6},i],{i,-3,-1}]//Column
(* Output *)
{{{x+6 y+5 z,2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z}}, {{2 x+y+6 z,3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z}}, {{3 x+2 y+z,4 x+3 y+2 z,5 x+4 y+3 z,6 x+5 y+4 z,x+6 y+5 z,2 x+y+6 z}}}
```

Varying alignments:

```wolfram
Table[ListConvolve[{1,1},{1,1,1,1},{i,j}],{i,4},{j,4}]//Grid
(* Output *)
{{{2,2,2,2}, {2,2,2,2,2}, {2,2,2,2,2,2}, {2,2,2,2,2,2,2}}, {{2,2,2}, {2,2,2,2}, {2,2,2,2,2}, {2,2,2,2,2,2}}, {{2,2}, {2,2,2}, {2,2,2,2}, {2,2,2,2,2}}, {{2}, {2,2}, {2,2,2}, {2,2,2,2}}}
```

```wolfram
Grid[Table[ListConvolve[{1,1},{1,1,1,1},{i,-j}],{i,4},{j,4}]]
(* Output *)
ListConvolve
(* Output *)
{{{2,2,2,2,2}, {2,2,2,2}, {2,2,2}, {2,2}}, {{2,2,2,2}, {2,2,2}, {2,2}, {2}}, {{2,2,2}, {2,2}, {2}, {}}, {{2,2}, {2}, {}, {}}}
```

```wolfram
Grid[Table[ListConvolve[{1,1},{1,1,1,1},{-i,j}],{i,4},{j,4}]]
(* Output *)
{{{2,2,2}, {2,2,2,2}, {2,2,2,2,2}, {2,2,2,2,2,2}}, {{2,2,2,2}, {2,2,2,2,2}, {2,2,2,2,2,2}, {2,2,2,2,2,2,2}}, {{2,2,2,2,2}, {2,2,2,2,2,2}, {2,2,2,2,2,2,2}, {2,2,2,2,2,2,2,2}}, {{2,2,2,2,2,2}, {2,2,2,2,2,2,2}, {2,2,2,2,2,2,2,2}, {2,2,2,2,2,2,2,2,2}}}
```

Varying alignments, with padding:

```wolfram
Table[ListConvolve[{1,1},{1,1,1,1},{i,j},0],{i,4},{j,4}]//Grid
(* Output *)
{{{1,2,2,2}, {1,2,2,2,1}, {1,2,2,2,1,0}, {1,2,2,2,1,0,0}}, {{2,2,2}, {2,2,2,1}, {2,2,2,1,0}, {2,2,2,1,0,0}}, {{2,2}, {2,2,1}, {2,2,1,0}, {2,2,1,0,0}}, {{2}, {2,1}, {2,1,0}, {2,1,0,0}}}
```

```wolfram
Table[ListConvolve[{1,1},{1,1,1,1},{i,-j},0],{i,4},{j,4}]//Grid
(* Output *)
ListConvolve
(* Output *)
{{{1,2,2,2,1}, {1,2,2,2}, {1,2,2}, {1,2}}, {{2,2,2,1}, {2,2,2}, {2,2}, {2}}, {{2,2,1}, {2,2}, {2}, {}}, {{2,1}, {2}, {}, {}}}
```

```wolfram
Table[ListConvolve[{1,1},{1,1,1,1},{-i,j},0],{i,4},{j,4}]//Grid
(* Output *)
{{{2,2,2}, {2,2,2,1}, {2,2,2,1,0}, {2,2,2,1,0,0}}, {{1,2,2,2}, {1,2,2,2,1}, {1,2,2,2,1,0}, {1,2,2,2,1,0,0}}, {{0,1,2,2,2}, {0,1,2,2,2,1}, {0,1,2,2,2,1,0}, {0,1,2,2,2,1,0,0}}, {{0,0,1,2,2,2}, {0,0,1,2,2,2,1}, {0,0,1,2,2,2,1,0}, {0,0,1,2,2,2,1,0,0}}}
```

### Possible Issues

By default, the output is not the same length as the input:

```wolfram
ListConvolve[{1,1},{1,2,3,4,5}]
(* Output *)
{3,5,7,9}
```

Use overhangs to make the output the same length as the input:

```wolfram
ListConvolve[{1,1},{1,2,3,4,5},1]
(* Output *)
{6,3,5,7,9}
```

Use overhangs to do the same in 2D:

```wolfram
ListConvolve[{{1,1},{1,1}},{{1,2,3},{4,5,6},{7,8,9}},1]
(* Output *)
{{20,18,22},{14,12,16},{26,24,28}}
```

## Tech Notes ▪Convolutions and Correlations ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Data Transforms and Smoothing ▪Signal Processing ▪Signal Filtering & Filter Design ▪Handling Arrays of Data ▪Linear and Nonlinear Filters ▪Image Filtering & Neighborhood Processing ▪Fourier Analysis ▪Numerical Data ▪Scientific Data Analysis ▪Summation Transforms

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=ListConvolve) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1999 (4.0)
