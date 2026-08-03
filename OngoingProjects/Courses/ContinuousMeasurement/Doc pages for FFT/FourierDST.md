# FourierDST | [SpanFromLeft]

> [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*]  — finds the Fourier discrete sine transform of a list of real numbers.
> [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*,*m*] — finds the Fourier discrete sine transform of type $m$.

## Details

Possible types $m$ of discrete sine transform for a list $u_{r}$ of length $n$ giving a result $v_{s}$ are:

1 (DST-I) | $\sqrt{\frac{2}{n+1}}\sum_{r=1}^{n}u_{r} sin(\frac{\pi}{n+1} r s)$
2 (DST-II) | $\frac{1}{\sqrt{n}}\sum_{r=1}^{n}u_{r} sin(\frac{\pi}{n} (r-\frac{1}{2}) s)$
3 (DST-III) | $\frac{1}{\sqrt{n}}(2 \sum_{r=1}^{n-1}u_{r} sin(\frac{\pi}{n} r (s-\frac{1}{2}))+(-1)^{s-1}u_{n} )$
4 (DST-IV) | $\sqrt{\frac{2}{n}}\sum_{r=1}^{n}u_{r} sin(\frac{\pi}{n} (r-\frac{1}{2}) (s-\frac{1}{2}))$

[FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*] is equivalent to [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*,2].

The inverse discrete sine transforms for types 1, 2, 3, 4 are types 1, 3, 2, 4, respectively.

The `*list*` given in [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*] can be nested to represent an array of data in any number of dimensions.

The array of data must be rectangular.

If the elements of `*list*` are exact numbers, [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html) begins by applying [N](https://reference.wolfram.com/language/ref/N.html) to them.

[FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html) can be used on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) objects.

## Examples

### Basic Examples

Find a discrete sine transform:

```wolfram
FourierDST[{0,0,1,0,1}]
(* Output *)
{0.5854101966249685,-0.2628655560595667,-0.08541019662496845,-0.42532540417601994,0.8944271909999159}
```

Find the inverse discrete sine transform:

```wolfram
FourierDST[%,3]
(* Output *)
{8.326672684688674×10^-17,2.220446049250313×10^-16,1.,5.551115123125783×10^-17,1.}
```

Find a discrete sine transform of type 1 (DST-I):

```wolfram
FourierDST[{1,0,0,1,2},1]
(* Output *)
{1.3660254037844388,-1.,1.7320508075688779,-8.326672684688674×10^-17,0.36602540378443893}
```

Find the inverse discrete sine transform:

```wolfram
FourierDST[%,1]
(* Output *)
{1.0000000000000004,-2.7755575615628914×10^-17,-1.6653345369377348×10^-16,1.,2.000000000000001}
```

### Scope

Use machine arithmetic to compute the discrete sine transform:

```wolfram
v = {0,1,2,3,4,3,2,1,0};
```

```wolfram
FourierDST[v]
(* Output *)
{4.567444499063787,4.199928165665723×10^-17,-0.9999999999999999,-2.8970898488349384×10^-16,0.06644681695159468,8.508276548579996×10^-17,-0.36610868398461766,-9.73613550891983×10^-17,0.}
```

Use 24[Hyphen]digit precision arithmetic:

```wolfram
FourierDST[N[v, 24]]
(* Output *)
{4.5674444990637874817114087553900185051,0`22.823846159168234,-0.99999999999999999999999999999999999999,0`22.88179328403019,0.0664468169515947902244210433617284481,0`22.834325851418544,-0.36610868398461772806417020124825304679,0`22.713935290882098,0`23.273001272063738}
```

Two-dimensional discrete sine transform:

```wolfram
m = RandomReal[1, {4,3}];
```

```wolfram
FourierDST[m]
(* Output *)
{{0.7213936819102706,0.19342213352502607,0.5113942721526287},{-0.05764909182384692,-0.09612452062873043,0.2548177028691185},{0.12725451865023552,-0.1930021353027679,0.07764989835205499},{0.37730944142035794,-0.33650863767682127,-0.1581014916257782}}
```

Four-dimensional discrete sine transform:

```wolfram
m = RandomReal[1, {4,4,4,4}];
```

```wolfram
Max[FourierDST[m]]
(* Output *)
1.416910689896044
```

### Generalizations & Extensions

The list may have complex values:

```wolfram
FourierDST[{1,2I,3,4I}]
(* Output *)
{1.577161014949475+1.6892463972414662 ⅈ,-0.7071067811865475-0.7071067811865475 ⅈ,-0.11208538229199139+1.4650756326574836 ⅈ,2.-3. ⅈ}
```

You can use `"I"`, `"II"`, `"III"`, or `"IV"` for the types `1, `2, `3, and `4 respectively:

```wolfram
FourierDST[{0,0,1,0,0},"IV"]
(* Output *)
{0.4472135954999578,0.44721359549995815,-0.44721359549995804,-0.44721359549995804,0.44721359549995804}
```

### Applications

#### Sine Series Expansion

Get an expansion for an odd function as a sum of sines:

```wolfram
f[x_]:=Exp[-100(x-1/2)^2];
```

The function values on a uniformly spaced grid with `*n*` points on `[-L,L)`:

```wolfram
n=20;
xg=N[Range[n-1]]/n;
fg=Map[f,xg];
fp=ListPlot[Transpose[{xg,fg}],PlotRange->All]
```

*([Graphics])*

Compute the DST-I and renormalize:

```wolfram
coef=FourierDST[fg,1]/Sqrt[n/2];
```

The function has, in effect, been periodized with a particular odd symmetry:

```wolfram
Show[fp,Plot[Sum[coef[[r]]*Sin[Pi r x],{r,n-1}],{x,-1,1}, PlotRange->All]]
```

*([Graphics])*

Plot the expansion error where the points are defined:

```wolfram
Plot[f[x]-Sum[coef[[r]]*Sin[Pi r x],{r,n-1}],{x,0,1-1/n},PlotRange->All]
```

*([Graphics])*

#### Pseudospectral PDE Discretization

Approximate the second derivative for a function with zero boundary conditions:

```wolfram
uxx[u_?VectorQ] := Module[{n = Length[u]},
FourierDST[-(Pi N[Range[n]/(n+1)])^2 FourierDST[u, 1], 1]]
```

```wolfram
n = 20;
xg = N[Range[n]]/(n+1);
u0 = 1 - 2 Abs[xg - .5];
ListPlot[{Transpose[{xg,u0}], Transpose[{xg, uxx[u0]}]}]
```

*([Graphics])*

Solve the wave equation for a plucked string:

```wolfram
sol = First[NDSolve[{u'[t] == v[t], v'[t] == uxx[u[t]], u[0] == u0, v[0] == ConstantArray[0.,n]},
{u,v}, {t,0,25}]]
(* Output *)
{u->InterpolatingFunction[...],v->InterpolatingFunction[...]}
```

Plot the solution $u(x,t)$ as a surface:

```wolfram
ListPlot3D[Table[u[t] /. sol,{t,0,25}], Mesh->False,DataRange->{{0,25},{0,1}}]
```

*([Graphics3D])*

### Properties & Relations

DST-I and DST-IV are their own inverses:

```wolfram
data=RandomReal[1,9];
```

```wolfram
Chop[FourierDST[FourierDST[data,1],1]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0}
```

```wolfram
Chop[FourierDST[FourierDST[data,4],4]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0}
```

DST-II and DST-III are inverses of each other:

```wolfram
data=RandomReal[1,15];
```

```wolfram
Chop[FourierDST[FourierDST[data,2],3]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
```

```wolfram
Chop[FourierDST[FourierDST[data,3],2]-data]
(* Output *)
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
```

The DST is equivalent to matrix multiplication:

```wolfram
dstII[n_]:=(1)/(Sqrt[n])Table[Sin[Pi(r-1/2)s/n],{s,n},{r,n}]
```

```wolfram
MatrixForm[dstII[7]]
(* Output *)
({{(Sin[(π)/(14)])/(Sqrt[7]), (Sin[(3 π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), (1)/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), (Sin[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(14)])/(Sqrt[7])}, {(Sin[(π)/(7)])/(Sqrt[7]), (Cos[(π)/(14)])/(Sqrt[7]), (Cos[(3 π)/(14)])/(Sqrt[7]), 0, -(Cos[(3 π)/(14)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7]), -(Sin[(π)/(7)])/(Sqrt[7])}, {(Sin[(3 π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), -(1)/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7]), (Sin[(3 π)/(14)])/(Sqrt[7])}, {(Cos[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7]), 0, (Cos[(π)/(14)])/(Sqrt[7]), -(Sin[(π)/(7)])/(Sqrt[7]), -(Cos[(3 π)/(14)])/(Sqrt[7])}, {(Cos[(π)/(7)])/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), (1)/(Sqrt[7]), -(Sin[(3 π)/(14)])/(Sqrt[7]), -(Sin[(π)/(14)])/(Sqrt[7]), (Cos[(π)/(7)])/(Sqrt[7])}, {(Cos[(π)/(14)])/(Sqrt[7]), -(Cos[(3 π)/(14)])/(Sqrt[7]), (Sin[(π)/(7)])/(Sqrt[7]), 0, -(Sin[(π)/(7)])/(Sqrt[7]), (Cos[(3 π)/(14)])/(Sqrt[7]), -(Cos[(π)/(14)])/(Sqrt[7])}, {(1)/(Sqrt[7]), -(1)/(Sqrt[7]), (1)/(Sqrt[7]), -(1)/(Sqrt[7]), (1)/(Sqrt[7]), -(1)/(Sqrt[7]), (1)/(Sqrt[7])}})
```

```wolfram
data=RandomReal[1,7];
```

```wolfram
Chop[dstII[7].data-FourierDST[data]]
(* Output *)
{0,0,0,0,0,0,0}
```

### Possible Issues

[FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html) always returns normalized results:

```wolfram
TableForm[Table[FourierDST[ConstantArray[1,5],type],{type,1,4}]]
(* Output *)
{{2.154700538379252, 0., 0.577350269189626, 0., 0.15470053837925146}, {1.4472135954999579, 1.2102406750369443×10^-16, 0.5527864045000419, -1.0039712117215771×10^-16, 0.4472135954999579}, {2.8235955159711317, 0.8777061007329484, 0.44721359549995776, 0.22786670826713562, 0.07083167502878462}, {2.021471201601977, 0.6965515053690705, 0.4472135954999581, 0.35491071886919645, 0.32016958489789715}}
```

To get unnormalized results, you can multiply by the normalization:

```wolfram
nc[n_,1]=Sqrt[(n+1)/2];
nc[n_,2|3]=Sqrt[n];
nc[n_,4]=Sqrt[n/2];
```

```wolfram
unnormalizedDST[data_,type_]:=FourierDST[data,type]*nc[Length[data],type]
```

```wolfram
TableForm[Table[unnormalizedDST[ConstantArray[1,5],type],{type,1,4}]]
(* Output *)
{{3.7320508075688776, 0., 1.0000000000000002, 0., 0.2679491924311226}, {3.23606797749979, 2.7061804185178403×10^-16, 1.2360679774997894, -2.24494787686228×10^-16, 1.}, {6.313751514675044, 1.9626105055051508, 0.9999999999999997, 0.5095254494944286, 0.1583844403245368}, {3.196226610749831, 1.1013446322926335, 0.7071067811865478, 0.5611631188171801, 0.5062325628940018}}
```

## Tech Notes ▪Discrete Fourier Transforms

## Related Guides ▪Fourier Analysis ▪Summation Transforms

## History Introduced in 2007 (6.0)
