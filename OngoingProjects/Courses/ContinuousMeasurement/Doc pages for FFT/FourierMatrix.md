# FourierMatrix | [SpanFromLeft]

> [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html)[*n*]  — returns an `*n*`×`*n*` Fourier matrix.

## Details and Options

[FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) of order `*n*` returns a list of the length-`*n*` discrete Fourier transform's basis sequences.

Each entry of the Fourier matrix is by default defined as $F[[r,s]]=\frac{1}{\sqrt{\mathit{n}}}\omega^{(\mathit{r}-1)(\mathit{s}-1)}$, where $\omega=e^{2 \pi \mathit{i}\mathit{ }/\mathit{n}}$.

$$
\begin{pmatrix}
\begin{pmatrix}
F=\frac{1}{\sqrt{n}}\begin{pmatrix}
1 & 1 & 1 & \cdot s & 1 \\
1 & \omega & \omega^{2} & \cdot s & \omega^{n-1} \\
1 & \omega^{2} & \omega^{4} & \cdot s & \omega^{2(n-1)} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & \omega^{n-1} & \omega^{2(n-1)} & \cdot s & \omega^{(n-1)(n-1)}
\end{pmatrix}
\end{pmatrix}
\end{pmatrix}
$$

Rows of the [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) are basis sequences of the discrete Fourier transform.

The result `*F*` of [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html)[*n*] is complex symmetric and unitary, meaning that `*F*^(-1)` is [Conjugate](https://reference.wolfram.com/language/ref/Conjugate.html)[*F*].

The following options can be given:

| [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) | {0,1} | parameters to define the Fourier transform |
| --- | --- | --- |
| [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the structure of the returned matrix |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [Infinity](https://reference.wolfram.com/language/ref/Infinity.html) | precision at which to create entries |

Different choices of definitions for the Fourier matrix can be specified using the option [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html). With the setting [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html)->{*a*,*b*}, the Fourier matrix has entries defined as $F[[r,s]]=\frac{1}{\mathit{n}^{(1-\mathit{a})/2}}\omega^{b(\mathit{r}-1)(\mathit{s}-1)}$, where $\omega=e^{2 \pi \mathit{i}\mathit{ }/\mathit{n}}$.

$$
\begin{pmatrix}
\begin{pmatrix}
F=\frac{1}{n^{(1-a)/2}}\begin{pmatrix}
1 & 1 & 1 & \cdot s & 1 \\
1 & \omega^{b} & \omega^{2 b} & \cdot s & \omega^{b(n-1)} \\
1 & \omega^{2 b} & \omega^{4 b} & \cdot s & \omega^{2 b(n-1)} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & \omega^{b(n-1)} & \omega^{2 b(n-1)} & \cdot s & \omega^{b(n-1)(n-1)}
\end{pmatrix}
\end{pmatrix}
\end{pmatrix}
$$

Some common choices for `{*a*,*b*}` are `{0,1}` (physics), `{-1,1}` (data analysis), `{1,-1}` (signal processing).

Possible settings for [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) include:

[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | automatically choose the representation returned
"Dense" | represent the matrix as a dense matrix
"Structured" | represent the matrix as a structured array
"Symmetric" | represent the matrix as a symmetric matrix
"Unitary" | represent the matrix as a unitary matrix

With [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html)[…,[TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)], a dense matrix is returned if the number of matrix entries is less than a preset threshold, and a structured array is returned otherwise.

The result of [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html)[*n*].*list* is equivalent to [Fourier](https://reference.wolfram.com/language/ref/Fourier.html)[*list*] when `*list*` has length `*n*`. However, the computation of [Fourier](https://reference.wolfram.com/language/ref/Fourier.html)[*list*] is much faster and has less numerical error, unless [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) is kept as a structured array.

For a structured [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) `*sa*`, the following properties `"StyleBox[prop, "TI"]"` can be accessed as `*sa*["StyleBox[prop, "TI"]"]`:

"FourierParameters" | parameters `{*a*,*b*}`
"WorkingPrecision" | precision used internally
"Properties" | list of supported properties
"Structure" | type of structured array
"StructuredData" | internal data stored by the structured array
"StructuredAlgorithms" | list of functions with special methods for the structured array
"Summary" | summary information, represented as a [Dataset](https://reference.wolfram.com/language/ref/Dataset.html)

## Examples

### Basic Examples

A $4 \times 4$ Fourier matrix:

```wolfram
FourierMatrix[4]//MatrixForm
(* Output *)
({{(1)/(2), (1)/(2), (1)/(2), (1)/(2)}, {(1)/(2), (ⅈ)/(2), -(1)/(2), -(ⅈ)/(2)}, {(1)/(2), -(1)/(2), (1)/(2), -(1)/(2)}, {(1)/(2), -(ⅈ)/(2), -(1)/(2), (ⅈ)/(2)}})
```

A large Fourier matrix:

```wolfram
FourierMatrix[800]
(* Output *)
FourierMatrix[...]
```

### Scope

The real and imaginary parts of the Fourier's basis sequences of length 128:

```wolfram
f=FourierMatrix[128];
{ArrayPlot[Re[f]],ArrayPlot[Im[f]]}
(* Output *)
{![image](img/image_001.png),![image](img/image_002.png)}
```

Construct a structured Fourier matrix using the option setting [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->"Structured":

```wolfram
ℱ = FourierMatrix[512, TargetStructure->"Structured"]
(* Output *)
FourierMatrix[...]
```

The structured representation saves a significant amount of memory for larger matrices:

```wolfram
{ByteCount[ℱ], ByteCount[Normal[ℱ]]}
(* Output *)
{488,119943424}
```

### Options

#### FourierParameters

The default definition of the Fourier matrix:

```wolfram
FourierMatrix[4,FourierParameters->{0,1}]//MatrixForm
(* Output *)
({{(1)/(2), (1)/(2), (1)/(2), (1)/(2)}, {(1)/(2), (ⅈ)/(2), -(1)/(2), -(ⅈ)/(2)}, {(1)/(2), -(1)/(2), (1)/(2), -(1)/(2)}, {(1)/(2), -(ⅈ)/(2), -(1)/(2), (ⅈ)/(2)}})
```

Use the definition of the Fourier matrix used in signal processing:

```wolfram
FourierMatrix[4,FourierParameters->{1,-1}]//MatrixForm
(* Output *)
({{1, 1, 1, 1}, {1, -ⅈ, -1, ⅈ}, {1, -1, 1, -1}, {1, ⅈ, -1, -ⅈ}})
```

Use the definition of the Fourier matrix used in data analysis:

```wolfram
FourierMatrix[4,FourierParameters->{-1,1}]//MatrixForm
(* Output *)
({{(1)/(4), (1)/(4), (1)/(4), (1)/(4)}, {(1)/(4), (ⅈ)/(4), -(1)/(4), -(ⅈ)/(4)}, {(1)/(4), -(1)/(4), (1)/(4), -(1)/(4)}, {(1)/(4), -(ⅈ)/(4), -(1)/(4), (ⅈ)/(4)}})
```

#### TargetStructure

Return the Fourier matrix as a dense matrix:

```wolfram
FourierMatrix[4,TargetStructure->"Dense"]
(* Output *)
{{(1)/(2),(1)/(2),(1)/(2),(1)/(2)},{(1)/(2),(ⅈ)/(2),-(1)/(2),-(ⅈ)/(2)},{(1)/(2),-(1)/(2),(1)/(2),-(1)/(2)},{(1)/(2),-(ⅈ)/(2),-(1)/(2),(ⅈ)/(2)}}
```

Return the Fourier matrix as a structured array:

```wolfram
FourierMatrix[4,TargetStructure->"Structured"]
(* Output *)
FourierMatrix[...]
```

Return the Fourier matrix as a symmetric matrix:

```wolfram
FourierMatrix[4,TargetStructure->"Symmetric"]
(* Output *)
SymmetricMatrix[...]
```

Return the Fourier matrix as a unitary matrix:

```wolfram
FourierMatrix[4,TargetStructure->"Unitary"]
(* Output *)
UnitaryMatrix[...]
```

#### WorkingPrecision

Use machine precision:

```wolfram
FourierMatrix[3,WorkingPrecision->MachinePrecision]//MatrixForm
(* Output *)
({{0.5773502691896258, 0.5773502691896258, 0.5773502691896258}, {0.5773502691896258+0. ⅈ, -0.2886751345948129+0.5 ⅈ, -0.2886751345948129-0.5 ⅈ}, {0.5773502691896258+0. ⅈ, -0.2886751345948129-0.5 ⅈ, -0.2886751345948129+0.5 ⅈ}})
```

Use arbitrary precision:

```wolfram
FourierMatrix[3,WorkingPrecision->20]//MatrixForm
(* Output *)
({{0.5773502691896257645091487805019574556476017512701268760187, 0.5773502691896257645091487805019574556476017512701268760187, 0.5773502691896257645091487805019574556476017512701268760187}, {0.5773502691896257645091487805019574556476017512701268760187, -0.2886751345948128822545743902509787278238008756350634379752+0.5000000000000000000000000000000000000000000000000000000198 ⅈ, -0.2886751345948128822545743902509787278238008756350634379752-0.5000000000000000000000000000000000000000000000000000000198 ⅈ}, {0.5773502691896257645091487805019574556476017512701268760187, -0.2886751345948128822545743902509787278238008756350634379752-0.5000000000000000000000000000000000000000000000000000000198 ⅈ, -0.2886751345948128822545743902509787278238008756350634379752+0.5000000000000000000000000000000000000000000000000000000198 ⅈ}})
```

### Applications

The efficiency of the fast Fourier transform (FFT) relies on being able to form a larger Fourier matrix from two smaller ones. Generate two small Fourier matrices of sizes `*p*` and `*q*`:

```wolfram
p=2;
MatrixForm[fmatp=FourierMatrix[p]]
(* Output *)
({{(1)/(Sqrt[2]), (1)/(Sqrt[2])}, {(1)/(Sqrt[2]), -(1)/(Sqrt[2])}})
```

```wolfram
q=3;
MatrixForm[fmatq=FourierMatrix[q]]
(* Output *)
({{(1)/(Sqrt[3]), (1)/(Sqrt[3]), (1)/(Sqrt[3])}, {(1)/(Sqrt[3]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[3]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[3])}, {(1)/(Sqrt[3]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[3]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[3])}})
```

The Fourier matrix of size `*p* *q*` can be expressed as a product of four simpler matrices:

```wolfram
fpkp=KroneckerProduct[fmatp,IdentityMatrix[q]];
diag=DiagonalMatrix[Flatten[Table[Exp[(2 π ⅈ j k)/(p q)],{k,0,p-1},{j,0,q-1}]]];
fqkp=BlockDiagonalMatrix[ConstantArray[fmatq,p]];
shuffle=PermutationMatrix[Flatten[Partition[Range[p q],p],{{2,1}}]];
MatrixForm[fmatpq=fpkp.diag.fqkp.shuffle]
(* Output *)
({{(1)/(Sqrt[6]), (1)/(Sqrt[6]), (1)/(Sqrt[6]), (1)/(Sqrt[6]), (1)/(Sqrt[6]), (1)/(Sqrt[6])}, {(1)/(Sqrt[6]), (ℯ^((ⅈ π)/(3)))/(Sqrt[6]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), -(1)/(Sqrt[6]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6]), (ℯ^(-(ⅈ π)/(3)))/(Sqrt[6])}, {(1)/(Sqrt[6]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6]), (1)/(Sqrt[6]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6])}, {(1)/(Sqrt[6]), -(1)/(Sqrt[6]), (1)/(Sqrt[6]), -(1)/(Sqrt[6]), (1)/(Sqrt[6]), -(1)/(Sqrt[6])}, {(1)/(Sqrt[6]), -(ℯ^((ⅈ π)/(3)))/(Sqrt[6]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), (1)/(Sqrt[6]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6]), -(ℯ^(-(ⅈ π)/(3)))/(Sqrt[6])}, {(1)/(Sqrt[6]), -(ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), (ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6]), -(1)/(Sqrt[6]), (ℯ^((2 ⅈ π)/(3)))/(Sqrt[6]), -(ℯ^(-(2 ⅈ π)/(3)))/(Sqrt[6])}})
```

Show that the resulting matrix is equivalent to the result of [FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html):

```wolfram
fmatpq==FourierMatrix[p q]
(* Output *)
True
```

The discrete Fourier transform of a vector can be computed by successively multiplying the factors of the Fourier matrix to the vector:

```wolfram
data=RandomComplex[1+ⅈ,p q];
res=fpkp.(diag.(fqkp.(shuffle.data)))
(* Output *)
{1.1254950907792824+1.558908283183154 ⅈ,0.05113395996477769-0.06080324723732458 ⅈ,0.3294711874505031+0.046545792527651564 ⅈ,-0.7108532475625455-0.03209811101182092 ⅈ,0.03399119127622028-0.2235596615429218 ⅈ,-0.02032086572141159-0.2898522662205433 ⅈ}
```

The result is equivalent to applying [Fourier](https://reference.wolfram.com/language/ref/Fourier.html) to the vector:

```wolfram
Chop[Max[Abs[res-Fourier[data]]]]==0
(* Output *)
True
```

Define a function for constructing a circulant matrix from a vector:

```wolfram
circulant[v_]:=ToeplitzMatrix[v,RotateRight[Reverse[v]]]
```

```wolfram
circulant[Array[C,4]]//MatrixForm
(* Output *)
({{1, 4, 3, 2}, {2, 1, 4, 3}, {3, 2, 1, 4}, {4, 3, 2, 1}})
```

Circulant matrices can be diagonalized by the Fourier matrix:

```wolfram
fm=FourierMatrix[4];
ConjugateTranspose[fm].circulant[Array[C,4]].fm//FullSimplify//MatrixForm
(* Output *)
({{1+2+3+4, 0, 0, 0}, {0, 1-ⅈ 2-3+ⅈ 4, 0, 0}, {0, 0, 1-2+3-4, 0}, {0, 0, 0, 1+ⅈ (2+ⅈ 3-4)}})
```

The diagonal elements of the resulting diagonal matrix are the same as the product of the Fourier matrix and the starting vector, up to a constant scaling factor:

```wolfram
Sqrt[4]Conjugate[fm].Array[C,4]//FullSimplify
(* Output *)
{1+2+3+4,1-ⅈ 2-3+ⅈ 4,1-2+3-4,1+ⅈ (2+ⅈ 3-4)}
```

A Fourier matrix with unit normalization:

```wolfram
fmat[n_Integer]:=FourierMatrix[n,FourierParameters->{1,1}]
```

For even dimensions, the permanent of the matrix is zero:

```wolfram
Table[Permanent[fmat[n]]//Simplify,{n,2,10,2}]
(* Output *)
{0,0,0,0,0}
```

For odd dimensions, the permanent of the matrix is always an integer:

```wolfram
Table[Permanent[fmat[n]]//RootReduce,{n,3,11,2}]
(* Output *)
{-3,-5,-105,81,6765}
```

For an odd prime `*p*>3, the permanent of the `*p*`×`*p*` matrix is congruent to `*p*!`, modulo `*p*^(3)`:

```wolfram
With[{p=13},Mod[{RootReduce[Permanent[fmat[p]]],p!},p^3]]
(* Output *)
{2184,2184}
```

### Properties & Relations

[FourierMatrix](https://reference.wolfram.com/language/ref/FourierMatrix.html) can be represented as a scaled [VandermondeMatrix](https://reference.wolfram.com/language/ref/VandermondeMatrix.html):

```wolfram
{n,a,b}={16,1,-1};
fm=FourierMatrix[n,FourierParameters->{a,b}];
fm==(1)/(n^((1-a)/2))VandermondeMatrix[Exp[(2π ⅈ b)/(n) Range[0,n-1]]]
(* Output *)
True
```

The Fourier transform of a vector is equivalent to the vector multiplied by a Fourier matrix:

```wolfram
data= RandomReal[1,6];
FourierMatrix[6] . data ==Fourier[data]
(* Output *)
True
```

The inverse Fourier transform is equivalent to multiplying by the conjugate transpose:

```wolfram
ConjugateTranspose[FourierMatrix[6]]. data ==InverseFourier[data]
(* Output *)
True
```

[Fourier](https://reference.wolfram.com/language/ref/Fourier.html) is much faster than the matrix-based computation:

```wolfram
data=RandomReal[1,1024];
{AbsoluteTiming[Do[Fourier[data],{10}];],
AbsoluteTiming[fm=FourierMatrix[1024,WorkingPrecision->MachinePrecision];Do[fm.data,{10}];]}
(* Output *)
{{0.0053491,Null},{9.6471255,Null}}
```

## Related Guides ▪Structured Arrays ▪Summation Transforms

## History Introduced in 2012 (9.0) | Updated in 2023 (13.3) ▪ 2024 (14.0)
