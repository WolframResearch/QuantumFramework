# FourierDCTMatrix | [SpanFromLeft]

> [FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[*n*]  — returns an `*n*`×`*n*` discrete cosine transform matrix of type 2.
> [FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[*n*,*m*]  — returns an `*n*`×`*n*` discrete cosine transform matrix of type `*m*`.

## Details and Options

Each entry `*F*_*rs*` of the discrete cosine transform matrix of type `*m*` is computed as:

| 1. | DCT-I | $\sqrt{\frac{2}{n-1}}\begin{pmatrix}
\begin{pmatrix}
\frac{1}{2} & r=1 \\
cos(\frac{\pi }{n-1}(r-1)(s-1)) & r=2,...,n-1 \\
(-1)^{s-1} & r=n
\end{pmatrix}
\end{pmatrix}$ |
| --- | --- | --- |
| 2. | DCT-II | $\frac{1}{\sqrt{n}}cos(\frac{\pi}{n}(r-\frac{1}{2})(s-1))$ |
| 3. | DCT-III | $\frac{1}{\sqrt{n}}\begin{pmatrix}
\begin{pmatrix}
1 & r=1 \\
2 cos(\frac{\pi}{n}(r-1)(s-\frac{1}{2})) & r=2,...,n
\end{pmatrix}
\end{pmatrix}$ |
| 4. | DCT-IV | $\sqrt{\frac{2}{n}} cos(\frac{\pi}{n}(r-\frac{1}{2})(s-\frac{1}{2}))$ |

The discrete cosine transform matrices of types 1, 2, 3 and 4 have inverses of type 1, 3, 2 and 4, respectively.

Rows of the [FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html) are basis sequences of the discrete cosine transform.

The result of [FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[*n*].*list* is equivalent to [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*] when `*list*` has length `*n*`. However, the computation of [FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html)[*list*] is much faster and has less numerical error.

For type 4, the option [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) is supported, which specifies the structure of the returned matrix. Possible settings for [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) include:

[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | automatically choose the representation returned
"Dense" | represent the matrix as a dense matrix
"Hermitian" | represent the matrix as a Hermitian matrix
"Orthogonal" | represent the matrix as an orthogonal matrix
"Symmetric" | represent the matrix as a symmetric matrix
"Unitary" | represent the matrix as a unitary matrix

[FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[…,[TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] is equivalent to [FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[…,[TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->"Dense"].

[FourierDCTMatrix](https://reference.wolfram.com/language/ref/FourierDCTMatrix.html)[…,[WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)->*p*] gives a matrix with entries of precision `*p*`.

## Examples

### Basic Examples

A 4×4 DCT matrix:

```wolfram
FourierDCTMatrix[4]//MatrixForm
(* Output *)
({{(1)/(2), (1)/(2) Cos[(π)/(8)], (1)/(2 Sqrt[2]), (1)/(2) Sin[(π)/(8)]}, {(1)/(2), (1)/(2) Sin[(π)/(8)], -(1)/(2 Sqrt[2]), -(1)/(2) Cos[(π)/(8)]}, {(1)/(2), -(1)/(2) Sin[(π)/(8)], -(1)/(2 Sqrt[2]), (1)/(2) Cos[(π)/(8)]}, {(1)/(2), -(1)/(2) Cos[(π)/(8)], (1)/(2 Sqrt[2]), -(1)/(2) Sin[(π)/(8)]}})
```

### Scope

The discrete cosine transform's basis sequences of length 128:

```wolfram
ArrayPlot[FourierDCTMatrix[128]]
```

![image](img/image_001.png)

### Options

#### TargetStructure

Return the DCT matrix as a dense matrix:

```wolfram
FourierDCTMatrix[4,"IV",TargetStructure->"Dense"]
(* Output *)
{{(Cos[(π)/(16)])/(Sqrt[2]),(Cos[(3 π)/(16)])/(Sqrt[2]),(Sin[(3 π)/(16)])/(Sqrt[2]),(Sin[(π)/(16)])/(Sqrt[2])},{(Cos[(3 π)/(16)])/(Sqrt[2]),-(Sin[(π)/(16)])/(Sqrt[2]),-(Cos[(π)/(16)])/(Sqrt[2]),-(Sin[(3 π)/(16)])/(Sqrt[2])},{(Sin[(3 π)/(16)])/(Sqrt[2]),-(Cos[(π)/(16)])/(Sqrt[2]),(Sin[(π)/(16)])/(Sqrt[2]),(Cos[(3 π)/(16)])/(Sqrt[2])},{(Sin[(π)/(16)])/(Sqrt[2]),-(Sin[(3 π)/(16)])/(Sqrt[2]),(Cos[(3 π)/(16)])/(Sqrt[2]),-(Cos[(π)/(16)])/(Sqrt[2])}}
```

Return the DCT matrix as an orthogonal matrix:

```wolfram
FourierDCTMatrix[4,"IV",TargetStructure->"Orthogonal"]
(* Output *)
OrthogonalMatrix[...]
```

Return the DCT matrix as a symmetric matrix:

```wolfram
FourierDCTMatrix[4,"IV",TargetStructure->"Symmetric"]
(* Output *)
SymmetricMatrix[...]
```

#### WorkingPrecision

Use machine precision:

```wolfram
FourierDCTMatrix[3,WorkingPrecision->MachinePrecision]//MatrixForm
(* Output *)
({{0.5773502691896258, 0.5000000000000001, 0.288675134594813}, {0.5773502691896258, -8.326672684688674×10^-17, -0.5773502691896257}, {0.5773502691896258, -0.5, 0.2886751345948128}})
```

Use arbitrary precision:

```wolfram
FourierDCTMatrix[3,WorkingPrecision->20]//MatrixForm
(* Output *)
({{0.5773502691896257645091487805019574556476017512701268760187, 0.5000000000000000000000000000000000000000000000000000000198, 0.2886751345948128822545743902509787278238008756350634379752}, {0.5773502691896257645091487805019574556476017512701268760187, 0`20.042440750329686, -0.5773502691896257645091487805019574556476017512701268760187}, {0.5773502691896257645091487805019574556476017512701268760187, -0.4999999999999999999999999999999999999998297907111648145392, 0.2886751345948128822545743902509787278238008756350634380093}})
```

### Applications

Define 2D discrete cosine transform of size 8×8 using matrix formulation:

```wolfram
DCT[x_?MatrixQ]:=Transpose[FourierDCTMatrix[8]].x.FourierDCTMatrix[8]
IDCT[x_?MatrixQ]:=Transpose[FourierDCTMatrix[8,3]].x.FourierDCTMatrix[8, 3]
img =![image](img/image_002.png);
```

Simplified JPEG compression algorithm:

```wolfram
jpeg=ImageAssemble[Map[Image[Threshold[DCT[ImageData[#]],{"LargestValues",8}]]&,ImagePartition[ImageSubtract[img, 0.5], 8],{2}]]
```

*([Graphics])*

Compare the original and compressed images:

```wolfram
Graphicsimg,ImageAdd[ImageAssemble[Map[Image[IDCT[ImageData[#]]]&,ImagePartition[jpeg, 8],{2],0.5]}, ImageSize->250]
```

![image](img/image_003.png)

### Properties & Relations

A DCT matrix multiplied by a vector is equivalent to the discrete cosine transform of that vector:

```wolfram
data= RandomReal[1,6];
data . FourierDCTMatrix[6]==FourierDCT[data,2]
(* Output *)
True
```

[FourierDCT](https://reference.wolfram.com/language/ref/FourierDCT.html) is much faster than the matrix-based computation:

```wolfram
data=RandomReal[1,1024];
{AbsoluteTiming[Do[FourierDCT[data],{10}];],
AbsoluteTiming[Do[FourierDCTMatrix[1024,WorkingPrecision->MachinePrecision].data,{10}];]}
(* Output *)
{{0.00046200000000000000831626434383281321,Null},{0.57981300000000002281552724525681696832,Null}}
```

A discrete cosine transform matrix of type 1 is its own inverse:

```wolfram
FourierDCTMatrix[4,1]. FourierDCTMatrix[4,1]==IdentityMatrix[4]
(* Output *)
True
```

A discrete cosine transform matrix of type 3 is an inverse of the type 2 matrix:

```wolfram
FourierDCTMatrix[4,2]. FourierDCTMatrix[4,3]==IdentityMatrix[4]
(* Output *)
True
```

A discrete cosine transform matrix of type 4 is its own inverse:

```wolfram
FullSimplify[FourierDCTMatrix[4,4]. FourierDCTMatrix[4,4]]==IdentityMatrix[4]
(* Output *)
True
```

## Related Guides ▪Summation Transforms ▪Structured Arrays

## History Introduced in 2012 (9.0) | Updated in 2024 (14.0)
