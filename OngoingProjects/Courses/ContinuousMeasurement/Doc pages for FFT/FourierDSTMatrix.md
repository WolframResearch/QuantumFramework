# FourierDSTMatrix | [SpanFromLeft]

> [FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[*n*]  — returns an `*n*`×`*n*` discrete sine transform matrix of type 2.
> [FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[*n*,*m*]  — returns an `*n*`×`*n*` discrete sine transform matrix of type `*m*`.

## Details and Options

Each entry `*F*_*rs*` of the discrete sine transform matrix of type `*m*` is computed as:

| 1. | DST-I | $\sqrt{\frac{2}{n+1}} sin(\frac{\pi}{n+1} r s)$ |
| --- | --- | --- |
| 2. | DST-II | $\frac{1}{\sqrt{n}}sin(\frac{\pi}{n} (r-\frac{1}{2}) s)$ |
| 3. | DST-III | $\frac{1}{\sqrt{n}}\begin{pmatrix}
\begin{pmatrix}
2 sin(\frac{\pi}{n} r (s-\frac{1}{2})), & r=1,...,n-1 \\
(-1)^{s-1}, & r=n
\end{pmatrix}
\end{pmatrix}$ |
| 4. | DST-IV | $\sqrt{\frac{2}{n}} sin(\frac{\pi}{n} (r-\frac{1}{2}) (s-\frac{1}{2}))$ |

The discrete sine transform matrices of types 1, 2, 3 and 4 have inverses of type 1, 3, 2 and 4, respectively.

Rows of the [FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html) are basis sequences of the discrete sine transform.

The result of [FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[*n*].*list* is equivalent to [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*] when `*list*` has length `*n*`. However, the computation of [FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html)[*list*] is much faster and has less numerical error.

For types 1 and 4, the option [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) is supported, which specifies the structure of the returned matrix. Possible settings for [TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html) include:

[Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | automatically choose the representation returned
"Dense" | represent the matrix as a dense matrix
"Hermitian" | represent the matrix as a Hermitian matrix
"Orthogonal" | represent the matrix as an orthogonal matrix
"Symmetric" | represent the matrix as a symmetric matrix
"Unitary" | represent the matrix as a unitary matrix

[FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[…,[TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->[Automatic](https://reference.wolfram.com/language/ref/Automatic.html)] is equivalent to [FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[…,[TargetStructure](https://reference.wolfram.com/language/ref/TargetStructure.html)->"Dense"].

[FourierDSTMatrix](https://reference.wolfram.com/language/ref/FourierDSTMatrix.html)[…,[WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)->*p*] gives a matrix with entries of precision `*p*`.

## Examples

### Basic Examples

A 4×4 DST matrix:

```wolfram
FourierDSTMatrix[4]//MatrixForm
(* Output *)
({{(1)/(2) Sin[(π)/(8)], (1)/(2 Sqrt[2]), (1)/(2) Cos[(π)/(8)], (1)/(2)}, {(1)/(2) Cos[(π)/(8)], (1)/(2 Sqrt[2]), -(1)/(2) Sin[(π)/(8)], -(1)/(2)}, {(1)/(2) Cos[(π)/(8)], -(1)/(2 Sqrt[2]), -(1)/(2) Sin[(π)/(8)], (1)/(2)}, {(1)/(2) Sin[(π)/(8)], -(1)/(2 Sqrt[2]), (1)/(2) Cos[(π)/(8)], -(1)/(2)}})
```

### Scope

The discrete sine transform's basis sequences of length 128:

```wolfram
ArrayPlot[FourierDSTMatrix[128]]
```

![image](img/image_001.png)

### Options

#### TargetStructure

Return the DST matrix as a dense matrix:

```wolfram
FourierDSTMatrix[4,"I",TargetStructure->"Dense"]
(* Output *)
{{Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))]},{Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))]},{Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))]},{Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],Sqrt[(2)/(5) ((5)/(8)+(Sqrt[5])/(8))],-Sqrt[(2)/(5) ((5)/(8)-(Sqrt[5])/(8))]}}
```

Return the DST matrix as an orthogonal matrix:

```wolfram
FourierDSTMatrix[4,"I",TargetStructure->"Orthogonal"]
(* Output *)
OrthogonalMatrix[...]
```

Return the DST matrix as a symmetric matrix:

```wolfram
FourierDSTMatrix[4,"I",TargetStructure->"Symmetric"]
(* Output *)
SymmetricMatrix[...]
```

#### WorkingPrecision

Use machine precision:

```wolfram
FourierDSTMatrix[3,WorkingPrecision->MachinePrecision]//MatrixForm
(* Output *)
({{0.28867513459481287, 0.5, 0.5773502691896258}, {0.5773502691896258, -5.551115123125783×10^-17, -0.5773502691896258}, {0.288675134594813, -0.5, 0.5773502691896258}})
```

Use arbitrary precision:

```wolfram
FourierDSTMatrix[3,WorkingPrecision->20]//MatrixForm
(* Output *)
({{0.2886751345948128822545743902509787278238008756350634379752, 0.5000000000000000000000000000000000000000000000000000000198, 0.5773502691896257645091487805019574556476017512701268760187}, {0.5773502691896257645091487805019574556476017512701268760187, 0`19.741410754665704, -0.5773502691896257645091487805019574556476017512701268760187}, {0.2886751345948128822545743902509787278240956867712461452763, -0.4999999999999999999999999999999999999995250378728584380749, 0.5773502691896257645091487805019574556476017512701268760187}})
```

### Applications

A tridiagonal Toeplitz matrix:

```wolfram
n=4;
td=Normal[SparseArray[{Band[{1,2}]->a,Band[{1,1}]->b,Band[{2,1}]->c},{n,n}]]
(* Output *)
{{b,a,0,0},{c,b,a,0},{0,c,b,a},{0,0,c,b}}
```

```wolfram
MatrixForm[%]
(* Output *)
({{b, a, 0, 0}, {c, b, a, 0}, {0, c, b, a}, {0, 0, c, b}})
```

The matrix of its eigenvectors can be expressed as a diagonal rescaling of the discrete sine transform matrix of type 1:

```wolfram
evecs=DiagonalMatrix[(Sqrt[c]/Sqrt[a])^Range[n]].FourierDSTMatrix[n,"I"];
```

Diagonalize the tridiagonal matrix:

```wolfram
Inverse[evecs].td.evecs//FullSimplify
(* Output *)
{{b+(1)/(2) (1+Sqrt[5]) Sqrt[a] Sqrt[c],0,0,0},{0,b+(1)/(2) (-1+Sqrt[5]) Sqrt[a] Sqrt[c],0,0},{0,0,b-(1)/(2) (-1+Sqrt[5]) Sqrt[a] Sqrt[c],0},{0,0,0,b-(1)/(2) (1+Sqrt[5]) Sqrt[a] Sqrt[c]}}
```

### Properties & Relations

A DST matrix multiplied by a vector is equivalent to the discrete sine transform of that vector:

```wolfram
data= RandomReal[1,6];
data . FourierDSTMatrix[6]==FourierDST[data,2]
(* Output *)
True
```

[FourierDST](https://reference.wolfram.com/language/ref/FourierDST.html) is much faster than the matrix-based computation:

```wolfram
data=RandomReal[1,1024];
{AbsoluteTiming[Do[FourierDST[data],{10}];],
AbsoluteTiming[Do[FourierDSTMatrix[1024,WorkingPrecision->MachinePrecision].data,{10}];]}
(* Output *)
{{0.00044400000000000000447211712106820869,Null},{0.54677100000000000701305680195218883455,Null}}
```

A discrete sine transform matrix of type 1 is its own inverse:

```wolfram
FourierDSTMatrix[4,1]. FourierDSTMatrix[4,1]==IdentityMatrix[4]
(* Output *)
True
```

A discrete sine transform matrix of type 3 is an inverse of the type 2 matrix:

```wolfram
FourierDSTMatrix[4,2]. FourierDSTMatrix[4,3]==IdentityMatrix[4]
(* Output *)
True
```

A discrete sine transform matrix of type 4 is its own inverse:

```wolfram
FullSimplify[FourierDSTMatrix[4,4]. FourierDSTMatrix[4,4]]==IdentityMatrix[4]
(* Output *)
True
```

## Related Guides ▪Summation Transforms

## History Introduced in 2012 (9.0) | Updated in 2024 (14.0)
