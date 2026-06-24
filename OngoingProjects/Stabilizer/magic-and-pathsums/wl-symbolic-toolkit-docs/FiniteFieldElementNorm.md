# FiniteFieldElementNorm | [SpanFromLeft]

> [FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*]  — gives the absolute norm of the finite field element `*a*`.
> [FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*,*k*]  — gives the norm of `*a*` relative to the $p^{k}$-element subfield of the ambient field of `*a*`.
> [FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*,*emb*]  — gives the norm of `*a*` relative to the finite field embedding `*emb*`.

## Details

For a finite field $\mathcal{F}$ with characteristic `*p*` and extension degree `*d*` over $\mathbb{Z}_{p}$, the absolute norm of `*a*` is given by $N_{\mathcal{F}}(a)=a \cdot a^{p}\cdot \cdots \cdot a^{p^{d-1}}$. $N_{\mathcal{F}}$ is a mapping from $\mathcal{F}$ to $\mathbb{Z}_{p}$ and $N_{\mathcal{F}}(a b)=N_{\mathcal{F}}(a) N_{\mathcal{F}}(b)$.

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $N_{\mathcal{F}}(a)=(-1)^{d}c_{0}^{d/n}$.

[FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*] gives an integer between $0$ and $p^{d}-1$.

For a finite field $\mathcal{F}$ with characteristic `*p*` and extension degree `*d*` over $\mathbb{Z}_{p}$, the norm of `*a*` relative to the $p^{k}$-element subfield $\mathcal{K}$ of $\mathcal{F}$ is given by $N_{\mathcal{F}/\mathcal{K}}(a)=a \cdot a^{q}\cdot \cdots \cdot a^{q^{d/k-1}}$, where $q=p^{k}$. $N_{\mathcal{F}/\mathcal{K}}$ is a mapping from $\mathcal{F}$ to $\mathcal{K}$ and $N_{\mathcal{F}/\mathcal{K}}(a b)=N_{\mathcal{F}/\mathcal{K}}(a) N_{\mathcal{F}/\mathcal{K}}(b)$. `*k*` needs to be a divisor of `*d*`.

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*,*k*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $N_{\mathcal{F}/K}(a)=(-1)^{d/k}c_{0}^{d/(k n)}$.

[FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*,*k*] gives an element of $\mathcal{K}$.

If `*emb*=[FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html)[*e*_1->*e*_2]`, then [FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html)[*a*,*emb*] effectively gives `*emb*["Projection"][FiniteFieldElementNorm[*a*,*k*]]`, where `*a*` belongs to the ambient field of `*e*_2 and `*k*` is the extension degree of the ambient field of `*e*_1

## Examples

### Basic Examples

Represent a finite field with characteristic $17$ and extension degree $6$:

```wolfram
ℱ=FiniteField[17,6]
(* Output *)
FiniteField[...]
```

Find the absolute norm of an element of the field:

```wolfram
FiniteFieldElementNorm[ℱ[123]]
(* Output *)
3
```

Find the norm relative to the $17^{2}$-element subfield:

```wolfram
FiniteFieldElementNorm[ℱ[123],2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[17, 3, +, 3, #, +, 10, #, ^, 2, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 3031306], index -> 8583881, shortIndex -> 8583881, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

### Scope

Find the absolute norm of a finite field element:

```wolfram
ℱ=FiniteField[127,12];
FiniteFieldElementNorm[ℱ[1234]]
(* Output *)
119
```

The absolute norm given as a finite field element:

```wolfram
FiniteFieldElementNorm[ℱ[1234],1]
(* Output *)
![image](img/image_001.png)
```

The norm relative to the $127^{3}$-element subfield:

```wolfram
FiniteFieldElementNorm[ℱ[1234],3]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[127, 3, +, 8, #, +, 99, #, ^, 2, +, 15, #, ^, 3, +, 97, #, ^, 4, +, 33, #, ^, 5, +, 25, #, ^, 6, +, 119, #, ^, 7, +, #, ^, 12, &, Polynomial], 11043120644152403225981999], index -> 13745438272967982333049823, shortIndex -> "13745<<16>>49823", indexShortened -> True, characteristic -> 127, shortCharacteristic -> 127, extensionDegree -> 12, field -> FiniteField[...], fieldDisplayed -> False|>
```

Compute the norm relative to a field embedding:

```wolfram
{𝒦,ℱ}={FiniteField[73,2],FiniteField[73,8]};
ℰ=FiniteFieldEmbedding[𝒦,ℱ];
FiniteFieldElementNorm[ℱ[1234],ℰ]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[73, 5, +, 70, #, +, #, ^, 2, &, Polynomial], 1636], index -> 2644, shortIndex -> 2644, indexShortened -> True, characteristic -> 73, shortCharacteristic -> 73, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>
```

The result is equivalent to computing the norm relative to $\mathcal{E}(\mathcal{K})$ and projecting it to $\mathcal{K}$:

```wolfram
FiniteFieldElementNorm[ℱ[1234],2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[73, 5, +, 18, #, +, 39, #, ^, 2, +, 53, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 8, &, Polynomial], 531143263267513], index -> 41000055982918, shortIndex -> 41000055982918, indexShortened -> True, characteristic -> 73, shortCharacteristic -> 73, extensionDegree -> 8, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
ℰ["Projection"][%]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[73, 5, +, 70, #, +, #, ^, 2, &, Polynomial], 1636], index -> 2644, shortIndex -> 2644, indexShortened -> True, characteristic -> 73, shortCharacteristic -> 73, extensionDegree -> 2, field -> FiniteField[...], fieldDisplayed -> False|>
```

### Applications

Define $\mathbb{Z}_{p}$-linear mappings $m_{b}:\mathcal{F}∋a \longrightarrow b a \in \mathcal{F}$:

```wolfram
ℱ=FiniteField[59,3];
b=ℱ[123];
m[b_][a_]:=a b
```

$N_{\mathcal{F}}(b)$ computes the determinant of $m_{b}$:

```wolfram
FiniteFieldElementNorm[b]
(* Output *)
5
```

Compute the determinant manually:

```wolfram
Mod[Det[{PadRight[ m[b][ℱ[{1,0,0}]]["Coefficients"],3],PadRight[m[b][ℱ[{0,1,0}]]["Coefficients"],3],PadRight[m[b][ℱ[{0,0,1}]]["Coefficients"],3]}],59]
(* Output *)
5
```

### Properties & Relations

$N_{\mathcal{F}}$ is a mapping from $\mathcal{F}$ to $\mathbb{Z}_{p}$ which preserves multiplication:

```wolfram
ℱ=FiniteField[79,3];
{a,b}={ℱ[123], ℱ[456]};
FiniteFieldElementNorm[a b]
(* Output *)
9
```

```wolfram
Mod[FiniteFieldElementNorm[a] FiniteFieldElementNorm[b],79]
(* Output *)
9
```

The absolute norm of `*a*` is equal to the product of all conjugates of `*a*`:

```wolfram
ℱ=FiniteField[17,4];
a=ℱ[123];
FiniteFieldElementNorm[a]
(* Output *)
9
```

Use [FrobeniusAutomorphism](https://reference.wolfram.com/language/ref/FrobeniusAutomorphism.html) to compute the conjugates of `*a*`:

```wolfram
Product[FrobeniusAutomorphism[a,k],{k,0,3}]
(* Output *)
![image](img/image_003.png)
```

The absolute norm of $a^{p}$ is equal to the absolute norm of $a$:

```wolfram
ℱ=FiniteField[101,3];
a=ℱ[123];
FiniteFieldElementNorm[a]
(* Output *)
10
```

```wolfram
FiniteFieldElementNorm[a^101]
(* Output *)
10
```

If $\mathcal{K}$ is the $p^{k}$-element subfield of $\mathcal{F}$, then $N_{\mathcal{F}/\mathcal{K}}$ is a mapping from $\mathcal{F}$ to $\mathcal{K}$, which preserves multiplication:

```wolfram
ℱ=FiniteField[11,6];
{a,b}={ℱ[123], ℱ[456]};
{c,d}={FiniteFieldElementNorm[a,3],FiniteFieldElementNorm[b,3]}
(* Output *)
![image](img/image_005.png)
```

Use [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html) to show that `*c*` and `*d*` belong to the $11^{3}$-element subfield $\mathcal{K}$ of $\mathcal{F}$:

```wolfram
Exponent[{MinimalPolynomial[c,x],MinimalPolynomial[d,x]},x]
(* Output *)
{3,3}
```

This illustrates the multiplication-preserving property of $N_{\mathcal{F}/\mathcal{K}}$:

```wolfram
FiniteFieldElementNorm[a b,3]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[11, 2, +, 7, #, +, 6, #, ^, 2, +, 4, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 6, &, Polynomial], 322225], index -> 837466, shortIndex -> 837466, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
c d
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[11, 2, +, 7, #, +, 6, #, ^, 2, +, 4, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 6, &, Polynomial], 322225], index -> 837466, shortIndex -> 837466, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

Construct field embeddings such that $\mathcal{E}_{21}(a)=\mathcal{E}_{2}(\mathcal{E}_{1}(a))$:

```wolfram
𝒦=FiniteField[71,3];
ℱ=FiniteField[71,6];
𝒢=FiniteField[71,12];
ℰ_1=FiniteFieldEmbedding[𝒦,ℱ];
ℰ_2=FiniteFieldEmbedding[ℱ,𝒢];
ℰ_21=ℰ_2@*ℰ_1;
ℰ_21[𝒦[123]]==ℰ_2[ℰ_1[𝒦[123]]]
(* Output *)
True
```

[FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html) satisfies a transitivity property:

```wolfram
FiniteFieldElementNorm[𝒢[1234],ℰ_21]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 4, #, +, #, ^, 3, &, Polynomial], 4434], index -> 2458, shortIndex -> 2458, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
FiniteFieldElementNorm[FiniteFieldElementNorm[𝒢[1234],ℰ_2],ℰ_1]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 4, #, +, #, ^, 3, &, Polynomial], 4434], index -> 2458, shortIndex -> 2458, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $N_{\mathcal{F}}(a)=(-1)^{d}c_{0}^{d/n}$:

```wolfram
ℱ=FiniteField[149,5];
a=ℱ[1234];
FiniteFieldElementNorm[a]
(* Output *)
53
```

```wolfram
MinimalPolynomial[a,x]
(* Output *)
96+46 x+97 x^2+58 x^3+88 x^4+x^5
```

```wolfram
Mod[-Coefficient[%,x,0],149]
(* Output *)
53
```

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*,*k*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $N_{\mathcal{F}/K}(a)=(-1)^{d/k}c_{0}^{d/(k n)}$:

```wolfram
ℱ=FiniteField[43,6];
a=ℱ[1234];
FiniteFieldElementNorm[a,2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 37261938517], index -> 2519295088, shortIndex -> 2519295088, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
MinimalPolynomial[a,x,2]
(* Output *)
![image](img/image_007.png)
```

```wolfram
-Coefficient[%,x,0]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 37261938517], index -> 2519295088, shortIndex -> 2519295088, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

## Related Guides ▪Finite Fields

## History Introduced in 2023 (13.3)
