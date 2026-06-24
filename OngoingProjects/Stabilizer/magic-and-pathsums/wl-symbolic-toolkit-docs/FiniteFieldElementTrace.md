# FiniteFieldElementTrace | [SpanFromLeft]

> [FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*]  — gives the absolute trace of the finite field element `*a*`.
> [FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*,*k*]  — gives the trace of `*a*` relative to the $p^{k}$-element subfield of the ambient field of `*a*`.
> [FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*,*emb*]  — gives the trace of `*a*` relative to the finite field embedding `*emb*`.

## Details

For a finite field $\mathcal{F}$ with characteristic `*p*` and extension degree `*d*` over $\mathbb{Z}_{p}$, the absolute trace of `*a*` is given by $Tr_{\mathcal{F}}(a)=a+a^{p}+\cdots+a^{p^{d-1}}$. $Tr_{\mathcal{F}}$ is a $\mathbb{Z}_{p}$-linear mapping from $\mathcal{F}$ to $\mathbb{Z}_{p}$.

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $Tr_{\mathcal{F}}(a)=-\frac{d }{n}c_{n-1}$.

[FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*] gives an integer between $0$ and $p^{d}-1$.

For a finite field $\mathcal{F}$ with characteristic `*p*` and extension degree `*d*` over $\mathbb{Z}_{p}$, the trace of `*a*` relative to the $p^{k}$-element subfield $\mathcal{K}$ of $\mathcal{F}$ is given by $Tr_{\mathcal{F}/\mathcal{K}}(a)=a+a^{q}+\cdots+a^{q^{d/k-1}}$, where $q=p^{k}$. $Tr_{\mathcal{F}/\mathcal{K}}$ is a $\mathcal{K}$-linear mapping from $\mathcal{F}$ to $\mathcal{K}$. `*k*` needs to be a divisor of `*d*`.

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*,*k*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $Tr_{\mathcal{F}/K}(a)=-\frac{d }{k n}c_{n-1}$.

[FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*,*k*] gives an element of $\mathcal{K}$.

If `*emb*=[FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html)[*e*_1->*e*_2]`, then [FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html)[*a*,*emb*] effectively gives `*emb*["Projection"][FiniteFieldElementTrace[*a*,*k*]]`, where `*a*` belongs to the ambient field of `*e*_2 and `*k*` is the extension degree of the ambient field of `*e*_1

## Examples

### Basic Examples

Represent a finite field with characteristic $17$ and extension degree $6$:

```wolfram
ℱ=FiniteField[17,6]
(* Output *)
FiniteField[...]
```

Find the absolute trace of an element of the field:

```wolfram
FiniteFieldElementTrace[ℱ[123]]
(* Output *)
7
```

Find the trace relative to the $17^{2}$-element subfield:

```wolfram
FiniteFieldElementTrace[ℱ[123],2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[17, 3, +, 3, #, +, 10, #, ^, 2, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 120516010], index -> 14278635, shortIndex -> 14278635, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

### Scope

Find the absolute trace of a finite field element:

```wolfram
ℱ=FiniteField[127,12];
FiniteFieldElementTrace[ℱ[1234]]
(* Output *)
76
```

The absolute trace given as a finite field element:

```wolfram
FiniteFieldElementTrace[ℱ[1234],1]
(* Output *)
![image](img/image_001.png)
```

The trace relative to the $127^{3}$-element subfield:

```wolfram
FiniteFieldElementTrace[ℱ[1234],3]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[127, 3, +, 8, #, +, 99, #, ^, 2, +, 15, #, ^, 3, +, 97, #, ^, 4, +, 33, #, ^, 5, +, 25, #, ^, 6, +, 119, #, ^, 7, +, #, ^, 12, &, Polynomial], 4310440611047703354113448], index -> 6691200733984981989115878, shortIndex -> "66912<<15>>15878", indexShortened -> True, characteristic -> 127, shortCharacteristic -> 127, extensionDegree -> 12, field -> FiniteField[...], fieldDisplayed -> False|>
```

Compute the trace relative to a field embedding:

```wolfram
{𝒦,ℱ}={FiniteField[73,2],FiniteField[73,8]};
ℰ=FiniteFieldEmbedding[𝒦,ℱ];
FiniteFieldElementTrace[ℱ[1234],ℰ]
(* Output *)
![image](img/image_003.png)
```

The result is equivalent to computing the trace relative to $\mathcal{E}(\mathcal{K})$ and projecting it to $\mathcal{K}$:

```wolfram
FiniteFieldElementTrace[ℱ[1234],2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[73, 5, +, 18, #, +, 39, #, ^, 2, +, 53, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 8, &, Polynomial], 626230474162270], index -> 776660870520844, shortIndex -> 776660870520844, indexShortened -> True, characteristic -> 73, shortCharacteristic -> 73, extensionDegree -> 8, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
ℰ["Projection"][%]
(* Output *)
![image](img/image_005.png)
```

### Applications

Define $\mathbb{Z}_{p}$-linear mappings $\mathcal{L}_{b}:\mathcal{F}∋a \longrightarrow Tr_{\mathcal{F}}(b a)\in \mathbb{Z}_{p}$. Every $\mathbb{Z}_{p}$-linear mapping from $\mathcal{F}$ to $\mathbb{Z}_{p}$ has this form:

```wolfram
ℱ=FiniteField[59,3];
b=ℱ[123];
ℒ[b_][a_]:=FiniteFieldElementTrace[a b]
```

Illustrate linearity of $\mathcal{L}_{b}$:

```wolfram
ℒ[b][15ℱ[456]+27 ℱ[789]]
(* Output *)
25
```

```wolfram
Mod[15 ℒ[b][ℱ[456]]+27ℒ[b][ℱ[789]],59]
(* Output *)
25
```

### Properties & Relations

$Tr_{\mathcal{F}}$ is a $\mathbb{Z}_{p}$-linear mapping from $\mathcal{F}$ to $\mathbb{Z}_{p}$:

```wolfram
ℱ=FiniteField[29,3];
{a,b}={ℱ[123], ℱ[456]};
FiniteFieldElementTrace[2a+3b]
(* Output *)
28
```

```wolfram
Mod[2FiniteFieldElementTrace[a]+3FiniteFieldElementTrace[b],29]
(* Output *)
28
```

The absolute trace of `*a*` is equal to the sum of all conjugates of `*a*`:

```wolfram
ℱ=FiniteField[17,4];
a=ℱ[123];
FiniteFieldElementTrace[a]
(* Output *)
16
```

Use [FrobeniusAutomorphism](https://reference.wolfram.com/language/ref/FrobeniusAutomorphism.html) to compute the conjugates of `*a*`:

```wolfram
Sum[FrobeniusAutomorphism[a,k],{k,0,3}]
(* Output *)
![image](img/image_007.png)
```

The absolute trace of $a^{p}$ is equal to the absolute trace of $a$:

```wolfram
ℱ=FiniteField[101,3];
a=ℱ[123];
FiniteFieldElementTrace[a]
(* Output *)
66
```

```wolfram
FiniteFieldElementTrace[a^101]
(* Output *)
66
```

If $\mathcal{K}$ is the $p^{k}$-element subfield of $\mathcal{F}$, then $Tr_{\mathcal{F}/\mathcal{K}}$ is a $\mathcal{K}$-linear mapping from $\mathcal{F}$ to $\mathcal{K}$:

```wolfram
ℱ=FiniteField[11,6];
{a,b}={ℱ[123], ℱ[456]};
{c,d}={FiniteFieldElementTrace[a,3],FiniteFieldElementTrace[b,3]}
(* Output *)
![image](img/image_009.png)
```

Use [FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html) to embed an $11^{3}$-element field $\mathcal{G}$ in $\mathcal{F}$:

```wolfram
𝒢=FiniteField[11,3];
ℰ=FiniteFieldEmbedding[𝒢,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[11, 9, +, 2, #, +, #, ^, 3, &, Polynomial], 01], index -> 11, shortIndex -> 11, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[11, 2, +, 7, #, +, 6, #, ^, 2, +, 4, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 6, &, Polynomial], 906331], index -> 209702, shortIndex -> 209702, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Since $\mathcal{E}(\mathcal{G})=\mathcal{K}$, this shows that `*c*` and `*d*` belong to $\mathcal{K}$:

```wolfram
{ℰ[ℰ["Projection"][c]]==c,ℰ[ℰ["Projection"][d]]==d}
(* Output *)
{True,True}
```

This illustrates $\mathcal{K}$-linearity of $Tr_{\mathcal{F}/\mathcal{K}}$:

```wolfram
FiniteFieldElementTrace[ℰ[𝒢[321]]a+ℰ[𝒢[654]]b,3]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[11, 2, +, 7, #, +, 6, #, ^, 2, +, 4, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 6, &, Polynomial], 45410105], index -> 965518, shortIndex -> 965518, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
ℰ[𝒢[321]]c+ℰ[𝒢[654]]d
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[11, 2, +, 7, #, +, 6, #, ^, 2, +, 4, #, ^, 3, +, 3, #, ^, 4, +, #, ^, 6, &, Polynomial], 45410105], index -> 965518, shortIndex -> 965518, indexShortened -> True, characteristic -> 11, shortCharacteristic -> 11, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
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

[FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html) satisfies a transitivity property:

```wolfram
FiniteFieldElementTrace[𝒢[1234],ℰ_21]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 4, #, +, #, ^, 3, &, Polynomial], 24484], index -> 23596, shortIndex -> 23596, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
FiniteFieldElementTrace[FiniteFieldElementTrace[𝒢[1234],ℰ_2],ℰ_1]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 4, #, +, #, ^, 3, &, Polynomial], 24484], index -> 23596, shortIndex -> 23596, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $Tr_{\mathcal{F}}(a)=-\frac{d }{n}c_{n-1}$:

```wolfram
ℱ=FiniteField[149,5];
a=ℱ[1234];
FiniteFieldElementTrace[a]
(* Output *)
61
```

```wolfram
MinimalPolynomial[a,x]
(* Output *)
96+46 x+97 x^2+58 x^3+88 x^4+x^5
```

```wolfram
Mod[-Coefficient[%,x,4],149]
(* Output *)
61
```

If [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*a*,*x*,*k*]==*x*^(*n*)+*c*_*n*-1*x*^(*n*-1)+⋯+*c*_0, then $Tr_{\mathcal{F}/K}(a)=-\frac{d }{k n}c_{n-1}$:

```wolfram
ℱ=FiniteField[43,6];
a=ℱ[1234];
FiniteFieldElementTrace[a,2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 2229363714], index -> 2185097697, shortIndex -> 2185097697, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
MinimalPolynomial[a,x,2]
(* Output *)
![image](img/image_011.png)
```

```wolfram
-Coefficient[%,x,2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 2229363714], index -> 2185097697, shortIndex -> 2185097697, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

## Related Guides ▪Finite Fields

## History Introduced in 2023 (13.3)
