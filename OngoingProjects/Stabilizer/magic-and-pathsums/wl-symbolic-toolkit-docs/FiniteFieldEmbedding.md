# FiniteFieldEmbedding | [SpanFromLeft]

> [FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html)[*ff*_1,*ff*_2]  — gives an embedding of the finite field `*ff*_1 in the finite field `*ff*_2
> [FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html)[*e*_1->*e*_2]  — represents the embedding of the ambient field of `*e*_1 in the ambient field of `*e*_2, which maps `*e*_1 to `*e*_2

## Details

Finite field embeddings are also known as Galois field embeddings or finite field monomorphisms.

Finite field embeddings are typically used to identify one finite field with a subfield of another.

If `ℰ=FiniteFieldEmbedding[*e*_1->*e*_2]`, where `*e*_1∈*ff*_1 and `*e*_2∈*ff*_2, then $\mathcal{E}$ maps `*ff*_1 into `*ff*_2, $\mathcal{E}(a+b)=\mathcal{E}(a)+\mathcal{E}(b)$, and $\mathcal{E}(a b)=\mathcal{E}(a) \mathcal{E}(b)$ for all $\mathit{a}\mathit{,}\mathit{b}\in \mathit{ff}_{1}$.

A finite field `*ff*_1 can be embedded in `*ff*_2 if it has the same characteristic as `*ff*_2 and its extension degree divides that of `*ff*_2

Finite field elements `*e*_1∈*ff*_1 and `*e*_2∈*ff*_2 define a field embedding of `*ff*_1 in `*ff*_2 iff they have the same [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html) and `*e*_1 generates `*ff*_1 The latter condition is satisfied iff the degree of the minimal polynomial of `*e*_1 is equal to the extension degree of `*ff*_1 over $\mathbb{Z}_{p}$.

For an embedding `ℰ=FiniteFieldEmbedding[*e*_1->*e*_2]`, `ℰ["Projection"]` represents a linear mapping $\mathcal{P}$ from the ambient field `*ff*_2 of `*e*_2 onto the ambient field `*ff*_1 of `*e*_1, treated as vector spaces over $\mathbb{Z}_{p}$, such that $\mathcal{P}(\mathcal{E}(a))=a$ for all $\mathit{a}\in \mathit{ff}_{1}$.

## Examples

### Basic Examples

Represent finite fields $\mathcal{K}$ and $\mathcal{F}$ with characteristic $17$ and extension degrees $3$ and $6$:

```wolfram
{𝒦,ℱ}={FiniteField[17,3],FiniteField[17,6]}
(* Output *)
{FiniteField[...],FiniteField[...]}
```

Find an embedding of $\mathcal{K}$ in $\mathcal{F}$:

```wolfram
ℰ=FiniteFieldEmbedding[𝒦,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[17, 14, +, #, +, #, ^, 3, &, Polynomial], 01], index -> 17, shortIndex -> 17, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[17, 3, +, 3, #, +, 10, #, ^, 2, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 146131612], index -> 18389764, shortIndex -> 18389764, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Map an element of $\mathcal{K}$ through the embedding:

```wolfram
ℰ[𝒦[123]]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[17, 3, +, 3, #, +, 10, #, ^, 2, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 08741016], index -> 23574733, shortIndex -> 23574733, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

Project the result back to $\mathcal{K}$:

```wolfram
ℰ["Projection"][%]
(* Output *)
![image](img/image_001.png)
```

### Scope

Represent finite fields $\mathcal{K}$ and $\mathcal{F}$ with characteristic $29$ and extension degrees $4$ and $12$:

```wolfram
{𝒦,ℱ}={FiniteField[29,4],FiniteField[29,12]}
(* Output *)
{FiniteField[...],FiniteField[...]}
```

Find an embedding of $\mathcal{K}$ in $\mathcal{F}$:

```wolfram
ℰ=FiniteFieldEmbedding[𝒦,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[29, 2, +, 15, #, +, 2, #, ^, 2, +, #, ^, 4, &, Polynomial], 01], index -> 29, shortIndex -> 29, indexShortened -> True, characteristic -> 29, shortCharacteristic -> 29, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[29, 2, +, #, +, #, ^, 2, +, 25, #, ^, 3, +, 16, #, ^, 4, +, 9, #, ^, 5, +, 28, #, ^, 6, +, 19, #, ^, 7, +, 3, #, ^, 8, +, #, ^, 12, &, Polynomial], 275243162209926226], index -> 318436523585826717, shortIndex -> "31843<<8>>26717", indexShortened -> True, characteristic -> 29, shortCharacteristic -> 29, extensionDegree -> 12, field -> FiniteField[...], fieldDisplayed -> False|>]
```

A field embedding preserves addition and multiplication:

```wolfram
{a,b}={𝒦[123],𝒦[456]}
(* Output *)
![image](img/image_003.png)
```

```wolfram
ℰ[a+b]==ℰ[a]+ℰ[b]
(* Output *)
True
```

```wolfram
ℰ[a b]==ℰ[a]ℰ[b]
(* Output *)
True
```

`ℰ["Projection"]` is a $\mathbb{Z}_{p}$-linear mapping but does not preserve multiplication:

```wolfram
{c,d}={ℱ[123],ℱ[456]}
(* Output *)
![image](img/image_005.png)
```

```wolfram
ℰ["Projection"][2c+3d]==2ℰ["Projection"][c]+3ℰ["Projection"][d]
(* Output *)
True
```

```wolfram
ℰ["Projection"][c d]==ℰ["Projection"][c]ℰ["Projection"][d]
(* Output *)
False
```

The composition of `ℰ["Projection"]` with $\mathcal{E}$ is the identity on $\mathcal{K}$:

```wolfram
ℰ["Projection"][ℰ[a]]==a
(* Output *)
True
```

The reverse composition is not the identity on $\mathcal{F}$:

```wolfram
ℰ[ℰ["Projection"][c]]==c
(* Output *)
False
```

Specify a field embedding by manually picking a generator and its value:

```wolfram
{𝒦,ℱ}={FiniteField[97,4],FiniteField[97,8]};
a=𝒦[123]
(* Output *)
![image](img/image_007.png)
```

`*a*` generates $\mathcal{K}$ if the degree of its minimal polynomial equals the extension degree of $\mathcal{K}$:

```wolfram
f=MinimalPolynomial[a,x]
(* Output *)
50+80 x+85 x^2+90 x^3+x^4
```

Find the roots of `*f*` in $\mathcal{F}$:

```wolfram
Factor[f,Extension->ℱ]
(* Output *)
![image](img/image_009.png)
```

Pick one of the roots:

```wolfram
b=-<|interpretation -> FiniteFieldElement[FiniteField[97, 5, +, 32, #, +, #, ^, 2, +, 65, #, ^, 3, +, #, ^, 8, &, Polynomial], 199024201195919], index -> 803801802386087, shortIndex -> 803801802386087, indexShortened -> True, characteristic -> 97, shortCharacteristic -> 97, extensionDegree -> 8, field -> FiniteField[...], fieldDisplayed -> False|>
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[97, 5, +, 32, #, +, #, ^, 2, +, 65, #, ^, 3, +, #, ^, 8, &, Polynomial], 7877377862688], index -> 7115271725265633, shortIndex -> "71152<<6>>65633", indexShortened -> True, characteristic -> 97, shortCharacteristic -> 97, extensionDegree -> 8, field -> FiniteField[...], fieldDisplayed -> False|>
```

Represent the embedding of $\mathcal{K}$ in $\mathcal{F}$ that maps `*a*` to `*b*`:

```wolfram
ℰ=FiniteFieldEmbedding[a->b]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[97, 5, +, 80, #, +, 6, #, ^, 2, +, #, ^, 4, &, Polynomial], 261], index -> 123, shortIndex -> 123, indexShortened -> True, characteristic -> 97, shortCharacteristic -> 97, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[97, 5, +, 32, #, +, #, ^, 2, +, 65, #, ^, 3, +, #, ^, 8, &, Polynomial], 7877377862688], index -> 7115271725265633, shortIndex -> "71152<<6>>65633", indexShortened -> True, characteristic -> 97, shortCharacteristic -> 97, extensionDegree -> 8, field -> FiniteField[...], fieldDisplayed -> False|>]
```

For the embedding to exist, both fields need to have the same characteristic:

```wolfram
FiniteFieldEmbedding[FiniteField[2,3],FiniteField[3,6]]
(* Output *)
FiniteFieldEmbedding
(* Output *)
FiniteFieldEmbedding[FiniteField[...],FiniteField[...]]
```

The extension degree of the first field needs to divide the extension degree of the second field:

```wolfram
FiniteFieldEmbedding[FiniteField[2,4],FiniteField[2,6]]
(* Output *)
FiniteFieldEmbedding
(* Output *)
FiniteFieldEmbedding[FiniteField[...],FiniteField[...]]
```

### Applications

Factor a polynomial in an algebraic extension of a finite field:

```wolfram
𝒦=FiniteField[59,3];
f=𝒦[123]x^2+𝒦[234]x+𝒦[345]
(* Output *)
![image](img/image_011.png)
```

```wolfram
IrreduciblePolynomialQ[f]
(* Output *)
True
```

Embed $\mathcal{K}$ in a finite field with $59^{6}$ elements:

```wolfram
ℱ=FiniteField[59,6];
ℰ=FiniteFieldEmbedding[𝒦,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[59, 57, +, 5, #, +, #, ^, 3, &, Polynomial], 01], index -> 59, shortIndex -> 59, indexShortened -> True, characteristic -> 59, shortCharacteristic -> 59, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[59, 2, +, 38, #, ^, 2, +, 18, #, ^, 3, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 30334435938], index -> 27283523017, shortIndex -> 27283523017, indexShortened -> True, characteristic -> 59, shortCharacteristic -> 59, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Map `*f*` through the embedding:

```wolfram
ff=f/.a_FiniteFieldElement:>ℰ[a]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[59, 2, +, 38, #, ^, 2, +, 18, #, ^, 3, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 234743574513], index -> 9851156214, shortIndex -> 9851156214, indexShortened -> True, characteristic -> 59, shortCharacteristic -> 59, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>+<|interpretation -> FiniteFieldElement[FiniteField[59, 2, +, 38, #, ^, 2, +, 18, #, ^, 3, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 294014462755], index -> 39657503749, shortIndex -> 39657503749, indexShortened -> True, characteristic -> 59, shortCharacteristic -> 59, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|> x+<|interpretation -> FiniteFieldElement[FiniteField[59, 2, +, 38, #, ^, 2, +, 18, #, ^, 3, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 6729111817], index -> 12374186118, shortIndex -> 12374186118, indexShortened -> True, characteristic -> 59, shortCharacteristic -> 59, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|> x^2
```

Factor the result:

```wolfram
Factor[ff]
(* Output *)
![image](img/image_013.png)
```

Use the [Extension](https://reference.wolfram.com/language/ref/Extension.html) option to combine the last two steps:

```wolfram
Factor[f,Extension->ℰ]
(* Output *)
![image](img/image_015.png)
```

### Properties & Relations

A field embedding preserves addition and multiplication:

```wolfram
{𝒦,ℱ}={FiniteField[7,3],FiniteField[7,9]};
ℰ=FiniteFieldEmbedding[𝒦,ℱ];
```

```wolfram
{a,b}={𝒦[123],𝒦[456]}
(* Output *)
![image](img/image_017.png)
```

```wolfram
ℰ[a+b]==ℰ[a]+ℰ[b]
(* Output *)
True
```

```wolfram
ℰ[a b]==ℰ[a]ℰ[b]
(* Output *)
True
```

`ℰ["Projection"]` is a $\mathbb{Z}_{p}$-linear mapping but does not preserve multiplication:

```wolfram
{c,d}={ℱ[123],ℱ[456]}
(* Output *)
![image](img/image_019.png)
```

```wolfram
ℰ["Projection"][2c+3d]==2ℰ["Projection"][c]+3ℰ["Projection"][d]
(* Output *)
True
```

```wolfram
ℰ["Projection"][c d]==ℰ["Projection"][c]ℰ["Projection"][d]
(* Output *)
False
```

The composition of `ℰ["Projection"]` with $\mathcal{E}$ is the identity on $\mathcal{K}$:

```wolfram
ℰ["Projection"][ℰ[a]]==a
(* Output *)
True
```

The reverse composition is not the identity on $\mathcal{F}$:

```wolfram
ℰ[ℰ["Projection"][c]]==c
(* Output *)
False
```

Find an automorphism of $\mathcal{F}$:

```wolfram
ℱ=FiniteField[71,5];
aut=FiniteFieldEmbedding[ℱ,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 18, #, +, #, ^, 5, &, Polynomial], 01], index -> 71, shortIndex -> 71, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 5, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 18, #, +, #, ^, 5, &, Polynomial], 5313561352], index -> 1326343527, shortIndex -> 1326343527, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 5, field -> FiniteField[...], fieldDisplayed -> False|>]
```

All finite field automorphisms are functional powers of the Frobenius automorphism:

```wolfram
Table[FrobeniusAutomorphism[ℱ[71],k],{k,4}]
(* Output *)
![image](img/image_021.png)
```

Here `*aut*[*a*]==[FrobeniusAutomorphism](https://reference.wolfram.com/language/ref/FrobeniusAutomorphism.html)[*a*,4]`:

```wolfram
aut[ℱ[1234]]==FrobeniusAutomorphism[ℱ[1234],4]
(* Output *)
True
```

An embedding allows identifying $\mathcal{K}$ with a subfield of $\mathcal{F}$:

```wolfram
{𝒦,ℱ}={FiniteField[43,3],FiniteField[43,6]};
ℰ=FiniteFieldEmbedding[𝒦,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[43, 40, +, #, +, #, ^, 3, &, Polynomial], 01], index -> 43, shortIndex -> 43, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 36102317221], index -> 3095409517, shortIndex -> 3095409517, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Use [FiniteFieldElementTrace](https://reference.wolfram.com/language/ref/FiniteFieldElementTrace.html) to compute $Tr_{\mathcal{F}/\mathcal{K}}(a)$:

```wolfram
FiniteFieldElementTrace[ℱ[123],ℰ]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 40, +, #, +, #, ^, 3, &, Polynomial], 22168], index -> 15502, shortIndex -> 15502, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

Use [FiniteFieldElementNorm](https://reference.wolfram.com/language/ref/FiniteFieldElementNorm.html) to compute $N_{\mathcal{F}/\mathcal{K}}(a)$:

```wolfram
FiniteFieldElementNorm[ℱ[123],ℰ]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 40, +, #, +, #, ^, 3, &, Polynomial], 201019], index -> 35581, shortIndex -> 35581, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

Use [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html) to find the minimal polynomial of an element of $\mathcal{F}$ over $\mathcal{K}$:

```wolfram
MinimalPolynomial[ℱ[123],x,ℰ]
(* Output *)
![image](img/image_023.png)
```

Use [Composition](https://reference.wolfram.com/language/ref/Composition.html) to compose finite field embeddings:

```wolfram
𝒦=FiniteField[71,3];
ℱ=FiniteField[71,6];
𝒢=FiniteField[71,12];
ℰ_1=FiniteFieldEmbedding[𝒦,ℱ];
ℰ_2=FiniteFieldEmbedding[ℱ,𝒢];
ℰ_21=ℰ_2@*ℰ_1
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[71, 64, +, 4, #, +, #, ^, 3, &, Polynomial], 01], index -> 71, shortIndex -> 71, indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[71, 7, +, 23, #, +, 58, #, ^, 2, +, 21, #, ^, 3, +, 55, #, ^, 4, +, 29, #, ^, 5, +, 28, #, ^, 6, +, 12, #, ^, 7, +, #, ^, 12, &, Polynomial], 66377049769336955615433], index -> 7805651698674017045103, shortIndex -> "78056<<12>>45103", indexShortened -> True, characteristic -> 71, shortCharacteristic -> 71, extensionDegree -> 12, field -> FiniteField[...], fieldDisplayed -> False|>]
```

```wolfram
ℰ_21[𝒦[123]]==ℰ_2[ℰ_1[𝒦[123]]]
(* Output *)
True
```

## Related Guides ▪Finite Fields ▪Polynomial Algebra

## History Introduced in 2023 (13.3)
