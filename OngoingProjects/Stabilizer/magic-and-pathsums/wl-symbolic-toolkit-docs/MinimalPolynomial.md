# MinimalPolynomial

> [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*s*,*x*]  — gives the minimal polynomial in `*x*` for which the algebraic number `*s*` is a root.
> [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*,*x*]  — gives the minimal polynomial of the finite field element `*u*` over $\mathbb{Z}_{p}$.
> [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*,*x*,*k*]  — gives the minimal polynomial of `*u*` over the $p^{k}$-element subfield of the ambient field of `*u*`.
> [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*,*x*,*emb*]  — gives the minimal polynomial of `*u*` relative to the finite field embedding `*emb*`.

## Details and Options

[MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*s*,*x*] gives the lowest-degree polynomial with integer coefficients, positive leading coefficient and the [GCD](https://reference.wolfram.com/language/ref/GCD.html) of all coefficients equal to $1$ for which the algebraic number `*s*` is a root.

[MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*s*] gives a pure function representation of the minimal polynomial of `*s*`.

[MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*s*,*x*,[Extension](https://reference.wolfram.com/language/ref/Extension.html)->*a*] finds the characteristic polynomial of $s \in \mathbb{Q}[a]$ over the field $\mathbb{Q}[a]$.

For a [FiniteFieldElement](https://reference.wolfram.com/language/ref/FiniteFieldElement.html) object `*u*` in a finite field $\mathcal{F}$ of characteristic $p$, [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*, *x*] gives the lowest-degree monic polynomial with integer coefficients between $0$ and $p-1$ for which `*u*` is a root.

[MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*,*x*,*k*] gives the lowest-degree monic polynomial with coefficients from the $p^{k}$-element subfield of $\mathcal{F}$ for which `*u*` is a root. `*k*` needs to be a divisor of the extension degree of $\mathcal{F}$ over $\mathbb{Z}_{p}$.

If `*emb*=[FiniteFieldEmbedding](https://reference.wolfram.com/language/ref/FiniteFieldEmbedding.html)[*e*_1->*e*_2]`, then [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html)[*u*,*x*,*emb*] gives the polynomial with coefficients in the ambient field of `*e*_1 that map through `*emb*` to the coefficients of the minimal polynomial of `*u*` over the image of `*emb*`.

## Examples

### Basic Examples

Minimal polynomials of algebraic numbers:

```wolfram
MinimalPolynomial[Sqrt[3],x]
(* Output *)
-3+x^2
```

```wolfram
MinimalPolynomial[Sqrt[2]+Sqrt[3],x]
(* Output *)
1-10 x^2+x^4
```

Minimal polynomials of finite field elements:

```wolfram
ℱ=FiniteField[2,12];
```

```wolfram
MinimalPolynomial[ℱ[123],x]
(* Output *)
1+x+x^3+x^4+x^5+x^7+x^10+x^11+x^12
```

```wolfram
MinimalPolynomial[ℱ[123],x,3]
(* Output *)
![image](img/image_001.png)
```

### Scope

#### Algebraic Numbers

Radical expressions:

```wolfram
MinimalPolynomial[Sqrt[1+Sqrt[3]],x]
(* Output *)
-2-2 x^2+x^4
```

```wolfram
MinimalPolynomial[(1+I)/Sqrt[2],x]
(* Output *)
1+x^4
```

[Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
MinimalPolynomial[Root[2#1^3-2 #1+7&,1]+17,x]
(* Output *)
-9785+1732 x-102 x^2+2 x^3
```

[AlgebraicNumber](https://reference.wolfram.com/language/ref/AlgebraicNumber.html) objects:

```wolfram
MinimalPolynomial[AlgebraicNumber[Root[-2+#1^3&,2],{2,2,1}],x]
(* Output *)
-4-6 x^2+x^3
```

[MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html) automatically threads over lists:

```wolfram
MinimalPolynomial[{Sqrt[3],Root[-2+#1^3&,2]+1},x]
(* Output *)
{-3+x^2,-3+3 x-3 x^2+x^3}
```

Pure function minimal polynomial:

```wolfram
MinimalPolynomial[Sqrt[2]+Sqrt[3]]
(* Output *)
1-10 #1^2+#1^4&
```

#### Finite Field Elements

Represent a finite field with characteristic $17$ and extension degree $6$:

```wolfram
ℱ=FiniteField[17,6]
(* Output *)
FiniteField[...]
```

Minimal polynomial over $\mathbb{Z}_{p}$:

```wolfram
MinimalPolynomial[ℱ[123],x]
(* Output *)
3+16 x+11 x^2+8 x^3+15 x^4+10 x^5+x^6
```

Minimal polynomial over $\mathbb{Z}_{p}$ with coefficients given as elements of $\mathcal{F}$:

```wolfram
MinimalPolynomial[ℱ[123],x,1]
(* Output *)
![image](img/image_003.png)
```

Minimal polynomial over the $17^{2}$-element subfield of $\mathcal{F}$:

```wolfram
MinimalPolynomial[ℱ[123],x,2]
(* Output *)
![image](img/image_005.png)
```

Embed a field $\mathcal{K}$ with $17^{3}$ elements in $\mathcal{F}$:

```wolfram
𝒦=FiniteField[17,3];
ℰ=FiniteFieldEmbedding[𝒦,ℱ]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[17, 14, +, #, +, #, ^, 3, &, Polynomial], 01], index -> 17, shortIndex -> 17, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[17, 3, +, 3, #, +, 10, #, ^, 2, +, 2, #, ^, 4, +, #, ^, 6, &, Polynomial], 146131612], index -> 18389764, shortIndex -> 18389764, indexShortened -> True, characteristic -> 17, shortCharacteristic -> 17, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Minimal polynomial relative to the finite field embedding $\mathcal{E}$:

```wolfram
MinimalPolynomial[ℱ[123],x,ℰ]
(* Output *)
![image](img/image_007.png)
```

Pure function minimal polynomial:

```wolfram
MinimalPolynomial[ℱ[123]]
(* Output *)
3+16 #1+11 #1^2+8 #1^3+15 #1^4+10 #1^5+#1^6&
```

### Options

#### Extension

Find the characteristic polynomial of $\sqrt{2}$ over the extension $\mathbb{Q}[e^{i \pi/4}]$ of $\mathbb{Q}$:

```wolfram
MinimalPolynomial[Sqrt[2],x,Extension -> E^( I Pi/4)]
(* Output *)
4-4 x^2+x^4
```

The characteristic polynomial is a power of the minimal polynomial of $\sqrt{2}$:

```wolfram
{Factor[%],MinimalPolynomial[Sqrt[2],x]}
(* Output *)
{(-2+x^2)^2,-2+x^2}
```

### Applications

Construct a polynomial with a root $\sqrt{1+\sqrt{2}}$:

```wolfram
MinimalPolynomial[Sqrt[1+Sqrt[2]],x]
(* Output *)
-1-2 x^2+x^4
```

```wolfram
Plot[%,{x,-2,2},Epilog->{Blue,PointSize[Large],Point[{Sqrt[1+Sqrt[2]],0}]}]
```

*([Graphics])*

The degree of the number field generated by `(2-[I](https://reference.wolfram.com/language/ref/I.html))/[Sqrt](https://reference.wolfram.com/language/ref/Sqrt.html)[5]`:

```wolfram
Exponent[MinimalPolynomial[(2-I)/Sqrt[5],x],x]
(* Output *)
4
```

Check whether a finite field element generates its ambient field:

```wolfram
ℱ=FiniteField[3,10];
```

```wolfram
Exponent[MinimalPolynomial[ℱ[123],x],x]==Information[ℱ,"ExtensionDegree"]
(* Output *)
True
```

```wolfram
Exponent[MinimalPolynomial[ℱ[636],x],x]==Information[ℱ,"ExtensionDegree"]
(* Output *)
False
```

### Properties & Relations

Compute the extension that defines the number field $F=\mathbb{Q}[\sqrt{3},e^{i \pi/4}]$:

```wolfram
F=ToNumberField[{Sqrt[3], E^(I Pi/4)}, All][[1,1]]
(* Output *)
Root
```

Find the characteristic polynomial of $\sqrt{2}+\sqrt{3}$ over $F$:

```wolfram
MinimalPolynomial[Sqrt[2]+Sqrt[3], x,Extension -> F]
(* Output *)
1-20 x^2+102 x^4-20 x^6+x^8
```

The characteristic polynomial is a power of the minimal polynomial of $\sqrt{2}+\sqrt{3}$:

```wolfram
{Factor[%],MinimalPolynomial[Sqrt[2]+Sqrt[3],x]}
(* Output *)
{(1-10 x^2+x^4)^2,1-10 x^2+x^4}
```

Use [FrobeniusAutomorphism](https://reference.wolfram.com/language/ref/FrobeniusAutomorphism.html) to find all conjugates of a finite field element `*a*`:

```wolfram
ℱ=FiniteField[7,5];
a=ℱ[123];
conj=Table[FrobeniusAutomorphism[a,k],{k,0,4}]
(* Output *)
![image](img/image_009.png)
```

The conjugates are roots of the minimal polynomial of `*a*`:

```wolfram
MinimalPolynomial[a]/@conj
(* Output *)
![image](img/image_011.png)
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
![image](img/image_013.png)
```

```wolfram
-Coefficient[%,x,2]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 2229363714], index -> 2185097697, shortIndex -> 2185097697, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
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
![image](img/image_015.png)
```

```wolfram
-Coefficient[%,x,0]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[43, 3, +, 21, #, +, 28, #, ^, 2, +, 19, #, ^, 3, +, #, ^, 6, &, Polynomial], 37261938517], index -> 2519295088, shortIndex -> 2519295088, indexShortened -> True, characteristic -> 43, shortCharacteristic -> 43, extensionDegree -> 6, field -> FiniteField[...], fieldDisplayed -> False|>
```

## Tech Notes ▪Algebraic Number Fields

## Related Guides ▪Algebraic Numbers ▪Polynomial Algebra ▪Number Theory ▪Finite Fields ▪Algebraic Number Theory ▪Number Recognition

## History Introduced in 2007 (6.0) | Updated in 2023 (13.3) ▪ 2026 (15.0)
