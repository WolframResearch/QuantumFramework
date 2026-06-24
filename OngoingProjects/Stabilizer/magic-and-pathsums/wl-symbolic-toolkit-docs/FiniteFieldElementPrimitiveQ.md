# FiniteFieldElementPrimitiveQ | [SpanFromLeft]

> [FiniteFieldElementPrimitiveQ](https://reference.wolfram.com/language/ref/FiniteFieldElementPrimitiveQ.html)[*a*]  — tests whether `*a*` is a primitive element of its ambient field.

## Details

`*a*` is a primitive element of a finite field $\mathcal{F}$ if the minimal polynomial of `*a*` is primitive and has degree equal to the extension degree of $\mathcal{F}$ over $\mathbb{Z}_{p}$.

A finite field with `*q*` elements contains [EulerPhi](https://reference.wolfram.com/language/ref/EulerPhi.html)[*q*-1] primitive elements.

## Examples

### Basic Examples

Test whether a finite field element is a primitive element of its ambient field:

```wolfram
ℱ=FiniteField[17,3];
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[123]]
(* Output *)
True
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[125]]
(* Output *)
False
```

### Scope

Test elements of a finite field in polynomial representation:

```wolfram
ℱ=FiniteField[179,2];
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[123]]
(* Output *)
False
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[234]]
(* Output *)
True
```

Test elements of a finite field in exponential representation:

```wolfram
ℱ=FiniteField[2,10,"Exponential"];
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[2]]
(* Output *)
True
```

```wolfram
FiniteFieldElementPrimitiveQ[ℱ[4]]
(* Output *)
False
```

Arguments that are not [FiniteFieldElement](https://reference.wolfram.com/language/ref/FiniteFieldElement.html) objects yield [False](https://reference.wolfram.com/language/ref/False.html):

```wolfram
FiniteFieldElementPrimitiveQ[17]
(* Output *)
False
```

```wolfram
FiniteFieldElementPrimitiveQ[x]
(* Output *)
False
```

### Applications

Implement a discrete logarithm function for a given base element:

```wolfram
ℱ=FiniteField[3,10];
b=ℱ[123]
(* Output *)
![image](img/image_001.png)
```

The base element needs to be primitive:

```wolfram
FiniteFieldElementPrimitiveQ[b]
(* Output *)
True
```

Construct a field with exponential element representation using the minimal polynomial of `*b*`:

```wolfram
𝒢=FiniteField[3,MinimalPolynomial[b],"Exponential"]
(* Output *)
FiniteField[...]
```

Construct a field isomorphism that maps `*b*` to the field generator of `𝒢`:

```wolfram
ℐ=FiniteFieldEmbedding[b->𝒢[2]]
(* Output *)
FiniteFieldEmbedding[<|interpretation -> FiniteFieldElement[FiniteField[3, 2, +, #, +, 2, #, ^, 4, +, 2, #, ^, 5, +, 2, #, ^, 6, +, #, ^, 10, &, Polynomial], 02111], index -> 123, shortIndex -> 123, indexShortened -> True, characteristic -> 3, shortCharacteristic -> 3, extensionDegree -> 10, field -> FiniteField[...], fieldDisplayed -> False|>-><|interpretation -> FiniteFieldElement[FiniteField[3, 2, +, #, +, 2, #, ^, 4, +, 2, #, ^, 8, +, 2, #, ^, 9, +, #, ^, 10, &, Exponential], 2], index -> 2, shortIndex -> 2, indexShortened -> True, characteristic -> 3, shortCharacteristic -> 3, extensionDegree -> 10, field -> FiniteField[...], fieldDisplayed -> False|>]
```

Define the discrete logarithm function:

```wolfram
dlog[a_]:=ℐ[a]["Index"]-1
```

Compute the discrete logarithm for a nonzero element of `ℱ`:

```wolfram
dlog[ℱ[456]]
(* Output *)
15981
```

```wolfram
b^%
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[3, 2, +, #, +, 2, #, ^, 4, +, 2, #, ^, 5, +, 2, #, ^, 6, +, #, ^, 10, &, Polynomial], 022121], index -> 456, shortIndex -> 456, indexShortened -> True, characteristic -> 3, shortCharacteristic -> 3, extensionDegree -> 10, field -> FiniteField[...], fieldDisplayed -> False|>
```

### Properties & Relations

`*a*` is a primitive element of its ambient field if every nonzero element of the field is an integer power of `*a*`:

```wolfram
ℱ=FiniteField[5,3];
a=ℱ[123];
FiniteFieldElementPrimitiveQ[a]
(* Output *)
True
```

```wolfram
Information[ℱ,"FieldSize"]
(* Output *)
125
```

```wolfram
Length[Union[a^Range[124]]]
(* Output *)
124
```

`*a*` is a primitive element of $\mathcal{F}$ if it is a generator of $\mathcal{F}$ and its minimal polynomial is primitive:

```wolfram
ℱ=FiniteField[59,6];
a=ℱ[123456];
FiniteFieldElementPrimitiveQ[a]
(* Output *)
True
```

Use [MinimalPolynomial](https://reference.wolfram.com/language/ref/MinimalPolynomial.html) to find the minimal polynomial of `*a*`:

```wolfram
f=MinimalPolynomial[a,x]
(* Output *)
31+55 x+57 x^2+3 x^3+45 x^4+31 x^5+x^6
```

This shows that `*a*` is a generator of $\mathcal{F}$:

```wolfram
Exponent[f,x]==Information[ℱ,"ExtensionDegree"]
(* Output *)
True
```

Use [PrimitivePolynomialQ](https://reference.wolfram.com/language/ref/PrimitivePolynomialQ.html) to show that `*f*` is primitive:

```wolfram
PrimitivePolynomialQ[f,59]
(* Output *)
True
```

The multiplicative order of a primitive element is one less than its ambient field size:

```wolfram
ℱ=FiniteField[1009,2];
a=ℱ[1017];
FiniteFieldElementPrimitiveQ[a]
(* Output *)
True
```

```wolfram
Information[ℱ,"FieldSize"]
(* Output *)
1018081
```

Use [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html) to compute the multiplicative order of `*a*`:

```wolfram
MultiplicativeOrder[a]
(* Output *)
1018080
```

A finite field with `*q*` elements contains [EulerPhi](https://reference.wolfram.com/language/ref/EulerPhi.html)[*q*-1] primitive elements:

```wolfram
ℱ=FiniteField[11,5];
Count[FiniteFieldElementPrimitiveQ[ℱ[#]]&/@Range[0,Information[ℱ,"FieldSize"]-1],True]
(* Output *)
64400
```

```wolfram
EulerPhi[Information[ℱ,"FieldSize"]-1]
(* Output *)
64400
```

## Related Guides ▪Finite Fields

## History Introduced in 2023 (13.3)
